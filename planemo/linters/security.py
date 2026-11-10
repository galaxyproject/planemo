"""Opt-in security checks for Galaxy tool wrappers.

These checks look for common patterns where user-controlled values reach a
tool's ``<command>`` template in a way that can carry shell metacharacters,
for instance a collection's ``element_identifier`` interpolated directly into
the command. The checks are informational: they only run when ``--security``
is passed and INFO findings never fail ``planemo lint``.

A finding does not prove a tool is exploitable. Cheetah logic is not evaluated
statically, so quoting or cleaning elsewhere may make an occurrence safe and
every finding should be reviewed. For raw collection identifiers Galaxy 26.2
offers the ``safe_element_identifier`` accessor (`Galaxy PR #23735
<https://github.com/galaxyproject/galaxy/pull/23735>`_); on older servers
cleaning a value with ``re.sub(r'[^\\w.-]', '_', ...)`` before interpolating
it is the recommended fallback.

Each check can be skipped individually (e.g. ``--skip SecurityRawIdentifier``)
and the whole module can be skipped with ``--skip security``.
"""

import re
import string
from typing import (
    Any,
    Callable,
    Dict,
    List,
    NamedTuple,
    Optional,
    Set,
    Tuple,
    TYPE_CHECKING,
)

from galaxy.tool_util.lint import Linter
from galaxy.util import mapped_chars

from planemo.linters.util import cached_lint_result

if TYPE_CHECKING:
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource

lint_tool_types = ["*"]

# Shell metacharacters that are dangerous when they survive sanitization and
# reach a <command>. Galaxy's default sanitization removes or rewrites all of
# these, so they only survive when sanitization is disabled or reconfigured.
SHELL_DANGEROUS_CHARS = "'\"$;|&`\\<>#"

# Characters that are safe to keep in a value interpolated into a shell
# command without further quoting. This mirrors the character set accepted by
# Galaxy's ``safe_filename_component`` helper.
SAFE_FILENAME_CHARS = frozenset(string.ascii_letters + string.digits + "-_.")

DEFAULT_VALID_CHARS = string.ascii_letters + string.digits + " -=_.()/+*^,:?!"

STRING_PRESET_ALIASES = {
    "letters": "ascii_letters",
    "lowercase": "ascii_lowercase",
    "uppercase": "ascii_uppercase",
}

MAPPING_PRESET: Dict[str, Dict[str, str]] = {
    "default": dict(mapped_chars),
    "none": {},
}

_TRUE_VALUES = frozenset(("true", "t", "yes", "y", "on", "1"))
_FALSE_VALUES = frozenset(("false", "f", "no", "n", "off", "0"))

# Cheetah directives that evaluate a condition without rendering it, so raw
# values on such lines do not flow into the rendered command.
_CONDITION_DIRECTIVES = frozenset(("if", "elif", "while", "unless"))

# Directive names Cheetah recognizes (Cheetah's Parser.directiveNamesAndParsers).
# A ``#`` followed by any other word is passed through literally and can then
# start a bash comment in the rendered command.
_CHEETAH_DIRECTIVE_NAMES = (
    "compiler-settings",
    "errorCatcher",
    "defmacro",
    "implements",
    "breakpoint",
    "assert",
    "capture",
    "closure",
    "extends",
    "filter",
    "import",
    "include",
    "silent",
    "return",
    "unless",
    "repeat",
    "while",
    "yield",
    "block",
    "break",
    "cache",
    "call",
    "else",
    "pass",
    "raise",
    "slurp",
    "stop",
    "attr",
    "echo",
    "from",
    "raw",
    "set",
    "try",
    "call",
    "def",
    "del",
    "end",
    "for",
    "if",
    "arg",
)
_DIRECTIVE = re.compile(r"#(" + "|".join(_CHEETAH_DIRECTIVE_NAMES) + r")(?![A-Za-z0-9_@-])")
_CLEANING_CALL = re.compile(r"\b(?:re\.sub|shlex\.quote)\s*\(")
_RAW_IDENTIFIER = re.compile(r"\$(?:\(|\{)?\s*\$?\s*([A-Za-z_][\w.]*)\s*\.\s*element_identifier\b")
_COLLECTION_NAME = re.compile(r"\$(?:\(|\{)?\s*\$?\s*([A-Za-z_][\w.]*)\s*\.\s*name\b")
_COLLECTION_KEYS = re.compile(r"\$(?:\(|\{)?\s*\$?\s*([A-Za-z_][\w.]*)\s*\.\s*keys\s*\(\s*\)?")
# Cheetah does not require the ``$`` on a loop variable, so both forms bind.
_FOR_SIMPLE = re.compile(r"#\s*for\s+\$?([A-Za-z_]\w*)\s+in\s+\$([A-Za-z_][\w.]*)")
_FOR_ENUMERATE = re.compile(
    r"#\s*for\s+\$?([A-Za-z_]\w*)\s*,\s*\$?([A-Za-z_]\w*)\s+in\s+enumerate\(\s*\$([A-Za-z_][\w.]*)"
)
# ``#for`` directives. Cheetah comments (``##`` lines and ``#*`` blocks) are
# skipped separately, and a ``#`` followed by whitespace is not a directive at
# all (Cheetah renders it literally). ``#end for`` is not matched: a loop
# variable outlives its loop, so where a loop ends does not matter.
_FOR_TOKEN = re.compile(r"(?<!#)#for\b")
_FOR_TARGETS = re.compile(r"#\s*for\s+(.+?)\s+in\b")
_LOOP_VAR = re.compile(r"\$?([A-Za-z_]\w*)")

# Types whose values are free text passed through parameter sanitization.
# ``hidden`` inputs are text carried in the tool state, so they go through the
# same sanitization.
_SANITIZED_PARAM_TYPES = frozenset(("text", "select", "hidden"))

# Types whose values are never sanitized by Galaxy when interpolated.
_RAW_VALUE_PARAM_TYPES = frozenset(("data_column", "group_tag"))

Emit = Callable[[str, str, str, Any], None]


class SecurityCheck(NamedTuple):
    """A single finding, consumed by the ``Linter`` classes below."""

    linter: str
    message: str
    node: Any


class _ParamInfo(NamedTuple):
    """Type and XML element of an input parameter, keyed by dotted name."""

    type: str
    node: Any


class _LoopScope(NamedTuple):
    """Loop variables in effect from the start of one ``#for`` body on.

    ``bindings`` also holds the variables of every loop whose body started
    earlier, so the latest scope started at a position gives every binding in
    effect there (see ``_loop_bindings``).
    """

    start: int
    bindings: Dict[str, str]


def _parse_name_from_argument(argument: Optional[str]) -> Optional[str]:
    """Derive a param name from its ``argument`` flag, as Galaxy does."""
    if not argument:
        return None
    return argument.lstrip("-").replace("-", "_")


def _index_params(root: Any) -> Dict[str, _ParamInfo]:
    """Map dotted parameter names (``repeat.param``, ``cond.param``) to types."""
    params: Dict[str, _ParamInfo] = {}
    inputs = root.find("inputs")
    if inputs is None:
        return params

    def walk(element: Any, prefix: str) -> None:
        for child in element:
            if child.tag == "param":
                # Galaxy derives the name from `argument` (dashes stripped,
                # inner dashes replaced by underscores) when `name` is absent.
                name = child.get("name") or _parse_name_from_argument(child.get("argument"))
                if not name:
                    continue
                dotted = f"{prefix}.{name}" if prefix else name
                params[dotted] = _ParamInfo(child.get("type") or "", child)
            elif child.tag in ("repeat", "conditional", "section"):
                name = child.get("name") or ""
                sub_prefix = f"{prefix}.{name}" if prefix else name
                walk(child, sub_prefix)
                if child.tag == "conditional":
                    for when in child.findall("when"):
                        walk(when, sub_prefix)

    walk(inputs, "")
    return params


def _resolve_loop_expr(expr: str, loops: Dict[str, str]) -> str:
    """Rewrite a leading loop variable of ``expr`` into the expression it iterates."""
    expr = expr.strip()
    head, _, rest = expr.partition(".")
    if head not in loops:
        return expr
    target = loops[head]
    return f"{target}.{rest}" if rest and target else target


def _for_header(text: str, pos: int, header_end: int, earlier: Dict[str, str]) -> Tuple[Dict[str, str], int]:
    """Bindings declared by the ``#for`` header at ``pos`` and where its body starts.

    ``earlier`` are the loop bindings in effect at the header. The result
    extends them, and iterated expressions are resolved against them, so
    bindings always refer to parameters. Variables of loop forms that are not
    recognized are bound to an unknown ("") expression, which still shadows
    an earlier binding of the same name.
    """
    bindings = dict(earlier)
    targets = _FOR_TARGETS.match(text, pos)
    if targets:
        bindings.update((name, "") for name in _LOOP_VAR.findall(targets.group(1)))
        header_end = targets.end()
    enumerate_match = _FOR_ENUMERATE.match(text, pos)
    simple_match = _FOR_SIMPLE.match(text, pos)
    if enumerate_match:
        bindings[enumerate_match.group(1)] = "__index__"
        bindings[enumerate_match.group(2)] = _resolve_loop_expr(enumerate_match.group(3), earlier)
        header_end = enumerate_match.end()
    elif simple_match:
        bindings[simple_match.group(1)] = _resolve_loop_expr(simple_match.group(2), earlier)
        header_end = simple_match.end()
    return bindings, header_end


def _cheetah_comment_spans(text: str) -> List[Tuple[int, int]]:
    """Spans of Cheetah comments: ``##`` to the end of its line, ``#*`` to ``*#``."""
    spans: List[Tuple[int, int]] = []
    pos = 0
    size = len(text)
    while pos < size:
        if text.startswith("##", pos):
            end = _line_end(text, pos)
            spans.append((pos, end))
            pos = end
            continue
        if text.startswith("#*", pos):
            block_end = text.find("*#", pos + 2)
            end = size if block_end == -1 else block_end + 2
            spans.append((pos, end))
            pos = end
            continue
        pos += 1
    return spans


def _loop_bindings(text: str) -> List[_LoopScope]:
    """Scopes of Cheetah loop variables, ordered by the start of their body.

    Cheetah compiles ``#for`` into a Python ``for`` of the template function,
    so a loop variable keeps its last value after ``#end for``, and nesting
    does not matter: a variable is bound from the start of a body that binds
    it until the start of a later body that binds it again. Before the first
    body no loop variable is bound. Cheetah comments are not part of the
    template, so a ``#for`` inside one binds nothing.
    """
    comments = _cheetah_comment_spans(text)
    scopes: List[_LoopScope] = []
    for token in _FOR_TOKEN.finditer(text):
        if _in_spans(token.start(), comments):
            continue
        earlier = scopes[-1].bindings if scopes else {}
        bindings, body_start = _for_header(text, token.start(), token.end(), earlier)
        scopes.append(_LoopScope(body_start, bindings))
    return scopes


def _scope_at(scopes: List[_LoopScope], pos: int) -> Optional[_LoopScope]:
    """Loop scope started last at or before ``pos``, or None before the first loop body."""
    latest = None
    for scope in scopes:
        if scope.start > pos:
            break
        latest = scope
    return latest


def _loops_at(scopes: List[_LoopScope], pos: int) -> Dict[str, str]:
    """Loop bindings in effect at ``pos``."""
    scope = _scope_at(scopes, pos)
    return scope.bindings if scope else {}


def _expr_type(expr: str, params: Dict[str, _ParamInfo], loops: Dict[str, str]) -> Optional[str]:
    """Best-effort static type of a Cheetah expression, or None if unknown.

    ``loops`` are the loop bindings in effect where the expression appears
    (see ``_loops_at``).
    """
    param = params.get(_resolve_loop_expr(expr, loops))
    return param.type if param else None


def _mark(raw_flow: List[bool], unquoted: List[bool], start: int, end: int, value: bool) -> None:
    """Set both classification lists to ``value`` for positions in [start, end)."""
    for pos in range(max(0, start), min(len(raw_flow), end)):
        raw_flow[pos] = value
        unquoted[pos] = value


def _apply_hash(text: str, pos: int, raw_flow: List[bool], unquoted: List[bool]) -> int:
    """Handle a ``#`` at a word boundary at position ``pos`` of a command.

    Returns the position to continue scanning from. Comments and condition
    directives only affect their own line: the rendered command is a shell
    script in which later lines still run.
    """
    line_end = _line_end(text, pos)
    directive = _DIRECTIVE.match(text, pos)
    if directive is not None:
        if directive.group(0)[1:] in _CONDITION_DIRECTIVES:
            _mark(raw_flow, unquoted, pos, line_end, False)
        else:
            # Other Cheetah directives (``#set``, ``#for``, ``#echo``...) are
            # not rendered, but can flow raw values into the command.
            for offset in range(pos, line_end):
                unquoted[offset] = False
        return line_end
    # A bash comment makes the rest of its line inert.
    _mark(raw_flow, unquoted, pos, line_end, False)
    return line_end


def _line_end(text: str, pos: int) -> int:
    """Position just past the end of the line containing ``pos``."""
    end = text.find("\n", pos)
    return len(text) if end == -1 else end


def _quoted_positions(text: str) -> List[bool]:
    """True for each position inside a single- or double-quoted span.

    Cheetah removes ``##`` lines and ``#*`` blocks wherever they appear and
    processes its directives everywhere too, so those spans are skipped
    without touching the quote state: a quote opened before them stays open,
    because Cheetah may have eaten its closing quote. A bash comment (a ``#``
    at a word boundary that starts no Cheetah directive) is only recognized
    outside quotes: inside quotes it is literal text.
    """
    size = len(text)
    quoted = [False] * size
    in_single = in_double = False
    pos = 0
    while pos < size:
        if text.startswith("##", pos):
            pos = _line_end(text, pos)
            continue
        if text.startswith("#*", pos):
            block_end = text.find("*#", pos + 2)
            pos = size if block_end == -1 else block_end + 2
            continue
        char = text[pos]
        if in_single or in_double:
            quoted[pos] = True
            close = "'" if in_single else '"'
            if char == close:
                in_single = in_double = False
            elif in_double and char == "\\":
                pos = min(pos + 2, size)
                continue
            elif char == "#" and (pos == 0 or text[pos - 1] in " \t\n;|&(") and _DIRECTIVE.match(text, pos) is not None:
                # A Cheetah directive is processed even inside quotes, so its
                # content never reaches the shell and cannot close the quote.
                pos = _line_end(text, pos)
                continue
            pos += 1
            continue
        if char in "'\"":
            in_single = char == "'"
            in_double = char == '"'
        elif char == "\\":
            pos = min(pos + 2, size)
            continue
        elif char == "#" and (pos == 0 or text[pos - 1] in " \t\n;|&("):
            # A directive or bash comment: its content never reaches the shell.
            pos = _line_end(text, pos)
            continue
        pos += 1
    return quoted


def _scan_command(
    text: str, quoted: List[bool], cleaning_spans: List[Tuple[int, int]]
) -> Tuple[List[bool], List[bool]]:
    """Classify each position of a ``<command>`` template.

    Returns two boolean lists: ``raw_flow`` (a value interpolated at this
    position reaches the shell, quoted or not) and ``unquoted`` (the position
    is also outside quotes, Cheetah directives, and shell comments).
    """
    size = len(text)
    raw_flow: List[bool] = [True] * size
    unquoted: List[bool] = [False] * size
    span_ends = {start: end for start, end in cleaning_spans}
    for start, end in cleaning_spans:
        _mark(raw_flow, unquoted, start, end, False)
    pos = 0
    while pos < size:
        if pos in span_ends:
            pos = span_ends.pop(pos)
            continue
        if text.startswith("##", pos):
            # A Cheetah single-line comment renders nothing on its line, even
            # inside quotes (Cheetah does not know about shell quotes).
            line_end = _line_end(text, pos)
            _mark(raw_flow, unquoted, pos, line_end, False)
            pos = line_end
            continue
        if text.startswith("#*", pos):
            # A Cheetah block comment renders nothing between ``#*`` and ``*#``.
            block_end = text.find("*#", pos + 2)
            if block_end == -1:
                _mark(raw_flow, unquoted, pos, size, False)
                pos = size
            else:
                _mark(raw_flow, unquoted, pos, block_end + 2, False)
                pos = block_end + 2
            continue
        char = text[pos]
        word_boundary = pos == 0 or text[pos - 1] in " \t\n;|&("
        if word_boundary and char == "#" and _DIRECTIVE.match(text, pos) is not None:
            pos = _apply_hash(text, pos, raw_flow, unquoted)
            continue
        if quoted[pos]:
            pos += 1
            continue
        if word_boundary and char == "#":
            # A bash comment makes the rest of its line inert; inside quotes
            # it is literal text and was handled above.
            pos = _apply_hash(text, pos, raw_flow, unquoted)
            continue
        if char not in "'\"":
            unquoted[pos] = True
        pos += 1
    return raw_flow, unquoted


def _matching_paren(text: str, opener: int) -> Optional[int]:
    """Index of the parenthesis closing the one at ``opener``, else None."""
    depth = 0
    quote: Optional[str] = None
    pos = opener
    while pos < len(text):
        char = text[pos]
        if quote is not None:
            if char == "\\":
                pos += 2
                continue
            if char == quote:
                quote = None
        elif char in "'\"":
            quote = char
        elif char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0:
                return pos
        pos += 1
    return None


def _matching_quote(text: str, opener: int) -> Optional[int]:
    """Index of the quote closing the one at ``opener``, else None."""
    pos = opener + 1
    while pos < len(text):
        char = text[pos]
        if char == "\\":
            pos += 2
            continue
        if char == text[opener]:
            return pos
        pos += 1
    return None


def _string_literals(text: str, start: int, end: int) -> List[str]:
    """Quoted string literals in ``text[start:end]``, unescaped."""
    literals: List[str] = []
    pos = start
    while pos < end:
        char = text[pos]
        if char not in "'\"":
            pos += 1
            continue
        prefix_start = pos
        while prefix_start > start and text[prefix_start - 1] in "rRuUbBfF":
            prefix_start -= 1
        is_raw = "r" in text[prefix_start:pos] or "R" in text[prefix_start:pos]
        close = _matching_quote(text, pos)
        if close is None or close >= end:
            return literals
        body = text[pos + 1 : close]
        if not is_raw:
            body = _unescape_python(body)
        literals.append(body)
        pos = close + 1
    return literals


def _unescape_python(body: str) -> str:
    """Undo the Python string escapes that change regex interpretation."""
    result: List[str] = []
    pos = 0
    while pos < len(body):
        char = body[pos]
        if char == "\\" and pos + 1 < len(body):
            escape = body[pos + 1]
            if escape == "\\":
                result.append("\\")
            elif escape == "'":
                result.append("'")
            elif escape == '"':
                result.append('"')
            elif escape == "n":
                result.append("\n")
            elif escape == "t":
                result.append("\t")
            else:
                result.extend((char, escape))
            pos += 2
            continue
        result.append(char)
        pos += 1
    return "".join(result)


def _class_contents_safe(contents: str) -> bool:
    """True if a regex class ``[contents]`` only contains safe characters.

    Whitespace escapes (``\\s``, ``\\n``, ``\\r``, ``\\t``) do not count: the
    whitespace they allow is exactly what unquoted interpolation lets through.
    """
    pos = 0
    size = len(contents)
    while pos < size:
        char = contents[pos]
        if char == "\\":
            if pos + 1 >= size or contents[pos + 1] not in "wd.-":
                return False
            pos += 2
            continue
        if pos + 2 < size and contents[pos + 1] == "-" and contents[pos + 2] != "\\":
            if not _range_safe(char, contents[pos + 2]):
                return False
            pos += 3
            continue
        if char not in SAFE_FILENAME_CHARS:
            return False
        pos += 1
    return True


def _range_safe(start: str, end: str) -> bool:
    """True if every character in the class range ``start-end`` is safe."""
    if ord(end) - ord(start) > 256:
        return False
    return all(chr(code) in SAFE_FILENAME_CHARS for code in range(ord(start), ord(end) + 1))


def _escapes_only(contents: str, categories: Set[str]) -> bool:
    """True if ``contents`` is only escapes of the given regex categories."""
    pos = 0
    while pos < len(contents):
        if contents[pos] != "\\":
            return False
        pos += 1
        if pos >= len(contents) or contents[pos] not in categories:
            return False
        pos += 1
    return True


def _re_sub_whitelist(text: str, start: int, end: int) -> bool:
    """True if a ``re.sub`` call removes everything outside a safe character set.

    The recognized idiom replaces the complement of a character class that
    only contains ``\\w``, ``\\d``, or other safe characters, with a safe
    replacement string. The pattern must be exactly that class: anything
    after it (an anchor, a literal) leaves other unsafe characters in place.
    """
    literals = _string_literals(text, start, end)
    if len(literals) < 2:
        return False
    pattern, replacement = literals[0], literals[1]
    class_match = re.fullmatch(r"\s*\[(\^?)([^\]]*)\]", pattern)
    if class_match is None:
        return False
    negated, contents = class_match.groups()
    if negated:
        if not _class_contents_safe(contents):
            return False
    elif not _escapes_only(contents, {"W", "D"}):
        # A positive class only keeps safe characters when it removes broad
        # unsafe categories (e.g. re.sub(r'[\W]', '_', value)).
        return False
    return "\\" not in replacement and all(char in SAFE_FILENAME_CHARS for char in replacement)


def _cleaning_spans(text: str, quoted: List[bool]) -> List[Tuple[int, int]]:
    """Spans of arguments passed to recognized shell-safe cleaning helpers."""
    spans: List[Tuple[int, int]] = []
    for match in _CLEANING_CALL.finditer(text):
        open_pos = match.end() - 1
        close_pos = _matching_paren(text, open_pos)
        if close_pos is None:
            continue
        if match.group(0).startswith("shlex"):
            # shlex.quote only makes a value shell-safe when the call itself is
            # outside quotes: inside quotes its output (which contains quotes)
            # breaks out of the surrounding context.
            if quoted[open_pos]:
                continue
        elif not _re_sub_whitelist(text, open_pos + 1, close_pos):
            continue
        spans.append((open_pos + 1, close_pos))
    return spans


def _in_spans(pos: int, spans: List[Tuple[int, int]]) -> bool:
    """True if ``pos`` falls inside one of ``spans``."""
    return any(start <= pos < end for start, end in spans)


def _reference_patterns(dotted: str, loops: Dict[str, str]) -> List[str]:
    """Regexes matching how a parameter value can be referenced in a template."""
    patterns = [
        r"\$" + re.escape(dotted) + r"(?!\w)",
        r"\$\{" + re.escape(dotted) + r"\}",
    ]
    for var, expr in loops.items():
        if expr == dotted:
            patterns.append(r"\$" + re.escape(var) + r"(?!\w)")
            patterns.append(r"\$\{" + re.escape(var) + r"\}")
        elif expr and dotted.startswith(expr + "."):
            suffix = re.escape(dotted[len(expr) + 1 :])
            patterns.append(r"\$" + re.escape(var) + r"\." + suffix + r"(?!\w)")
            patterns.append(r"\$\{" + re.escape(var) + r"\." + suffix + r"\}")
    return patterns


def _preset_chars(name: str) -> str:
    """Characters allowed by a ``<valid initial=...>`` sanitizer preset."""
    if name == "default":
        return DEFAULT_VALID_CHARS
    if name == "none":
        return ""
    result = []
    for part in name.split(","):
        part = part.split(";", 1)[0].strip()
        if not part.startswith("string."):
            continue
        constant_name = STRING_PRESET_ALIASES.get(part[7:], part[7:])
        constant = getattr(string, constant_name, None)
        if isinstance(constant, str):
            result.append(constant)
    return "".join(result)


def _allowed_chars(sanitizer_el: Any) -> str:
    """Characters a ``<sanitizer>`` lets through after add/remove actions."""
    valid = ""
    for valid_el in sanitizer_el.findall("valid"):
        # Each <valid> resets the set of allowed characters.
        valid = _preset_chars(valid_el.get("initial", "default"))
        for action in valid_el:
            preset = _preset_chars(action.get("preset", "none"))
            value = action.get("value") or ""
            if action.tag == "add":
                valid = valid + "".join(char for char in preset + value if char not in valid)
            elif action.tag == "remove":
                valid = "".join(char for char in valid if char not in preset and char not in value)
    return valid


def _mapped_chars(sanitizer_el: Any) -> Dict[str, str]:
    """Character mapping a ``<sanitizer>`` applies to disallowed characters."""
    mapped: Dict[str, str] = {}
    for mapping_el in sanitizer_el.findall("mapping"):
        # Each <mapping> resets the character mapping.
        mapped = dict(MAPPING_PRESET.get(mapping_el.get("initial", "default"), {}))
        for action in mapping_el:
            source, target = action.get("source"), action.get("target")
            if action.tag == "add":
                mapped.update(MAPPING_PRESET.get(action.get("preset", "none"), {}))
                if source is not None and target is not None:
                    mapped[source] = target
            elif action.tag == "remove":
                for key in MAPPING_PRESET.get(action.get("preset", "none"), {}):
                    mapped.pop(key, None)
                if source is not None:
                    mapped.pop(source, None)
    return mapped


def _bool_value(value: Optional[str], default: bool) -> bool:
    """Parse a PermissiveBoolean XML attribute."""
    if value is None:
        return default
    lowered = value.strip().lower()
    if lowered in _TRUE_VALUES:
        return True
    if lowered in _FALSE_VALUES:
        return False
    return default


def _tool_sanitize_off(root: Any, emit: Emit) -> None:
    """Flag tools that disable sanitization for all inputs."""
    for options_el in root.findall("options"):
        if not _bool_value(options_el.get("sanitize", "true"), True):
            emit(
                "tool_sanitize_off",
                "tool",
                '<options sanitize="false"/> disables sanitization for all tool inputs; any character '
                "can reach the command, remove it or clean values at use",
                options_el,
            )


def _param_is_sanitized(param_el: Any, ptype: str) -> bool:
    """True if parameter sanitization applies to values of this parameter."""
    if ptype == "select":
        # Values of static selects come from the tool XML, not from users.
        return bool(param_el.find("options") is not None or param_el.get("dynamic_options"))
    return True


def _has_safe_validator(param_el: Any) -> bool:
    """True if a regex validator restricts the param to safe characters."""
    for validator in param_el.findall("validator"):
        if validator.get("type") != "regex" or _bool_value(validator.get("negate", "false"), False):
            continue
        expression = (validator.get("expression") or validator.text or "").strip()
        class_match = re.match(r"^\^\[([^\]]*)\](?:[+*?]|\{[^}]*\})?\$$", expression)
        if class_match is not None and _class_contents_safe(class_match.group(1)):
            return True
    return False


def _param_sanitizer_checks(params: Dict[str, _ParamInfo], emit: Emit) -> None:
    """Find disabled or too permissive parameter sanitization."""
    for dotted, param_info in params.items():
        if param_info.type not in _SANITIZED_PARAM_TYPES or not _param_is_sanitized(param_info.node, param_info.type):
            continue
        sanitizer_el = param_info.node.find("sanitizer")
        if sanitizer_el is None:
            continue
        if not _bool_value(sanitizer_el.get("sanitize", "true"), True):
            emit(
                "sanitizer_disabled",
                dotted,
                f"Sanitization is disabled for param '{dotted}' (<sanitizer sanitize=\"false\"/>); any "
                "character can reach the command, restrict the valid characters or clean the value at use",
                sanitizer_el,
            )
            continue
        valid_chars = set(_allowed_chars(sanitizer_el))
        mapping = _mapped_chars(sanitizer_el)
        # Galaxy sanitizes character by character: a valid character passes
        # through unchanged, then a mapped character is rewritten to its
        # target, and anything else becomes ``invalid_char`` -- all three
        # outputs are emitted verbatim -- so simulate that order instead of
        # assuming that mapped characters are gone.
        surviving = sorted(
            {
                *(c for c in SHELL_DANGEROUS_CHARS if c in valid_chars),
                *(
                    d
                    for source, target in mapping.items()
                    if source not in valid_chars
                    for d in target
                    if d in SHELL_DANGEROUS_CHARS
                ),
            }
        )
        if surviving:
            emit(
                "permissive_sanitizer",
                dotted,
                f"Sanitizer for param '{dotted}' allows shell metacharacters ('{''.join(surviving)}'); "
                "restrict the valid characters or clean the value at use",
                sanitizer_el,
            )
        invalid_char = sanitizer_el.get("invalid_char", "X")
        if any(c in SHELL_DANGEROUS_CHARS for c in invalid_char):
            emit(
                "dangerous_invalid_char",
                dotted,
                f"invalid_char for param '{dotted}' contains a shell metacharacter; invalid input "
                "characters are rewritten to it",
                sanitizer_el,
            )


def _raw_flow_checks(
    text: str,
    params: Dict[str, _ParamInfo],
    scopes: List[_LoopScope],
    raw_flow: List[bool],
    cleaning: List[Tuple[int, int]],
    node: Any,
    emit: Emit,
) -> None:
    """Flag raw collection identifiers and names that reach the command."""
    for match in _RAW_IDENTIFIER.finditer(text):
        base = match.group(1)
        if not raw_flow[match.start()] or _in_spans(match.start(), cleaning):
            continue
        message = (
            f"Raw value of ${base}.element_identifier reaches the command; consider "
            f"${base}.safe_element_identifier (Galaxy 26.2+) or clean it, e.g. "
            "re.sub(r'[^\\w.-]', '_', ...), before interpolating"
        )
        emit("raw_identifier", base, message, node)
    for match in _COLLECTION_NAME.finditer(text):
        base = match.group(1)
        if _expr_type(base, params, _loops_at(scopes, match.start())) != "data_collection":
            continue
        if not raw_flow[match.start()] or _in_spans(match.start(), cleaning):
            continue
        message = (
            f"${base}.name is the raw element identifier of the collection; consider "
            f"${base}.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating"
        )
        emit("raw_collection_name", base, message, node)
    for match in _COLLECTION_KEYS.finditer(text):
        base = match.group(1)
        if _expr_type(base, params, _loops_at(scopes, match.start())) != "data_collection":
            continue
        if not raw_flow[match.start()] or _in_spans(match.start(), cleaning):
            continue
        message = f"${base}.keys() returns raw element identifiers; clean each key before interpolating it"
        emit("raw_collection_keys", base, message, node)


def _param_usage_checks(
    text: str,
    params: Dict[str, _ParamInfo],
    scopes: List[_LoopScope],
    raw_flow: List[bool],
    unquoted: List[bool],
    cleaning: List[Tuple[int, int]],
    node: Any,
    emit: Emit,
) -> None:
    """Flag parameter values that reach the command unsafely."""
    for dotted, param_info in params.items():
        if param_info.type in _RAW_VALUE_PARAM_TYPES:
            message = (
                f"Value of {param_info.type} param '{dotted}' is not sanitized and reaches the command; "
                "clean the value before interpolating it"
            )
            for match in _find_reference(dotted, scopes, text):
                if raw_flow[match.start()] and not _in_spans(match.start(), cleaning):
                    emit(param_info.type, dotted, message, node)
                    break
        elif param_info.type == "text" and not _has_safe_validator(param_info.node):
            message = (
                f"Value of text param '{dotted}' is interpolated unquoted in the command; default "
                "sanitization still allows spaces, quote it or restrict the valid characters"
            )
            for match in _find_reference(dotted, scopes, text):
                if unquoted[match.start()] and not _in_spans(match.start(), cleaning):
                    emit("unquoted_param", dotted, message, node)
                    break


def _command_checks(text: str, node: Any, params: Dict[str, _ParamInfo], emit: Emit) -> None:
    """Analyze one ``<command>`` template for raw value flows."""
    quoted = _quoted_positions(text)
    cleaning = _cleaning_spans(text, quoted)
    scopes = _loop_bindings(text)
    raw_flow, unquoted = _scan_command(text, quoted, cleaning)
    _raw_flow_checks(text, params, scopes, raw_flow, cleaning, node, emit)
    _param_usage_checks(text, params, scopes, raw_flow, unquoted, cleaning, node, emit)


def _find_reference(dotted: str, scopes: List[_LoopScope], text: str) -> List[Any]:
    """All matches of the possible template references to a parameter.

    A reference through a loop variable only counts where that variable is
    bound to the parameter, so each region of the template (before the first
    loop body, or from the start of one loop body to the start of the next)
    is searched with the bindings in effect there.
    """
    matches: List[Any] = []
    for scope in [None, *scopes]:
        loops = scope.bindings if scope else {}
        for pattern in _reference_patterns(dotted, loops):
            matches.extend(match for match in re.finditer(pattern, text) if _scope_at(scopes, match.start()) is scope)
    return matches


def _security_checks(tool_source: "ToolSource") -> List[SecurityCheck]:
    """Compute the security findings for a tool source (cached per source)."""

    def calculate() -> List[SecurityCheck]:
        tool_xml = getattr(tool_source, "xml_tree", None)
        if tool_xml is None:
            return []
        root = tool_xml.getroot()
        checks: List[SecurityCheck] = []
        seen: Set[Tuple[str, str]] = set()

        def emit(linter: str, key: str, message: str, node: Any) -> None:
            if (linter, key) not in seen:
                seen.add((linter, key))
                checks.append(SecurityCheck(linter=linter, message=message, node=node))

        params = _index_params(root)
        _tool_sanitize_off(root, emit)
        _param_sanitizer_checks(params, emit)
        for command_element in root.findall("command"):
            text = "".join(command_element.itertext())
            if text.strip():
                _command_checks(text, command_element, params, emit)
        return checks

    return cached_lint_result(tool_source, "security", calculate)


class SecurityRawIdentifier(Linter):
    """Raw collection identifiers interpolated into a command."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter == "raw_identifier":
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)


class SecurityRawCollectionName(Linter):
    """Collection names or keys interpolated into a command."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter in ("raw_collection_name", "raw_collection_keys"):
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)


class SecurityGroupTag(Linter):
    """Group tag values interpolated into a command."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter == "group_tag":
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)


class SecurityDataColumn(Linter):
    """Data column (sample sheet) values interpolated into a command."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter == "data_column":
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)


class SecuritySanitizerDisabled(Linter):
    """Parameters with sanitization disabled entirely."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter == "sanitizer_disabled":
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)


class SecurityToolSanitizeOff(Linter):
    """Tools that disable sanitization for all inputs."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter == "tool_sanitize_off":
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)


class SecurityPermissiveSanitizer(Linter):
    """Sanitizer declarations that let shell metacharacters through."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter in ("permissive_sanitizer", "dangerous_invalid_char"):
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)


class SecurityUnquotedParam(Linter):
    """Text parameter values interpolated without quotes."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter == "unquoted_param":
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)
