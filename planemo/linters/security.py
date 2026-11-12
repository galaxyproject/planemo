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

Sanitizer findings report the shell metacharacters that can escape the
context of the value's worst use. Only a use in a plain shell string narrows
the set: inside single quotes only ``'`` is reported, inside double quotes
``"``, ``$``, a backtick and a backslash. Every other use gets the full set:
unquoted uses, Cheetah directive lines, strings that another shell or
interpreter parses again (``sh -c``, ``eval``, ``\\$(...)``, ``python -c``,
awk programs...), ``<configfile>`` contents, references the linter cannot
follow, and templates whose quotes do not balance. A value that never
reaches a rendered template is not reported. Line numbers are counted from
the ``<command>`` (or ``<configfile>``) tag after macro expansion and after
XML comments are removed, so they can differ from the source file.

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

from Cheetah.Parser import directiveNamesAndParsers
from galaxy.tool_util.lint import Linter
from galaxy.util import mapped_chars

from planemo.linters.util import cached_lint_result

if TYPE_CHECKING:
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource

lint_tool_types = ["*"]

# Shell metacharacters that are dangerous when they survive sanitization and
# reach a <command> unquoted. Galaxy's default sanitization removes or rewrites
# all of these, so they only survive when sanitization is disabled or
# reconfigured. A newline is not one: Galaxy replaces every newline of the
# rendered command line with a space before running it.
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
# Cheetah's inline condition (``#if $c then $a else $b#``) renders a value.
_INLINE_CONDITION = re.compile(r"[^\n]*?\bthen\b")

# Directive names Cheetah recognizes. A ``#`` followed by any other word is
# passed through literally and can then start a bash comment in the rendered
# command. ``#@`` (a decorator) is a directive only when an identifier
# follows it.
_CHEETAH_DIRECTIVE_NAMES = sorted((name for name in directiveNamesAndParsers if name != "@"), key=len, reverse=True)
_DIRECTIVE = re.compile(
    r"#(?:(?:" + "|".join(map(re.escape, _CHEETAH_DIRECTIVE_NAMES)) + r")(?![A-Za-z0-9_@-])|@(?=[A-Za-z_]))"
)
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

# Select-like types whose plain placeholder (``$p``) Galaxy sanitizes, but
# whose raw value can still be read: through ``.value`` or another attribute,
# by iterating a multiple param, or in a Python expression. Galaxy only strips
# whitespace from tags, so a group tag can hold any shell metacharacter.
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
    effect there (see ``_loop_bindings``). ``depths`` counts, for each bound
    variable, how many loops separate it from the expression in ``bindings``:
    1 for an element of a collection, 2 for an element of that element.
    """

    start: int
    bindings: Dict[str, str]
    depths: Dict[str, int]


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


def _loop_depth(expr: str, loops: Dict[str, str], depths: Dict[str, int]) -> int:
    """How many loops separate ``expr`` from the expression it resolves to.

    A bare loop variable has the depth of its binding. Any attribute access
    (``$r.coll`` for a repeat ``$r``) reaches a parameter again, at depth 0.
    """
    head, _, rest = expr.strip().partition(".")
    return depths.get(head, 0) if head in loops and not rest else 0


def _for_header(text: str, pos: int, header_end: int, earlier: _LoopScope) -> _LoopScope:
    """Bindings declared by the ``#for`` header at ``pos``, from where its body starts.

    ``earlier`` holds the loop bindings in effect at the header. The result
    extends them, and iterated expressions are resolved against them, so
    bindings always refer to parameters. Variables of loop forms that are not
    recognized are bound to an unknown ("") expression, which still shadows
    an earlier binding of the same name.
    """
    bindings = dict(earlier.bindings)
    depths = dict(earlier.depths)

    def bind(name: str, expr: str) -> None:
        bindings[name] = _resolve_loop_expr(expr, earlier.bindings)
        depths[name] = _loop_depth(expr, earlier.bindings, earlier.depths) + 1

    targets = _FOR_TARGETS.match(text, pos)
    if targets:
        for name in _LOOP_VAR.findall(targets.group(1)):
            bindings[name] = ""
            depths.pop(name, None)
        header_end = targets.end()
    enumerate_match = _FOR_ENUMERATE.match(text, pos)
    simple_match = _FOR_SIMPLE.match(text, pos)
    if enumerate_match:
        bindings[enumerate_match.group(1)] = "__index__"
        bind(enumerate_match.group(2), enumerate_match.group(3))
        header_end = enumerate_match.end()
    elif simple_match:
        bind(simple_match.group(1), simple_match.group(2))
        header_end = simple_match.end()
    return _LoopScope(header_end, bindings, depths)


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
        earlier = scopes[-1] if scopes else _LoopScope(0, {}, {})
        scopes.append(_for_header(text, token.start(), token.end(), earlier))
    return scopes


def _scope_at(scopes: List[_LoopScope], pos: int) -> Optional[_LoopScope]:
    """Loop scope started last at or before ``pos``, or None before the first loop body."""
    latest = None
    for scope in scopes:
        if scope.start > pos:
            break
        latest = scope
    return latest


def _is_collection(expr: str, params: Dict[str, _ParamInfo], scope: Optional[_LoopScope]) -> bool:
    """True if the Cheetah expression may be a dataset collection.

    ``scope`` is the loop scope in effect where the expression appears (see
    ``_scope_at``). A ``data_collection`` parameter is a collection. Each loop
    over it goes one level down its ``collection_type``: the elements of a
    ``list`` are datasets, those of a ``list:paired`` are ``paired``
    collections. Without a ``collection_type`` the input may be nested, so
    every level may be a collection.
    """
    loops = scope.bindings if scope else {}
    depths = scope.depths if scope else {}
    param = params.get(_resolve_loop_expr(expr, loops))
    if param is None or param.type != "data_collection":
        return False
    depth = _loop_depth(expr, loops, depths)
    collection_types = [kind.strip() for kind in (param.node.get("collection_type") or "").split(",")]
    return depth == 0 or any(not kind or kind.count(":") >= depth for kind in collection_types)


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
        end = _directive_end(text, pos)
        inline = _INLINE_CONDITION.match(text, directive.end(), end) is not None
        if directive.group(0)[1:] in _CONDITION_DIRECTIVES and not inline:
            _mark(raw_flow, unquoted, pos, end, False)
        else:
            # Other Cheetah directives (``#set``, ``#for``, ``#echo``...) are
            # not rendered, but can flow raw values into the command.
            for offset in range(pos, end):
                unquoted[offset] = False
        return end
    # A bash comment makes the rest of its line inert.
    _mark(raw_flow, unquoted, pos, line_end, False)
    return line_end


def _directive_end(text: str, pos: int) -> int:
    """End of the Cheetah directive starting at ``pos``.

    A directive normally runs to the end of its line, but the one-line form of
    a condition (``#if $x: echo $x``) renders the text after its colon, so the
    directive ends just past that colon. The colon is the last one outside
    strings and brackets, and only counts when text other than a comment
    follows it.
    """
    line_end = _line_end(text, pos)
    directive = _DIRECTIVE.match(text, pos)
    if directive is None or directive.group(0)[1:] not in _CONDITION_DIRECTIVES:
        return line_end
    quote: Optional[str] = None
    depth = 0
    colon = -1
    for index in range(directive.end(), line_end):
        char = text[index]
        if quote is not None:
            if char == quote:
                quote = None
        elif char in "'\"":
            quote = char
        elif char in "([{":
            depth += 1
        elif char in ")]}":
            depth -= 1
        elif text.startswith("##", index):
            break
        elif char == ":" and depth == 0:
            colon = index
    if colon != -1:
        rest = text[colon + 1 : line_end].strip()
        if rest and not rest.startswith("##"):
            return colon + 1
    return line_end


def _line_end(text: str, pos: int) -> int:
    """Position just past the end of the line containing ``pos``."""
    end = text.find("\n", pos)
    return len(text) if end == -1 else end


def _quoted_positions(text: str) -> Tuple[List[str], List[int], bool]:
    """The quote enclosing each position, where that quoted span opens, and balance.

    The first list holds ``'``, ``"``, or ``""`` when unquoted, so it is
    truthy exactly for positions inside a quoted span. The second holds the
    position of the opening quote (-1 when unquoted). The flag is False when
    a quote is still open at the end, which means the quote model is wrong
    somewhere in the template.

    The bodies of Cheetah placeholders with brackets (``${...}``,
    ``$t.replace('"', '')``) are Python code, so their quotes are skipped and
    the body takes the shell context of the placeholder.

    Cheetah removes ``##`` lines and ``#*`` blocks wherever they appear and
    processes its directives everywhere too, so those spans are skipped
    without touching the quote state: a quote opened before them stays open,
    because Cheetah may have eaten its closing quote. A bash comment (a ``#``
    at a word boundary that starts no Cheetah directive) is only recognized
    outside quotes: inside quotes it is literal text.
    """
    size = len(text)
    quoted = [""] * size
    opens = [-1] * size
    in_single = in_double = False
    opened = -1
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
        if char == "$" and (pos == 0 or text[pos - 1] != "\\"):
            # ``\$`` is Cheetah's escape for a literal ``$``, not a placeholder.
            placeholder_end = _placeholder_end(text, pos)
            if placeholder_end is not None:
                kind = "'" if in_single else ('"' if in_double else "")
                for index in range(pos, placeholder_end):
                    quoted[index] = kind
                    opens[index] = opened if kind else -1
                pos = placeholder_end
                continue
        if in_single or in_double:
            close = "'" if in_single else '"'
            quoted[pos] = close
            opens[pos] = opened
            if char == close:
                in_single = in_double = False
            elif in_double and char == "\\":
                pos = min(pos + 2, size)
                continue
            elif char == "#" and (pos == 0 or text[pos - 1] in " \t\n;|&(") and _DIRECTIVE.match(text, pos) is not None:
                # A Cheetah directive is processed even inside quotes, so its
                # content never reaches the shell and cannot close the quote.
                pos = _directive_end(text, pos)
                continue
            pos += 1
            continue
        if char in "'\"":
            in_single = char == "'"
            in_double = char == '"'
            opened = pos
        elif char == "\\":
            pos = min(pos + 2, size)
            continue
        elif char == "#" and (pos == 0 or text[pos - 1] in " \t\n;|&("):
            # A directive or bash comment: its content never reaches the shell.
            pos = _directive_end(text, pos)
            continue
        pos += 1
    return quoted, opens, not (in_single or in_double)


def _scan_command(text: str, quoted: List[str], cleaning_spans: List[Tuple[int, int]]) -> Tuple[List[bool], List[bool]]:
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
    """Index of the bracket closing the ``(``, ``[`` or ``{`` at ``opener``, else None."""
    open_char = text[opener]
    close_char = {"(": ")", "[": "]", "{": "}"}[open_char]
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
        elif char == open_char:
            depth += 1
        elif char == close_char:
            depth -= 1
            if depth == 0:
                return pos
        pos += 1
    return None


_PLACEHOLDER_NAME = re.compile(r"[A-Za-z_][\w.]*")
_PLACEHOLDER_ATTRIBUTE = re.compile(r"\.[A-Za-z_]\w*")


def _placeholder_end(text: str, pos: int) -> Optional[int]:
    """End of a Cheetah placeholder at ``pos`` whose body has brackets, else None.

    ``${...}``, ``$(...)`` and calls or subscripts such as ``$t.replace('"', '')``
    or ``$cond['x']`` are Python expressions: Cheetah evaluates them, so the
    quotes inside never reach the shell. Plain ``$name.attr`` has no quotes and
    is not skipped.
    """
    start = pos + 1
    if start < len(text) and text[start] in "{(":
        close = _matching_paren(text, start)
        return None if close is None else close + 1
    name = _PLACEHOLDER_NAME.match(text, start)
    if name is None:
        return None
    end = name.end()
    while end < len(text):
        if text[end] in "([":
            close = _matching_paren(text, end)
            if close is None:
                break
            end = close + 1
            continue
        attribute = _PLACEHOLDER_ATTRIBUTE.match(text, end)
        if attribute is None or end == name.end():
            break
        end = attribute.end()
    return end if end > name.end() else None


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


def _cleaning_spans(text: str, quoted: List[str]) -> List[Tuple[int, int]]:
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


def _reference_paths(dotted: str, loops: Dict[str, str]) -> List[str]:
    """Regexes of the dotted paths that name a parameter, directly or through a loop variable."""
    paths = [re.escape(dotted)]
    for var, expr in loops.items():
        if expr == dotted:
            paths.append(re.escape(var))
        elif expr and dotted.startswith(expr + "."):
            paths.append(re.escape(var) + r"\." + re.escape(dotted[len(expr) + 1 :]))
    return paths


def _reference_patterns(dotted: str, loops: Dict[str, str]) -> List[str]:
    """Regexes matching how a parameter value can be referenced in a template."""
    patterns = []
    for path in _reference_paths(dotted, loops):
        patterns.append(r"\$" + path + r"(?!\w)")
        patterns.append(r"\$\{" + path + r"\}")
    return patterns


def _attribute_patterns(dotted: str, loops: Dict[str, str]) -> List[str]:
    """Regexes matching an attribute, call or subscript on a parameter (``$p.value``, ``${p[0]}``)."""
    return [r"\$[{(]?\s*" + path + r"(?:\s*\.\s*[A-Za-z_]|\s*[(\[])" for path in _reference_paths(dotted, loops)]


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


# Contexts a parameter value can be used in, with the characters that let
# the value escape each of them. Inside single quotes only a quote ends the
# string; inside double quotes expansions still run. Every other context gets
# the full set: unquoted uses, Cheetah directive lines (``#set`` can flow the
# value anywhere), strings that another shell or interpreter parses again,
# ``<configfile>`` contents (which may be any language), and references the
# linter cannot follow.
_SINGLE, _DOUBLE, _UNQUOTED, _DIRECTIVE_LINE, _REPARSED, _CONFIGFILE, _UNKNOWN, _UNBALANCED = (
    "'",
    '"',
    "",
    "#",
    "reparsed",
    "configfile",
    "unknown",
    "unbalanced",
)
_CONTEXT_UNSAFE_CHARS = {
    _SINGLE: "'",
    _DOUBLE: '"$`\\',
}
_CONTEXT_LABELS = {
    _SINGLE: "inside single quotes",
    _DOUBLE: "inside double quotes",
    _UNQUOTED: "unquoted",
    _DIRECTIVE_LINE: "in a Cheetah directive",
    _REPARSED: "in a string that another shell or interpreter parses again",
    _CONFIGFILE: "in a configfile",
    _UNKNOWN: "in a form the linter cannot follow",
    _UNBALANCED: "in a template whose quotes do not balance",
}

# Text before a quoted use (from the start of the shell command where the
# quote opens) that makes the quoted string code for a second parser: command substitution,
# ``eval``, ``sh -c`` and the like, inline interpreter scripts, awk programs.
_REPARSING = re.compile(
    r"\\\$\(|`|\beval\b|\b(?:ba|da|k|z)?sh\b[^\n]*?\s-\w*c\b"
    r"|\b(?:python[\d.]*|perl|ruby|node|Rscript|R)\b[^\n]*?\s-\w*[ce]\b|\b[gmn]?awk\b"
)


class _CommandScan(NamedTuple):
    """Classification of one ``<command>`` or ``<configfile>`` template."""

    text: str
    node: Any
    label: str
    configfile: bool
    quoted: List[str]
    opens: List[int]
    balanced: bool
    raw_flow: List[bool]
    unquoted: List[bool]
    cleaning: List[Tuple[int, int]]
    scopes: List[_LoopScope]


class _Use(NamedTuple):
    """One place where a parameter value reaches a rendered template."""

    context: str
    location: str


def _scan(text: str, node: Any, label: str = "<command>", configfile: bool = False) -> _CommandScan:
    """Classify a template once for all checks."""
    quoted, opens, balanced = _quoted_positions(text)
    cleaning = _cleaning_spans(text, quoted)
    raw_flow, unquoted = _scan_command(text, quoted, cleaning)
    return _CommandScan(
        text, node, label, configfile, quoted, opens, balanced, raw_flow, unquoted, cleaning, _loop_bindings(text)
    )


def _location(text: str, pos: int, label: str = "<command>") -> str:
    """Human-readable position of ``pos``; line 1 is the line of the template's tag."""
    return f"line {text.count(chr(10), 0, pos) + 1} of {label}"


def _command_start(scan: _CommandScan, pos: int) -> int:
    """Start of the shell command containing ``pos``.

    Galaxy joins the lines of the rendered command before running it, so a
    command runs back to the last unquoted ``;``, ``|`` or ``&``, across lines.
    """
    start = pos
    while start > 0 and not (scan.text[start - 1] in ";|&" and not scan.quoted[start - 1]):
        start -= 1
    return start


def _context(scan: _CommandScan, pos: int) -> str:
    """Context of a rendered reference at ``pos``."""
    if scan.configfile:
        return _CONFIGFILE
    if scan.quoted[pos]:
        if not scan.balanced:
            return _UNBALANCED
        if _REPARSING.search(scan.text, _command_start(scan, scan.opens[pos]), pos):
            return _REPARSED
        return scan.quoted[pos]
    return _UNQUOTED if scan.unquoted[pos] else _DIRECTIVE_LINE


def _mention_pattern(dotted: str) -> "re.Pattern[str]":
    """Loose matches of anything that may refer to the parameter ``dotted``.

    This catches the forms ``_find_reference`` does not resolve, such as
    ``${x.strip()}``, an alias (``#set $r = $cond`` then ``$r.x``) or a
    lookup by name (``$cond['x']``, ``$getVar('x')``).
    """
    leaf = re.escape(dotted.rsplit(".", 1)[-1])
    return re.compile(r"\$\{?\s*(?:[A-Za-z_][\w.]*\.)?" + leaf + r"(?!\w)|(['\"])" + leaf + r"\1")


def _uses(dotted: str, scans: List[_CommandScan], references: Callable[[_CommandScan], Set[int]]) -> List[_Use]:
    """Where a parameter's value reaches the rendered templates, in order.

    References that do not render (comments, condition directives) or that
    are passed through a recognized cleaning helper are left out. A rendered
    mention of the parameter that is not one of the recognized references
    (``references`` gives the positions of all of them, for any parameter)
    counts as a use in an unknown context.
    """
    uses: List[_Use] = []
    for scan in scans:
        found = set()
        for match in _find_reference(dotted, scan.scopes, scan.text):
            pos = match.start()
            found.add(pos)
            if scan.raw_flow[pos] and not _in_spans(pos, scan.cleaning):
                uses.append(_Use(_context(scan, pos), _location(scan.text, pos, scan.label)))
        for match in _mention_pattern(dotted).finditer(scan.text):
            pos = match.start()
            if pos in found or not scan.raw_flow[pos] or _in_spans(pos, scan.cleaning) or pos in references(scan):
                continue
            uses.append(_Use(_UNKNOWN, _location(scan.text, pos, scan.label)))
    return uses


def _show_chars(chars: Set[str]) -> str:
    """Characters for a message."""
    return "".join(sorted(chars))


def _worst_use(chars: Set[str], uses: List[_Use]) -> Tuple[Set[str], Optional[_Use]]:
    """The use where most of ``chars`` are unsafe (the first on ties), with those characters."""
    worst: Set[str] = set()
    worst_use: Optional[_Use] = None
    for use in uses:
        unsafe = set(_CONTEXT_UNSAFE_CHARS.get(use.context, SHELL_DANGEROUS_CHARS))
        here = {char for char in chars if char in unsafe}
        if len(here) > len(worst):
            worst, worst_use = here, use
    return worst, worst_use


def _param_sanitizer_checks(params: Dict[str, _ParamInfo], scans: List[_CommandScan], emit: Emit) -> None:
    """Find disabled or too permissive parameter sanitization.

    A permissive sanitizer or a dangerous ``invalid_char`` is reported with
    the characters that can escape the context of its worst use. Only uses
    in a plain single- or double-quoted shell string narrow the set; values
    that never reach a rendered template are not reported.
    """
    reference_cache: Dict[int, Set[int]] = {}

    def references(scan: _CommandScan) -> Set[int]:
        key = id(scan)
        if key not in reference_cache:
            reference_cache[key] = {
                match.start() for name in params for match in _find_reference(name, scan.scopes, scan.text)
            }
        return reference_cache[key]

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
        surviving = {
            *(c for c in SHELL_DANGEROUS_CHARS if c in valid_chars),
            *(
                d
                for source, target in mapping.items()
                if source not in valid_chars
                for d in target
                if d in SHELL_DANGEROUS_CHARS
            ),
        }
        invalid_chars = set(sanitizer_el.get("invalid_char", "X"))
        if not (surviving | invalid_chars) & set(SHELL_DANGEROUS_CHARS):
            continue
        uses = _uses(dotted, scans, references)
        unsafe, use = _worst_use(surviving, uses)
        if use is not None:
            emit(
                "permissive_sanitizer",
                dotted,
                f"Sanitizer for param '{dotted}' allows shell metacharacters ('{_show_chars(unsafe)}') and the "
                f"value is used {_CONTEXT_LABELS[use.context]} ({use.location}); restrict the valid "
                "characters or clean the value at use",
                sanitizer_el,
            )
        unsafe, use = _worst_use(invalid_chars, uses)
        if use is not None:
            emit(
                "dangerous_invalid_char",
                dotted,
                f"invalid_char for param '{dotted}' contains a shell metacharacter and the value is used "
                f"{_CONTEXT_LABELS[use.context]} ({use.location}); invalid input characters are rewritten to it",
                sanitizer_el,
            )


def _raw_flow_checks(scan: _CommandScan, params: Dict[str, _ParamInfo], emit: Emit) -> None:
    """Flag raw collection identifiers and names that reach the command."""
    text, scopes, raw_flow, cleaning, node = scan.text, scan.scopes, scan.raw_flow, scan.cleaning, scan.node
    for match in _RAW_IDENTIFIER.finditer(text):
        base = match.group(1)
        if not raw_flow[match.start()] or _in_spans(match.start(), cleaning):
            continue
        message = (
            f"Raw value of ${base}.element_identifier reaches the command ({_location(text, match.start())}); "
            f"consider ${base}.safe_element_identifier (Galaxy 26.2+) or clean it, e.g. "
            "re.sub(r'[^\\w.-]', '_', ...), before interpolating"
        )
        emit("raw_identifier", base, message, node)
    for match in _COLLECTION_NAME.finditer(text):
        base = match.group(1)
        if not _is_collection(base, params, _scope_at(scopes, match.start())):
            continue
        if not raw_flow[match.start()] or _in_spans(match.start(), cleaning):
            continue
        message = (
            f"${base}.name is the raw element identifier of the collection ({_location(text, match.start())}); "
            f"consider ${base}.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating"
        )
        emit("raw_collection_name", base, message, node)
    for match in _COLLECTION_KEYS.finditer(text):
        base = match.group(1)
        if not _is_collection(base, params, _scope_at(scopes, match.start())):
            continue
        if not raw_flow[match.start()] or _in_spans(match.start(), cleaning):
            continue
        message = (
            f"${base}.keys() returns raw element identifiers ({_location(text, match.start())}); clean each "
            "key before interpolating it"
        )
        emit("raw_collection_keys", base, message, node)


_STR_CALL_OPEN = re.compile(r"(?<![\w.])str\(\s*$")
_CALL_CLOSE = re.compile(r"\s*\)")


def _expression_spans(text: str) -> List[Tuple[int, int]]:
    """Spans of the outermost Cheetah placeholders whose body is a Python expression."""
    spans: List[Tuple[int, int]] = []
    pos = text.find("$")
    while pos != -1:
        end = None if pos > 0 and text[pos - 1] == "\\" else _placeholder_end(text, pos)
        if end is not None:
            spans.append((pos, end))
        pos = text.find("$", end if end is not None else pos + 1)
    return spans


def _raw_reads(scan: _CommandScan, dotted: str) -> List[int]:
    """Positions where a template reads the raw value of a select-like parameter.

    Galaxy sanitizes a plain ``$p`` or ``${p}`` rendered as text, and
    ``str($p)`` anywhere, since both call the wrapper's ``__str__``. Any other
    read is raw: an attribute such as ``.value``, a use inside a Python
    expression (``${",".join($p)}``) and a use on a directive line, where
    ``#for`` iterates the raw values and ``#set`` may alias them.
    """
    text = scan.text
    expressions = _expression_spans(text)
    positions = {match.start() for match in _find_reference(dotted, scan.scopes, text, _attribute_patterns)}
    for match in _find_reference(dotted, scan.scopes, text):
        pos = match.start()
        if _STR_CALL_OPEN.search(text, max(0, pos - 32), pos) and _CALL_CLOSE.match(text, match.end()):
            continue
        in_expression = any(start < pos < end for start, end in expressions)
        if in_expression or not (scan.quoted[pos] or scan.unquoted[pos]):
            positions.add(pos)
    return sorted(pos for pos in positions if scan.raw_flow[pos] and not _in_spans(pos, scan.cleaning))


def _param_usage_checks(scan: _CommandScan, params: Dict[str, _ParamInfo], emit: Emit) -> None:
    """Flag parameter values that reach the command unsafely."""
    text, unquoted, cleaning = scan.text, scan.unquoted, scan.cleaning
    for dotted, param_info in params.items():
        if param_info.type in _RAW_VALUE_PARAM_TYPES:
            for pos in _raw_reads(scan, dotted)[:1]:
                message = (
                    f"Raw value of {param_info.type} param '{dotted}' reaches the command "
                    f"({_location(text, pos)}); Galaxy sanitizes only a plain placeholder, not .value, iterated "
                    "values or Python expressions: clean the value before interpolating it"
                )
                emit(param_info.type, dotted, message, scan.node)
        elif param_info.type == "text" and not _has_safe_validator(param_info.node):
            for match in _find_reference(dotted, scan.scopes, text):
                if unquoted[match.start()] and not _in_spans(match.start(), cleaning):
                    message = (
                        f"Value of text param '{dotted}' is interpolated unquoted in the command "
                        f"({_location(text, match.start())}); default sanitization still allows spaces, quote "
                        "it or restrict the valid characters"
                    )
                    emit("unquoted_param", dotted, message, scan.node)
                    break


def _find_reference(
    dotted: str,
    scopes: List[_LoopScope],
    text: str,
    patterns: Callable[[str, Dict[str, str]], List[str]] = _reference_patterns,
) -> List[Any]:
    """All matches of the possible template references to a parameter, in order.

    A reference through a loop variable only counts where that variable is
    bound to the parameter, so each region of the template (before the first
    loop body, or from the start of one loop body to the start of the next)
    is searched with the bindings in effect there. ``patterns`` builds the
    regexes of a reference from the parameter and those bindings.
    """
    matches: List[Any] = []
    for scope in [None, *scopes]:
        loops = scope.bindings if scope else {}
        for pattern in patterns(dotted, loops):
            matches.extend(match for match in re.finditer(pattern, text) if _scope_at(scopes, match.start()) is scope)
    matches.sort(key=lambda match: match.start())
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
        scans = []
        for command_element in root.findall("command"):
            text = "".join(command_element.itertext())
            if text.strip():
                scans.append(_scan(text, command_element))
        command_scans = list(scans)
        for configfile_element in root.findall("configfiles/configfile"):
            text = "".join(configfile_element.itertext())
            if text.strip():
                label = f"<configfile name=\"{configfile_element.get('name', '')}\">"
                scans.append(_scan(text, configfile_element, label, configfile=True))
        _tool_sanitize_off(root, emit)
        _param_sanitizer_checks(params, scans, emit)
        for scan in command_scans:
            _raw_flow_checks(scan, params, emit)
            _param_usage_checks(scan, params, emit)
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
    """Raw group tag values (``.value``, iteration) interpolated into a command."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext") -> None:
        for check in _security_checks(tool_source):
            if check.linter == "group_tag":
                lint_ctx.info(check.message, linter=cls.name(), node=check.node)


class SecurityDataColumn(Linter):
    """Raw data column values (``.value``, iteration) interpolated into a command."""

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
