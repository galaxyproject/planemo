from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest
from galaxy.tool_util.lint import (
    lint_tool_source_with_modules,
    LintContext,
)
from galaxy.tool_util.parser import get_tool_source

from planemo.lint import build_lint_args
from planemo.linters import (
    biocontainer_registered,
    conda_requirements,
    doi,
    security,
    urls,
)
from planemo.shed_lint import (
    lint_shed_version,
    lint_tool_dependencies_urls,
    REPOSITORY_LINTER_NAMES,
)
from .test_utils import (
    PROJECT_TEMPLATES_DIR,
    TEST_REPOS_DIR,
    TEST_TOOLS_DIR,
)


def _tool_source(directory, name):
    return get_tool_source(str(Path(directory) / name))


def _messages(lint_ctx):
    return [(message.level, message.message, message.linter) for message in lint_ctx.message_list]


def test_doi_linter_checks_each_citation_once():
    tool_source = _tool_source(TEST_TOOLS_DIR, "invalid_doi.xml")
    responses = [mock.Mock(status_code=200), mock.Mock(status_code=404)]
    lint_ctx = LintContext("all")

    with mock.patch("planemo.linters.doi.requests.get", side_effect=responses) as get:
        lint_tool_source_with_modules(lint_ctx, tool_source, [doi])

    assert get.call_count == 2
    assert all(call.kwargs["timeout"] == 5 for call in get.call_args_list)
    assert _messages(lint_ctx) == [
        ("error", "doi:10.1101/666666 is not a valid DOI", "DoiInvalid"),
        (
            "error",
            "doi:10.1101/014043 is valid, but Galaxy expects DOI without 'doi:' prefix",
            "DoiPrefix",
        ),
    ]


def test_doi_linter_can_be_skipped_by_module_name():
    tool_source = _tool_source(TEST_TOOLS_DIR, "invalid_doi.xml")
    lint_ctx = LintContext("all", skip_types=["doi"])

    with mock.patch("planemo.linters.doi.requests.get") as get:
        lint_tool_source_with_modules(lint_ctx, tool_source, [doi])

    get.assert_not_called()
    assert lint_ctx.message_list == []


def test_doi_linter_reports_request_failures_without_crashing():
    tool_source = _tool_source(TEST_TOOLS_DIR, "invalid_doi.xml")
    lint_ctx = LintContext("all")

    with mock.patch(
        "planemo.linters.doi.requests.get",
        side_effect=doi.requests.ConnectionError("offline"),
    ) as get:
        lint_tool_source_with_modules(lint_ctx, tool_source, [doi])

    assert get.call_count == 2
    assert _messages(lint_ctx) == [
        (
            "warning",
            "Error 'offline' accessing https://doi.org/10.1101/014043",
            "DoiUnexpectedResponse",
        ),
        (
            "warning",
            "Error 'offline' accessing https://doi.org/10.1101/666666",
            "DoiUnexpectedResponse",
        ),
    ]


def test_url_linter_checks_each_help_url_once():
    tool_source = _tool_source(TEST_TOOLS_DIR, "url.xml")
    response = mock.Mock(status_code=200)
    response.iter_content.side_effect = lambda _: iter([b"content"])
    lint_ctx = LintContext("all")

    with mock.patch("planemo.lint.requests.get", return_value=response) as get:
        lint_tool_source_with_modules(lint_ctx, tool_source, [urls])

    assert get.call_count == 2
    assert all(call.kwargs["timeout"] == 5 for call in get.call_args_list)
    assert _messages(lint_ctx) == [
        ("info", "URL OK http://galaxyproject.org/", "URLValid"),
        ("info", "URL OK https://galaxyproject.org/", "URLValid"),
    ]


def test_tool_dependency_url_linter_reuses_hardened_url_check():
    realized_repository = SimpleNamespace(real_path=str(Path(TEST_REPOS_DIR) / "package_1"))
    response = mock.Mock(status_code=200)
    response.iter_content.side_effect = lambda _: iter([b"content"])
    lint_ctx = LintContext("all")

    with mock.patch("planemo.lint.requests.get", return_value=response) as get:
        lint_tool_dependencies_urls(realized_repository, lint_ctx)

    assert get.call_count == 5
    assert all(call.kwargs["timeout"] == 5 for call in get.call_args_list)
    assert all(message.level == "info" for message in lint_ctx.message_list)


def test_repository_linter_names_are_valid_skip_targets():
    ctx = mock.Mock(global_config={})

    with mock.patch("planemo.lint.error") as error:
        lint_args = build_lint_args(
            ctx,
            skip=["shed_remote_repository_url"],
            extra_linter_names=REPOSITORY_LINTER_NAMES,
        )

    error.assert_not_called()
    assert lint_args["skip_types"] == ["shed_remote_repository_url"]


def test_conda_linter_checks_each_requirement_once():
    tool_source = _tool_source(TEST_TOOLS_DIR, "bwa_invalid_version.xml")
    lint_ctx = LintContext("all")

    with mock.patch("planemo.linters.conda_requirements.best_practice_search", return_value=(None, None)) as search:
        lint_tool_source_with_modules(lint_ctx, tool_source, [conda_requirements])

    search.assert_called_once()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "Requirement [bwa@0.4.12] doesn't match any recipe in a best practice "
            "Conda channel [['conda-forge', 'bioconda']].",
            "CondaRequirementMissing",
        )
    ]


def test_conda_linter_distinguishes_inexact_from_missing_requirements():
    tool_source = _tool_source(TEST_TOOLS_DIR, "bwa_invalid_version.xml")
    lint_ctx = LintContext("all")
    best_hit = {"channel": "bioconda", "version": "0.7.10"}

    with mock.patch(
        "planemo.linters.conda_requirements.best_practice_search",
        return_value=(best_hit, False),
    ) as search:
        lint_tool_source_with_modules(lint_ctx, tool_source, [conda_requirements])

    search.assert_called_once()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "Requirement [bwa@0.4.12] doesn't exactly match available version "
            "[0.7.10] in best practice Conda channel [bioconda].",
            "CondaRequirementInexact",
        )
    ]


def test_conda_linter_reports_missing_requirements_without_searching():
    tool_source = _tool_source(TEST_TOOLS_DIR, "bwa_without_requirements.xml")
    lint_ctx = LintContext("all")

    with mock.patch("planemo.linters.conda_requirements.best_practice_search") as search:
        lint_tool_source_with_modules(lint_ctx, tool_source, [conda_requirements])

    search.assert_not_called()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "No valid package requirement tags found to check against Conda.",
            "CondaRequirementsMissing",
        )
    ]


def test_biocontainer_linter_resolves_requirement_set_once():
    tool_source = _tool_source(Path(PROJECT_TEMPLATES_DIR) / "seqtk_complete", "seqtk_seq.xml")
    lint_ctx = LintContext("all")

    with mock.patch(
        "planemo.linters.biocontainer_registered.mulled_container_name",
        return_value="quay.io/biocontainers/seqtk:1.2",
    ) as resolve:
        lint_tool_source_with_modules(lint_ctx, tool_source, [biocontainer_registered])

    resolve.assert_called_once()
    assert _messages(lint_ctx) == [
        (
            "info",
            "BioContainer best-practice container found [quay.io/biocontainers/seqtk:1.2].",
            "BiocontainerValid",
        )
    ]


def test_biocontainer_linter_reports_missing_requirements_without_resolving():
    tool_source = _tool_source(TEST_TOOLS_DIR, "bwa_without_requirements.xml")
    lint_ctx = LintContext("all")

    with mock.patch("planemo.linters.biocontainer_registered.mulled_container_name") as resolve:
        lint_tool_source_with_modules(lint_ctx, tool_source, [biocontainer_registered])

    resolve.assert_not_called()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "No valid package requirement tags found to infer BioContainer from.",
            "BiocontainerRequirementsMissing",
        )
    ]


def test_biocontainer_linter_reports_missing_container_once():
    tool_source = _tool_source(Path(PROJECT_TEMPLATES_DIR) / "seqtk_complete", "seqtk_seq.xml")
    lint_ctx = LintContext("all")

    with mock.patch(
        "planemo.linters.biocontainer_registered.mulled_container_name",
        return_value=None,
    ) as resolve:
        lint_tool_source_with_modules(lint_ctx, tool_source, [biocontainer_registered])

    resolve.assert_called_once()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "Failed to find a BioContainer registered for these requirements.",
            "BiocontainerMissing",
        )
    ]


def test_security_linter_flags_raw_identifier():
    tool_source = _tool_source(TEST_TOOLS_DIR, "security_raw_identifier.xml")
    lint_ctx = LintContext("all")

    lint_tool_source_with_modules(lint_ctx, tool_source, [security])

    assert _messages(lint_ctx) == [
        (
            "info",
            "$collection1.name is the raw element identifier of the collection; consider "
            "$collection1.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating",
            "SecurityRawCollectionName",
        ),
        (
            "info",
            "$pair.name is the raw element identifier of the collection; consider "
            "$pair.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating",
            "SecurityRawCollectionName",
        ),
        (
            "info",
            "$collection1.keys() returns raw element identifiers; clean each key before interpolating it",
            "SecurityRawCollectionName",
        ),
        (
            "info",
            r"Raw value of $collection1.element_identifier reaches the command; consider "
            r"$collection1.safe_element_identifier (Galaxy 26.2+) or clean it, e.g. "
            r"re.sub(r'[^\w.-]', '_', ...), before interpolating",
            "SecurityRawIdentifier",
        ),
    ]


def test_security_linter_flags_disabled_sanitization():
    tool_source = _tool_source(TEST_TOOLS_DIR, "security_sanitize_disabled.xml")
    lint_ctx = LintContext("all")

    lint_tool_source_with_modules(lint_ctx, tool_source, [security])

    assert _messages(lint_ctx) == [
        (
            "info",
            "Sanitizer for param 'permissive_text' allows shell metacharacters "
            "('\"#$&';<>\\`|'); restrict the valid characters or clean the value at use",
            "SecurityPermissiveSanitizer",
        ),
        (
            "info",
            "Sanitization is disabled for param 'free_text' (<sanitizer sanitize=\"false\"/>); any "
            "character can reach the command, restrict the valid characters or clean the value at use",
            "SecuritySanitizerDisabled",
        ),
        (
            "info",
            '<options sanitize="false"/> disables sanitization for all tool inputs; any character '
            "can reach the command, remove it or clean values at use",
            "SecurityToolSanitizeOff",
        ),
    ]


def test_security_linter_flags_column_and_tag_values():
    tool_source = _tool_source(TEST_TOOLS_DIR, "security_data_column_group_tag.xml")
    lint_ctx = LintContext("all")

    lint_tool_source_with_modules(lint_ctx, tool_source, [security])

    assert _messages(lint_ctx) == [
        (
            "info",
            "Value of data_column param 'columnvalue' is not sanitized and reaches the command; "
            "clean the value before interpolating it",
            "SecurityDataColumn",
        ),
        (
            "info",
            "Value of group_tag param 'tagvalue' is not sanitized and reaches the command; "
            "clean the value before interpolating it",
            "SecurityGroupTag",
        ),
    ]


def test_security_linter_flags_unquoted_text():
    tool_source = _tool_source(TEST_TOOLS_DIR, "security_unquoted_param.xml")
    lint_ctx = LintContext("all")

    lint_tool_source_with_modules(lint_ctx, tool_source, [security])

    # 'checked_text' is restricted by a regex validator and not reported.
    assert _messages(lint_ctx) == [
        (
            "info",
            "Value of text param 'textparam' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_stays_quiet_on_clean_tool():
    tool_source = _tool_source(TEST_TOOLS_DIR, "security_clean_tool.xml")
    lint_ctx = LintContext("all")

    lint_tool_source_with_modules(lint_ctx, tool_source, [security])

    assert lint_ctx.message_list == []


def test_security_linter_can_be_skipped_by_module_name():
    tool_source = _tool_source(TEST_TOOLS_DIR, "security_raw_identifier.xml")
    lint_ctx = LintContext("all", skip_types=["security"])

    lint_tool_source_with_modules(lint_ctx, tool_source, [security])

    assert lint_ctx.message_list == []


SECURITY_LOOP_TOOL = """<tool id="security_loop_scope" name="Security Loop Scope" version="1.0">
    <command><![CDATA[
{command}
    ]]></command>
    <inputs>
        <param name="coll" type="data_collection" collection_type="list"/>
        <repeat name="rep" title="r">
            <param name="txt" type="text"/>
        </repeat>
    </inputs>
    <outputs />
</tool>
"""


def _lint_loop_tool(tmp_path, command):
    (tmp_path / "loop_tool.xml").write_text(SECURITY_LOOP_TOOL.format(command=command))
    lint_ctx = LintContext("all")
    lint_tool_source_with_modules(lint_ctx, _tool_source(tmp_path, "loop_tool.xml"), [security])
    return _messages(lint_ctx)


@pytest.mark.parametrize(
    "command",
    [
        pytest.param(
            "#for $i in $coll:\necho '${i.name}';\n#end for\n#for $i in $rep:\ntool --opt $i.txt;\n#end for",
            id="collection-loop-then-repeat-loop",
        ),
        pytest.param(
            "#for $i in $rep:\ntool --opt $i.txt;\n#end for\n#for $i in $coll:\necho '${i.name}';\n#end for",
            id="repeat-loop-then-collection-loop",
        ),
        pytest.param(
            "#for $n, $i in enumerate($coll):\necho '${i.name}' $n;\n#end for\n"
            "#for $i in $rep:\ntool --opt $i.txt;\n#end for",
            id="enumerate-loop-then-repeat-loop",
        ),
    ],
)
def test_security_linter_scopes_loop_variables_to_their_body(tmp_path, command):
    assert _lint_loop_tool(tmp_path, command) == [
        (
            "info",
            "$i.name is the raw element identifier of the collection; consider "
            "$i.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating",
            "SecurityRawCollectionName",
        ),
        (
            "info",
            "Value of text param 'rep.txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


# Cheetah compiles #for into a Python for, so a loop variable keeps its last
# value after #end for, until a later loop binds the same name again.
@pytest.mark.parametrize(
    "command",
    [
        pytest.param(
            "#for $i in $coll:\necho hi;\n#end for\necho '$i.name';",
            id="reference-after-loop",
        ),
        pytest.param(
            "#for $i in $coll:\necho hi;\n#end for\n#for $j in $rep:\necho hi;\n#end for\necho '$i.name';",
            id="reference-after-loop-of-another-variable",
        ),
        pytest.param(
            "#for $i in $rep:\necho hi;\n#end for\n#for $i in $coll:\necho hi;\n#end for\necho '$i.name';",
            id="later-loop-rebinds-the-name",
        ),
        pytest.param(
            # After the inner #end for, $i holds the last collection element,
            # not the outer repeat, so $i.txt no longer reaches rep.txt.
            "#for $i in $rep:\n#for $i in $coll:\necho '${i.name}';\n#end for\ntool --opt $i.txt;\n#end for",
            id="inner-loop-shadows-outer-variable",
        ),
    ],
)
def test_security_linter_keeps_loop_variable_bound_after_its_loop(tmp_path, command):
    assert _lint_loop_tool(tmp_path, command) == [
        (
            "info",
            "$i.name is the raw element identifier of the collection; consider "
            "$i.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating",
            "SecurityRawCollectionName",
        ),
    ]


@pytest.mark.parametrize(
    "command",
    [
        pytest.param(
            "#for $i in $coll:\n#for $i in range(3):\necho '${i.name}';\n#end for\n#end for",
            id="inner-loop-of-unknown-form-shadows-outer-variable",
        ),
        pytest.param(
            "#for $i in $rep:\ntool --opt '$i.txt';\n#end for\n#for $i in $coll:\ntool --opt $i.txt;\n#end for",
            id="same-name-in-another-loop-is-another-value",
        ),
        pytest.param(
            "#for $i in $coll:\necho hi;\n#end for\n#for $i in $rep:\necho hi;\n#end for\necho '$i.name';",
            id="later-loop-rebinds-the-name-to-a-repeat",
        ),
        pytest.param(
            "echo '$i.name';\n#for $i in $coll:\necho hi;\n#end for",
            id="reference-before-loop",
        ),
        pytest.param(
            "#* #for $i in $coll: *#\n## #for $i in $coll:\n# for $i in $coll:\necho '$i.name';",
            id="commented-out-or-literal-for-binds-nothing",
        ),
    ],
)
def test_security_linter_resolves_loop_variables_only_where_bound(tmp_path, command):
    assert _lint_loop_tool(tmp_path, command) == []


def test_security_loop_bindings_persist_after_their_loop():
    # Loop variables keep their last value after #end for, so past the inner
    # loop's body start every position sees the bindings of both loops.
    text = "$a\n#for $a in $coll:\n## end for comment\n#for $b in $a:\n$b\n#end for\n$a\n#end for\n$a\n"

    scopes = security._loop_bindings(text)

    assert security._loops_at(scopes, text.index("$a\n")) == {}
    assert security._loops_at(scopes, text.index("## end")) == {"a": "coll"}
    assert security._loops_at(scopes, text.index("$b\n")) == {"a": "coll", "b": "coll"}
    assert security._loops_at(scopes, text.index("$a\n#end for")) == {"a": "coll", "b": "coll"}
    assert security._loops_at(scopes, text.rindex("$a")) == {"a": "coll", "b": "coll"}


SECURITY_TEMPLATE_TOOL = """<tool id="security_template_tool" name="Security Template Tool" version="1.0">
    <command><![CDATA[
{command}
    ]]></command>
    <inputs>
{inputs}
    </inputs>
    <outputs />
</tool>
"""


def _lint_template_tool(tmp_path, command, inputs):
    (tmp_path / "template_tool.xml").write_text(SECURITY_TEMPLATE_TOOL.format(command=command, inputs=inputs))
    lint_ctx = LintContext("all")
    lint_tool_source_with_modules(lint_ctx, _tool_source(tmp_path, "template_tool.xml"), [security])
    return _messages(lint_ctx)


def test_security_linter_ignores_name_on_plain_data_input(tmp_path):
    # .name on a plain data param is the sanitized upload name, not a raw
    # element identifier, so it is not flagged.
    messages = _lint_template_tool(
        tmp_path,
        "ln -s '$input1.name' 'out/'",
        '<param name="input1" type="data" format="txt"/>',
    )

    assert messages == []


def test_security_linter_ignores_if_conditions_but_not_bodies(tmp_path):
    # References in a Cheetah #if condition are not shell commands, but the
    # body of the directive is. The condition holds a collection name so that
    # treating the condition as shell text fails this test.
    messages = _lint_template_tool(
        tmp_path,
        "#if $coll.name:\ntool --opt $txt;\n#end if",
        '<param name="coll" type="data_collection" collection_type="list"/>' '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_cheetah_line_comment_silences_only_its_line(tmp_path):
    # A ## comment renders nothing on its own line; the next line still runs.
    messages = _lint_template_tool(
        tmp_path,
        "## note\ntool --opt $txt",
        '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_block_comment_is_inert(tmp_path):
    # A #* ... *# comment renders nothing, so a reference inside it is not
    # reported, while the line after it is.
    messages = _lint_template_tool(
        tmp_path,
        "#* '$coll.element_identifier' *#\ntool --opt $txt",
        '<param name="coll" type="data_collection" collection_type="list"/>' '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_only_treats_word_boundary_hash_as_comment(tmp_path):
    # A # in the middle of a word is a literal character, not a comment.
    messages = _lint_template_tool(
        tmp_path,
        "tool --comment=# --opt $txt",
        '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_bash_comment_silences_rest_of_line_only(tmp_path):
    # A bash comment hides what follows it on the same line, but the next
    # line of the rendered script still runs.
    messages = _lint_template_tool(
        tmp_path,
        "echo hi # '$coll.element_identifier'\ntool --opt $txt",
        '<param name="coll" type="data_collection" collection_type="list"/>' '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_flags_brace_reference(tmp_path):
    # ${txt} is another way to reference a parameter value.
    messages = _lint_template_tool(
        tmp_path,
        "tool --opt ${txt}",
        '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_flags_shlex_quote_inside_quotes(tmp_path):
    # Inside quotes, shlex.quote's output breaks out of the surrounding
    # context, so the call is not recognized as cleaning.
    messages = _lint_template_tool(
        tmp_path,
        "echo '$( shlex.quote($coll.element_identifier) )'",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == [
        (
            "info",
            r"Raw value of $coll.element_identifier reaches the command; consider "
            r"$coll.safe_element_identifier (Galaxy 26.2+) or clean it, e.g. "
            r"re.sub(r'[^\w.-]', '_', ...), before interpolating",
            "SecurityRawIdentifier",
        ),
    ]


def test_security_linter_flags_disabled_sanitization_on_hidden_param(tmp_path):
    # Hidden params carry user-controlled values, so their sanitization is
    # checked like any other text-like param.
    messages = _lint_template_tool(
        tmp_path,
        "tool --opt $h",
        '<param name="h" type="hidden"><sanitizer sanitize="false"/></param>',
    )

    assert messages == [
        (
            "info",
            "Sanitization is disabled for param 'h' (<sanitizer sanitize=\"false\"/>); any "
            "character can reach the command, restrict the valid characters or clean the value at use",
            "SecuritySanitizerDisabled",
        ),
    ]


def test_security_linter_binds_enumerate_variables(tmp_path):
    # The enumerate() form binds the value variable to the collection on its
    # own, without a second loop.
    messages = _lint_loop_tool(tmp_path, "#for $n, $i in enumerate($coll):\necho '$i.name' '$n';\n#end for")

    assert messages == [
        (
            "info",
            "$i.name is the raw element identifier of the collection; consider "
            "$i.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating",
            "SecurityRawCollectionName",
        ),
    ]


def test_security_linter_derives_param_name_from_argument(tmp_path):
    # A param declared only with argument= is registered under the name
    # Galaxy derives from the flag.
    messages = _lint_template_tool(
        tmp_path,
        "tool --opt $my_opt",
        '<param argument="--my-opt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'my_opt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_flags_valid_character_that_is_also_mapped(tmp_path):
    # A character in the valid set passes through even when the mapping would
    # rewrite it: Galaxy checks the valid set first.
    messages = _lint_template_tool(
        tmp_path,
        "echo '$txt'",
        '<param name="txt" type="text"><sanitizer>\n'
        '<valid initial="default"><add value="&apos;"/></valid>\n'
        '<mapping initial="default"/>\n'
        "</sanitizer></param>",
    )

    assert messages == [
        (
            "info",
            "Sanitizer for param 'txt' allows shell metacharacters ('''); restrict the valid "
            "characters or clean the value at use",
            "SecurityPermissiveSanitizer",
        ),
    ]


def test_security_linter_flags_mapping_target_metacharacter(tmp_path):
    # The mapping target is emitted verbatim for values that are not in the
    # valid set, so a metacharacter target reaches the command.
    messages = _lint_template_tool(
        tmp_path,
        "echo '$txt'",
        '<param name="txt" type="text"><sanitizer>\n'
        '<mapping initial="none"><add source="&apos;" target="&apos;"/></mapping>\n'
        "</sanitizer></param>",
    )

    assert messages == [
        (
            "info",
            "Sanitizer for param 'txt' allows shell metacharacters ('''); restrict the valid "
            "characters or clean the value at use",
            "SecurityPermissiveSanitizer",
        ),
    ]


def test_security_linter_flags_dangerous_invalid_char_despite_mapping(tmp_path):
    # invalid_char is emitted for every character that is neither valid nor
    # mapped, so it reaches the command regardless of the mapping.
    messages = _lint_template_tool(
        tmp_path,
        "echo '$txt'",
        '<param name="txt" type="text"><sanitizer invalid_char="&apos;">\n'
        '<mapping initial="default"/>\n'
        "</sanitizer></param>",
    )

    assert messages == [
        (
            "info",
            "invalid_char for param 'txt' contains a shell metacharacter; invalid input characters "
            "are rewritten to it",
            "SecurityPermissiveSanitizer",
        ),
    ]


def test_security_linter_accepts_word_class_validator_with_dots(tmp_path):
    # Word and digit escapes, dots and dashes stay safe in a validator class,
    # while whitespace escapes do not (see the backslash-s case above).
    messages = _lint_template_tool(
        tmp_path,
        "tool --opt $txt",
        r'<param name="txt" type="text"><validator type="regex">^[\w.-]+$</validator></param>',
    )

    assert messages == []


def test_security_linter_flags_validator_with_whitespace_escape(tmp_path):
    # \s allows the spaces that unquoted interpolation relies on, so it does
    # not restrict the param to safe characters.
    messages = _lint_template_tool(
        tmp_path,
        "tool --opt $txt",
        r'<param name="txt" type="text"><validator type="regex">^[\w\s]+$</validator></param>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_flags_unsafe_re_sub_class(tmp_path):
    # A re.sub class that keeps shell metacharacters behind (;) is not
    # recognized as cleaning.
    messages = _lint_template_tool(
        tmp_path,
        "tool --opt '$( re.sub(r'[^\\w;]', '_', $coll.element_identifier) )'",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == [
        (
            "info",
            r"Raw value of $coll.element_identifier reaches the command; consider "
            r"$coll.safe_element_identifier (Galaxy 26.2+) or clean it, e.g. "
            r"re.sub(r'[^\w.-]', '_', ...), before interpolating",
            "SecurityRawIdentifier",
        ),
    ]


def test_security_linter_cheetah_comments_inside_quotes_do_not_flip_quote_state(tmp_path):
    # Cheetah removes ## and #* comments wherever they appear, even between
    # shell quotes, so their quotes cannot pair with the shell's: an
    # apostrophe inside a comment must not close the surrounding quote.
    # P1: the ## comment's apostrophe must not close the awk script's quote,
    # which would leave $txt quoted and hide the finding.
    messages = _lint_template_tool(
        tmp_path,
        "awk '\n## skip the header, it's not data\nNR > 1 { print }\n' $txt > out.txt",
        '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]

    # Same for a #* block comment.
    messages = _lint_template_tool(
        tmp_path,
        "awk '\n#* it's a block *#\nNR > 1\n' $txt",
        '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_comment_apostrophe_keeps_shlex_quote_untrusted(tmp_path):
    # The comment's apostrophe must not make a later shlex.quote(...) call
    # look like it sits outside quotes: it does not, so the call is not
    # recognized as cleaning and the raw identifier is reported.
    messages = _lint_template_tool(
        tmp_path,
        "python -c '\n## it's fine\nprint(1)\n' && echo '$( shlex.quote($coll.element_identifier) )'",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == [
        (
            "info",
            r"Raw value of $coll.element_identifier reaches the command; consider "
            r"$coll.safe_element_identifier (Galaxy 26.2+) or clean it, e.g. "
            r"re.sub(r'[^\w.-]', '_', ...), before interpolating",
            "SecurityRawIdentifier",
        ),
    ]


def test_security_linter_directive_line_inside_quotes_does_not_flip_quote_state(tmp_path):
    # Cheetah processes its directives everywhere, so a #set line between
    # shell quotes (holding an apostrophe) renders nothing: the following
    # $txt stays unquoted and the directive's quotes cannot pair with it.
    messages = _lint_template_tool(
        tmp_path,
        "python3 -c '\n#set $q = \"it's\"\nprint(1)\n' $txt",
        '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_cheetah_comment_inside_quotes_hides_reference(tmp_path):
    # A ## comment between double quotes renders nothing (Cheetah does not
    # know about shell quotes), so the reference inside it is not reported.
    messages = _lint_template_tool(
        tmp_path,
        'echo "x ## $coll.element_identifier"',
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == []


def test_security_linter_unknown_hash_word_is_bash_comment(tmp_path):
    # A word after # that is no Cheetah directive name is passed through by
    # Cheetah and starts a bash comment, so the rest of the line is inert.
    messages = _lint_template_tool(
        tmp_path,
        "tool #x '$coll.element_identifier'",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == []


def test_security_linter_binds_loop_variables_without_dollar(tmp_path):
    # Cheetah accepts '#for i in $coll:' without the $ on the loop variable.
    messages = _lint_template_tool(
        tmp_path,
        "#for i in $coll:\necho '$i.name';\n#end for",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == [
        (
            "info",
            "$i.name is the raw element identifier of the collection; consider "
            "$i.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating",
            "SecurityRawCollectionName",
        ),
    ]


def test_security_linter_binds_enumerate_variables_without_dollar(tmp_path):
    # The enumerate() form binds both loop variables without the $.
    messages = _lint_template_tool(
        tmp_path,
        "#for n, i in enumerate($coll):\necho '$i.name' $n;\n#end for",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == [
        (
            "info",
            "$i.name is the raw element identifier of the collection; consider "
            "$i.safe_element_identifier (Galaxy 26.2+) or clean it before interpolating",
            "SecurityRawCollectionName",
        ),
    ]


def test_security_linter_midword_cheetah_comment_silences_its_line(tmp_path):
    # Cheetah starts a ## comment wherever one appears, not only at a word
    # boundary: everything after it on the line renders nothing.
    messages = _lint_template_tool(
        tmp_path,
        "tool --x=1## '$coll.element_identifier'",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == []

    # ... and its quotes cannot pair with later shell quotes.
    messages = _lint_template_tool(
        tmp_path,
        "tool --x=1## don't\ntool --opt $txt",
        '<param name="txt" type="text"/>',
    )

    assert messages == [
        (
            "info",
            "Value of text param 'txt' is interpolated unquoted in the command; default "
            "sanitization still allows spaces, quote it or restrict the valid characters",
            "SecurityUnquotedParam",
        ),
    ]


def test_security_linter_flags_mapping_target_metacharacter_from_safe_source(tmp_path):
    # A mapped character whose source is not itself a metacharacter still
    # emits its target verbatim for values that are not in the valid set.
    messages = _lint_template_tool(
        tmp_path,
        "echo '$txt'",
        '<param name="txt" type="text"><sanitizer><mapping initial="none">'
        '<add source="@" target=";"/></mapping></sanitizer></param>',
    )

    assert messages == [
        (
            "info",
            "Sanitizer for param 'txt' allows shell metacharacters (';'); restrict the valid "
            "characters or clean the value at use",
            "SecurityPermissiveSanitizer",
        ),
    ]


def test_security_linter_flags_multicharacter_invalid_char(tmp_path):
    # An invalid_char holding several characters is dangerous when any of
    # them is a shell metacharacter.
    messages = _lint_template_tool(
        tmp_path,
        "echo '$txt'",
        '<param name="txt" type="text"><sanitizer invalid_char="&apos;;"/></param>',
    )

    assert messages == [
        (
            "info",
            "invalid_char for param 'txt' contains a shell metacharacter; invalid input "
            "characters are rewritten to it",
            "SecurityPermissiveSanitizer",
        ),
    ]


def test_security_linter_flags_re_sub_with_trailing_anchor(tmp_path):
    # A re.sub pattern that only starts with the safe class keeps other
    # characters (here through the $ anchor) in place, so it is not cleaning.
    messages = _lint_template_tool(
        tmp_path,
        "tool --opt '$( re.sub(r'[^\\w.-]+$', '', $coll.element_identifier) )'",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == [
        (
            "info",
            r"Raw value of $coll.element_identifier reaches the command; consider "
            r"$coll.safe_element_identifier (Galaxy 26.2+) or clean it, e.g. "
            r"re.sub(r'[^\w.-]', '_', ...), before interpolating",
            "SecurityRawIdentifier",
        ),
    ]


def test_security_linter_bash_comment_apostrophe_keeps_shlex_quote_untrusted(tmp_path):
    # A bash comment only affects its own line, and its quotes are not
    # rendered: the apostrophe in '# don't touch' must not pair with the
    # quote around the shlex.quote call on the next line, which stays inside
    # quotes and is not recognized as cleaning.
    messages = _lint_template_tool(
        tmp_path,
        "# don't touch\necho '$( shlex.quote($coll.element_identifier) )'",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == [
        (
            "info",
            r"Raw value of $coll.element_identifier reaches the command; consider "
            r"$coll.safe_element_identifier (Galaxy 26.2+) or clean it, e.g. "
            r"re.sub(r'[^\w.-]', '_', ...), before interpolating",
            "SecurityRawIdentifier",
        ),
    ]


def test_security_linter_accepts_unquoted_shlex_quote(tmp_path):
    # An unquoted shlex.quote call is the recognized cleaning idiom.
    messages = _lint_template_tool(
        tmp_path,
        "tool --opt $( shlex.quote($coll.element_identifier) )",
        '<param name="coll" type="data_collection" collection_type="list"/>',
    )

    assert messages == []


SHED_VERSION_TOOL = """<tool id="shed_version_tool" name="Shed Version Tool" version="1.0+galaxy10">
    <command>echo hi</command>
    <inputs />
    <outputs />
</tool>
"""


@dataclass
class _ShedRepository:
    path: str
    owner: str = "iuc"
    name: str = "shed_version"


class _FakeShedRepositories:
    def __init__(self, ts_version):
        self.ts_version = ts_version

    def get_repositories(self, name, owner):
        return [{"name": name, "owner": owner}]

    def get_ordered_installable_revisions(self, name, owner):
        return ["abc123"]

    def get_repository_revision_install_info(self, name, owner, revision):
        return {}, {"valid_tools": [{"id": "shed_version_tool", "version": self.ts_version}]}, {}


class _FakeToolShed:
    def __init__(self, ts_version):
        self.repositories = _FakeShedRepositories(ts_version)


@pytest.mark.parametrize(
    "ts_version, expected",
    [
        ("1.0+galaxy9", []),
        (
            "1.0+galaxy10",
            [
                (
                    "error",
                    "shed_version_tool: version 1.0+galaxy10 is less or equal than version of the latest installable revision 1.0+galaxy10",
                    "ShedVersion",
                )
            ],
        ),
    ],
)
def test_shed_version_compares_galaxy_suffix_numerically(tmp_path, ts_version, expected):
    (tmp_path / "shed_version_tool.xml").write_text(SHED_VERSION_TOOL)
    lint_ctx = LintContext("all")

    with mock.patch("planemo.shed_lint.tool_shed_instance", return_value=_FakeToolShed(ts_version)):
        lint_ctx.lint("shed_version", lint_shed_version, _ShedRepository(str(tmp_path)))

    assert _messages(lint_ctx) == expected
