import os
import tempfile

from .test_utils import (
    CliTestCase,
    TEST_TOOLS_DIR,
)


UNFORMATTED_XML = """\
<tool id="test" name="Test" version="1.0">
<description>A test tool</description>
<command>echo hello</command>
<!-- A comment that should survive -->
<inputs>
<param name="input1" type="data" format="txt"
       label="Input file" help="Pick a file"/>
</inputs>
<outputs>
<data name="output1" format="txt"/>
</outputs>
<help><![CDATA[
Some help text.
]]></help>
</tool>
"""

FORMATTED_XML = """\
<tool id="test" name="Test" version="1.0">
    <description>A test tool</description>
    <command>echo hello</command>
    <!-- A comment that should survive -->
    <inputs>
        <param name="input1" type="data" format="txt" label="Input file" help="Pick a file"/>
    </inputs>
    <outputs>
        <data name="output1" format="txt"/>
    </outputs>
    <help><![CDATA[
Some help text.
]]></help>
</tool>
"""


class FormatTestCase(CliTestCase):
    def test_format_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, "tool.xml")
            with open(xml_path, "w") as f:
                f.write(UNFORMATTED_XML)
            self._check_exit_code(["format", xml_path])
            with open(xml_path) as f:
                result = f.read()
            assert result == FORMATTED_XML, f"Formatted output doesn't match expected:\n{result}"

    def test_format_dry_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, "tool.xml")
            with open(xml_path, "w") as f:
                f.write(UNFORMATTED_XML)
            self._check_exit_code(["format", "--dry-run", xml_path])
            with open(xml_path) as f:
                result = f.read()
            assert result == UNFORMATTED_XML, "Dry run should not modify the file"

    def test_format_idempotent(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, "tool.xml")
            with open(xml_path, "w") as f:
                f.write(FORMATTED_XML)
            self._check_exit_code(["format", xml_path])
            with open(xml_path) as f:
                result = f.read()
            assert result == FORMATTED_XML, "Already-formatted file should not change"

    def test_format_preserves_comments(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, "tool.xml")
            with open(xml_path, "w") as f:
                f.write(UNFORMATTED_XML)
            self._check_exit_code(["format", xml_path])
            with open(xml_path) as f:
                result = f.read()
            assert "<!-- A comment that should survive -->" in result

    def test_format_preserves_cdata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, "tool.xml")
            with open(xml_path, "w") as f:
                f.write(UNFORMATTED_XML)
            self._check_exit_code(["format", xml_path])
            with open(xml_path) as f:
                result = f.read()
            assert "<![CDATA[" in result
            assert "Some help text." in result

    def test_format_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for name in ["a.xml", "b.xml"]:
                with open(os.path.join(tmpdir, name), "w") as f:
                    f.write(UNFORMATTED_XML)
            # also a non-xml file that should be ignored
            with open(os.path.join(tmpdir, "readme.txt"), "w") as f:
                f.write("not xml")
            self._check_exit_code(["format", tmpdir])
            for name in ["a.xml", "b.xml"]:
                with open(os.path.join(tmpdir, name)) as f:
                    assert f.read() == FORMATTED_XML

    def test_format_skips_invalid_xml(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, "bad.xml")
            with open(xml_path, "w") as f:
                f.write("<tool><unclosed>")
            self._check_exit_code(["format", xml_path])
            with open(xml_path) as f:
                assert f.read() == "<tool><unclosed>"

    def test_format_existing_tool(self):
        """Format an existing test tool and verify it still parses."""
        ok_tool = os.path.join(TEST_TOOLS_DIR, "ok_conditional.xml")
        with tempfile.TemporaryDirectory() as tmpdir:
            copy_path = os.path.join(tmpdir, "ok_conditional.xml")
            with open(ok_tool) as src, open(copy_path, "w") as dst:
                dst.write(src.read())
            self._check_exit_code(["format", copy_path])
            # verify the result is valid XML
            from lxml import etree
            with open(copy_path) as f:
                etree.fromstring(f.read())
