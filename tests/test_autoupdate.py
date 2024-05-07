import json
import os

from galaxy.util import parse_xml_string

from planemo.autoupdate import get_newest_tool_id
from planemo.galaxy.config import get_shed_tools_conf_string_for_tool_ids
from planemo.galaxy.workflows import (
    get_tool_ids_for_workflow,
    get_toolshed_url_for_tool_id,
)
from .test_utils import TEST_DATA_DIR


def test_get_newest_tool_id():
    tool_ids = [
        "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy0",
        "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a",
    ]
    assert get_newest_tool_id(tool_ids) == tool_ids[0]


def test_get_toolshed_url_for_tool_id():
    assert get_toolshed_url_for_tool_id("github.com/repos/iuc/tool/1") == "https://github.com"


def test_get_shed_tools_conf_string_for_workflow():
    test_artifact = os.path.join(TEST_DATA_DIR, "wf_same_name_outputs.ga")
    with open(test_artifact) as wf:
        wf_dict = json.load(wf)
        wf_dict["steps"]["3"]["tool_id"] = wf_dict["steps"]["3"]["tool_id"].replace("toolshed", "testtoolshed")
    tool_ids = get_tool_ids_for_workflow(wf_dict=wf_dict)
    xml_string = get_shed_tools_conf_string_for_tool_ids(tool_ids)
    element = parse_xml_string(xml_string)
    tool_shed_entries = element.findall("tool_shed")
    assert len(tool_shed_entries) == 2
    assert tool_shed_entries[0].attrib["url"] == "https://toolshed.g2.bx.psu.edu"
    assert tool_shed_entries[0].attrib["name"] == "toolshed.g2.bx.psu.edu"
    assert tool_shed_entries[1].attrib["url"] == "https://testtoolshed.g2.bx.psu.edu"
    assert tool_shed_entries[1].attrib["name"] == "testtoolshed.g2.bx.psu.edu"
