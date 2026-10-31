import json
import os
from types import SimpleNamespace

from galaxy.util import parse_xml_string
from galaxy.util.unittest_utils import skip_if_github_down

from planemo import autoupdate
from planemo.autoupdate import get_newest_tool_id
from planemo.galaxy.config import get_shed_tools_conf_string_for_tool_ids
from planemo.galaxy.workflows import (
    get_tool_ids_for_workflow,
    get_toolshed_url_for_tool_id,
)
from .test_utils import (
    skip_if_testtoolshed_down,
    TEST_DATA_DIR,
)


def test_get_newest_tool_id():
    tool_ids = [
        "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy0",
        "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a",
    ]
    assert get_newest_tool_id(tool_ids) == tool_ids[0]


def test_get_newest_tool_id_galaxy_suffix_is_numeric():
    tool_ids = [
        "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy10",
        "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy9",
    ]
    assert get_newest_tool_id(tool_ids) == tool_ids[0]


@skip_if_github_down
def test_get_toolshed_url_for_tool_id():
    assert get_toolshed_url_for_tool_id("github.com/repos/iuc/tool/1") == "https://github.com"


@skip_if_testtoolshed_down
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


COREOGRAPH = "toolshed.g2.bx.psu.edu/repos/perssond/coreograph/unet_coreograph"


def _revision(numeric_revision, version, **kwds):
    return {"numeric_revision": numeric_revision, "tools": [{"guid": f"{COREOGRAPH}/{version}"}], **kwds}


def _find_updated(current, repos):
    return autoupdate.find_updated_tool_id(f"{COREOGRAPH}/{current}", repos)


def test_find_updated_tool_id_uses_latest_revision():
    # an older revision carries a version that sorts higher than the one in the latest revision
    repos = {
        "1:57f1260ca94e": _revision(1, "2.2.8.3"),
        "2:224e0cf4aaeb": _revision(2, "2.2.8+galaxy0"),
        "3:ee92746d141a": _revision(3, "2.2.8+galaxy1"),
    }
    assert _find_updated("2.2.8+galaxy1", repos) is None
    assert _find_updated("2.2.8+galaxy0", repos) == f"{COREOGRAPH}/2.2.8+galaxy1"


def test_outdated_tools(monkeypatch):
    repos = {"2:224e0cf4aaeb": _revision(2, "2.2.8+galaxy0"), "3:ee92746d141a": _revision(3, "2.2.8+galaxy1")}
    monkeypatch.setattr(autoupdate, "find_repos_from_tool_id", lambda tool_id, ts: ("", repos))
    wf_dict = {"steps": {"0": {"type": "tool", "tool_id": f"{COREOGRAPH}/2.2.8+galaxy0"}}}
    assert autoupdate.outdated_tools(SimpleNamespace(log=print), wf_dict, []) == {
        COREOGRAPH: {"current": f"{COREOGRAPH}/2.2.8+galaxy0", "updated": f"{COREOGRAPH}/2.2.8+galaxy1"}
    }


def test_find_updated_tool_id_skips_uninstallable_revision():
    repos = {
        "2:224e0cf4aaeb": _revision(2, "2.2.8+galaxy0"),
        "3:ee92746d141a": _revision(3, "2.2.8+galaxy1", downloadable=False),
    }
    assert _find_updated("2.2.8+galaxy0", repos) is None


def _wf(*versions, subworkflow_versions=()):
    steps = {str(i): {"type": "tool", "tool_id": f"{COREOGRAPH}/{v}"} for i, v in enumerate(versions)}
    if subworkflow_versions:
        steps["sub"] = {"type": "subworkflow", "subworkflow": _wf(*subworkflow_versions)}
    return {"steps": steps}


def test_updated_tools():
    assert autoupdate.updated_tools(_wf("2.2.8+galaxy0"), _wf("2.2.8+galaxy0")) == {}
    assert autoupdate.updated_tools(_wf("2.2.8+galaxy0"), _wf("2.2.8+galaxy1")) == {
        COREOGRAPH: {"current": f"{COREOGRAPH}/2.2.8+galaxy0", "updated": f"{COREOGRAPH}/2.2.8+galaxy1"}
    }
    assert autoupdate.updated_tools(
        _wf(subworkflow_versions=["2.2.8+galaxy0"]), _wf(subworkflow_versions=["2.2.8+galaxy1"])
    ) == {COREOGRAPH: {"current": f"{COREOGRAPH}/2.2.8+galaxy0", "updated": f"{COREOGRAPH}/2.2.8+galaxy1"}}


def test_updated_tools_to_version_used_by_another_step():
    assert autoupdate.updated_tools(_wf("2.2.8+galaxy0", "2.2.8+galaxy1"), _wf("2.2.8+galaxy1", "2.2.8+galaxy1")) == {
        COREOGRAPH: {"current": f"{COREOGRAPH}/2.2.8+galaxy0", "updated": f"{COREOGRAPH}/2.2.8+galaxy1"}
    }
    assert autoupdate.updated_tools(_wf("2.2.8+galaxy0", "2.2.8+galaxy1"), _wf("2.2.8+galaxy0", "2.2.8+galaxy1")) == {}


def test_find_updated_tool_id_never_downgrades():
    repos = {"1:57f1260ca94e": _revision(1, "2.2.8.3"), "3:ee92746d141a": _revision(3, "2.2.8+galaxy1")}
    assert _find_updated("2.2.8.3", repos) is None


def test_find_updated_tool_id_skips_malicious_revision():
    repos = {
        "2:224e0cf4aaeb": _revision(2, "2.2.8+galaxy0"),
        "3:ee92746d141a": _revision(3, "2.2.8+galaxy1", malicious=True),
    }
    assert _find_updated("2.2.8+galaxy0", repos) is None


def test_find_updated_tool_id_orders_revisions_numerically():
    repos = {"10:ee92746d141a": _revision(10, "2.2.8+galaxy1"), "9:224e0cf4aaeb": _revision(9, "2.2.8+galaxy0")}
    assert _find_updated("2.2.8+galaxy0", repos) == f"{COREOGRAPH}/2.2.8+galaxy1"


def test_updated_tools_from_gxformat2_original():
    original = {"class": "GalaxyWorkflow", "steps": [{"tool_id": f"{COREOGRAPH}/2.2.8+galaxy0"}]}
    assert autoupdate.updated_tools(original, _wf("2.2.8+galaxy1")) == {
        COREOGRAPH: {"current": f"{COREOGRAPH}/2.2.8+galaxy0", "updated": f"{COREOGRAPH}/2.2.8+galaxy1"}
    }
