"""Training:tutorial functions."""
import os
import shutil

import pytest

from planemo.engine import engine_context
from planemo.engine.galaxy import LocalManagedGalaxyEngine
from planemo.training import Training
from planemo.training.topic import Topic
from planemo.training.tutorial import (
    get_galaxy_datatype,
    get_hands_on_boxes_from_local_galaxy,
    get_hands_on_boxes_from_running_galaxy,
    get_wf_inputs,
    get_wf_param_values,
    get_zenodo_record,
    Tutorial,
)
from planemo.training.utils import save_to_yaml
from .test_training import (
    create_existing_tutorial,
    CTX,
    datatype_fp,
    KWDS,
    RUNNABLE,
    tuto_fp,
    wf,
    WF_FP,
    wf_param_values,
    zenodo_link,
)
from .test_utils import skip_if_environ

topic = Topic()
training = Training(KWDS)


def test_get_galaxy_datatype() -> None:
    """Test :func:`planemo.training.tutorial.get_galaxy_datatype`."""
    assert get_galaxy_datatype("csv", datatype_fp) == "csv"
    assert get_galaxy_datatype("test", datatype_fp) == "strange_datatype"
    assert "# Please add" in get_galaxy_datatype("unknown", datatype_fp)


def test_get_zenodo_record() -> None:
    """Test :func:`planemo.training.tutorial.get_zenodo_record`."""
    z_record, req_res = get_zenodo_record(zenodo_link)
    file_link_prefix = "https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f"
    assert z_record == "1321885"
    assert "files" in req_res
    assert req_res["files"][0]["type"] in ["rdata", "csv"]
    assert file_link_prefix in req_res["files"][0]["links"]["self"]
    # check with wrong zenodo link
    z_record, req_res = get_zenodo_record("https://zenodo.org/api/records/zenodooo")
    assert z_record is None
    assert "files" in req_res
    assert len(req_res["files"]) == 0
    # using DOI
    z_link = "https://doi.org/10.5281/zenodo.1321885"
    z_record, req_res = get_zenodo_record(z_link)
    file_link_prefix = "https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f"
    assert z_record == "1321885"
    assert "files" in req_res
    assert req_res["files"][0]["type"] in ["rdata", "csv"]
    assert file_link_prefix in req_res["files"][0]["links"]["self"]


def test_get_wf_inputs() -> None:
    """Test :func:`planemo.training.tutorial.get_wf_inputs`."""
    step_inp = {
        "tables_1|table": {"output_name": "output", "id": 2},
        "add_to_database|withdb": {"output_name": "output", "id": 0},
        "tables_0|table": {"output_name": "output", "id": 1},
        "add_to_database|tab_0|tt": {"output_name": "output", "id": 0},
        "tables_2|section|sect": {"output_name": "output", "id": 1},
        "tables_3|tables_0|sect": {"output_name": "output", "id": 1},
    }
    step_inputs = get_wf_inputs(step_inp)
    assert "tables" in step_inputs
    assert "0" in step_inputs["tables"]
    assert "table" in step_inputs["tables"]["0"]
    assert "2" in step_inputs["tables"]
    assert "section" in step_inputs["tables"]["2"]
    assert "sect" in step_inputs["tables"]["2"]["section"]
    assert "output_name" in step_inputs["tables"]["2"]["section"]["sect"]
    assert "add_to_database" in step_inputs
    assert "withdb" in step_inputs["add_to_database"]
    assert "tab" in step_inputs["add_to_database"]
    assert "0" in step_inputs["add_to_database"]["tab"]
    assert "tt" in step_inputs["add_to_database"]["tab"]["0"]


def test_get_wf_param_values() -> None:
    """Test :func:`planemo.training.tutorial.get_wf_param_values`."""
    wf_step = wf["steps"]["3"]
    wf_param_value_tests = get_wf_param_values(wf_step["tool_state"], get_wf_inputs(wf_step["input_connections"]))
    assert isinstance(wf_param_value_tests, dict)
    for k in wf_param_values:
        assert k in wf_param_value_tests


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_get_hands_on_boxes_from_local_galaxy() -> None:
    """Test :func:`planemo.training.tutorial.get_hands_on_boxes_from_local_galaxy`."""
    tuto_body = get_hands_on_boxes_from_local_galaxy(KWDS, WF_FP, CTX)
    assert_body_contains(tuto_body, "## Sub-step with **FastQC**")
    assert_body_contains(tuto_body, "## Sub-step with **Query Tabular**")
    assert_body_contains(tuto_body, "## Sub-step with **Select first**")


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_get_hands_on_boxes_from_running_galaxy() -> None:
    """Test :func:`planemo.training.tutorial.get_hands_on_boxes_from_running_galaxy`."""
    galaxy_url = f"http://{KWDS['host']}:{KWDS['port']}"
    with engine_context(CTX, **KWDS) as galaxy_engine:
        assert isinstance(galaxy_engine, LocalManagedGalaxyEngine)
        with galaxy_engine.ensure_runnables_served([RUNNABLE]) as config:
            wf_id = config.workflow_id(WF_FP)
            tuto_body = get_hands_on_boxes_from_running_galaxy(wf_id, galaxy_url, config.user_api_key)

    assert_body_contains(tuto_body, "## Sub-step with **FastQC**")
    assert_body_contains(tuto_body, "## Sub-step with **Query Tabular**")
    assert_body_contains(tuto_body, "## Sub-step with **Select first**")


def test_tutorial_init() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.init`."""
    # with default parameter
    tuto = Tutorial(training=training, topic=topic)
    assert tuto.name == "new_tuto"
    assert tuto.title == "The new tutorial"
    assert tuto.zenodo_link == ""
    assert tuto.hands_on
    assert not tuto.slides
    assert tuto.init_wf_id is None
    assert tuto.init_wf_fp is None
    assert tuto.datatype_fp == ""
    assert "new_tuto" in tuto.dir
    assert "## Sub-step with **My Tool**" in tuto.body
    assert tuto.data_lib
    # with non default parameter
    tuto = Tutorial(training=training, topic=topic, name="my_tuto", title="The tutorial", zenodo_link="URL")
    assert tuto.name == "my_tuto"
    assert tuto.title == "The tutorial"
    assert tuto.zenodo_link == "URL"
    assert "my_tuto" in tuto.dir


def test_tutorial_init_from_kwds() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.init_from_kwds`."""
    kwds = {
        "tutorial_name": "my_tuto",
        "tutorial_title": "Title of tuto",
        "hands_on": True,
        "slides": True,
        "workflow": WF_FP,
        "workflow_id": "id",
        "zenodo_link": None,
        "datatypes": datatype_fp,
    }
    tuto = Tutorial(training=training, topic=topic)
    tuto.init_from_kwds(kwds)
    assert tuto.name == "my_tuto"
    assert tuto.title == "Title of tuto"
    assert tuto.zenodo_link == ""
    assert "Which biological questions are addressed by the tutorial?" in tuto.questions
    assert tuto.hands_on
    assert tuto.slides
    assert tuto.init_wf_id == "id"
    assert tuto.init_wf_fp == WF_FP
    assert tuto.datatype_fp == datatype_fp
    assert "my_tuto" in tuto.dir


def test_tutorial_init_from_existing_tutorial() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.init_from_existing_tutorial`."""
    tuto = Tutorial(training=training, topic=topic)
    # non existing tutorial
    with pytest.raises(Exception, match="The tutorial existing_tutorial does not exists. It should be created"):
        tuto.init_from_existing_tutorial("existing_tutorial")
    # existing tutorial
    create_existing_tutorial("existing_tutorial", tuto_fp, tuto.topic)
    tuto.init_from_existing_tutorial("existing_tutorial")
    assert tuto.title == "A tutorial to test"
    assert "A learning objective" in tuto.objectives
    assert tuto.time_estimation == "1H"
    assert "the_best_contributor" in tuto.contributors
    assert "# First section" in tuto.body
    shutil.rmtree("topics")


def test_tutorial_init_data_lib() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.init_data_lib`."""
    tuto = Tutorial(training=training, topic=topic)
    tuto.init_data_lib()
    assert tuto.data_lib["destination"]["type"] == "library"
    assert tuto.data_lib["items"][0]["name"] == topic.title
    assert tuto.data_lib["items"][0]["items"][0]["name"] == tuto.title
    # from existing data library file
    os.makedirs(tuto.dir)
    tuto.data_lib = {}
    tuto.init_data_lib()
    assert tuto.data_lib["items"][0]["name"] == topic.title
    assert tuto.data_lib["items"][0]["items"][0]["name"] == tuto.title
    # other tutorial already there and add the new one
    tuto.data_lib["items"][0]["items"][0]["name"] = "Different tutorial"
    save_to_yaml(tuto.data_lib, tuto.data_lib_fp)
    tuto.init_data_lib()
    assert tuto.data_lib["items"][0]["items"][0]["name"] == "Different tutorial"
    assert tuto.data_lib["items"][0]["items"][1]["name"] == tuto.title
    shutil.rmtree("topics")


def test_tutorial_get_tuto_metata() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.get_tuto_metata`."""
    tuto = Tutorial(training=training, topic=topic)
    tuto.questions = ["q1", "q2"]
    metadata = tuto.get_tuto_metata()
    assert "title: The new tutorial" in metadata
    assert "- q1" in metadata


def test_tutorial_set_dir_name() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.set_dir_name`."""
    tuto = Tutorial(training=training, topic=topic)
    tuto.name = "the_tuto"
    tuto.set_dir_name()
    assert tuto.name in tuto.dir
    assert tuto.name in tuto.tuto_fp
    assert tuto.name in tuto.slide_fp
    assert tuto.name in tuto.data_lib_fp
    assert tuto.name in tuto.wf_dir
    assert tuto.name in tuto.wf_fp


def test_tutorial_exists() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.exists`."""
    # default
    tuto = Tutorial(training=training, topic=topic)
    assert not tuto.exists()
    # after dir creation
    os.makedirs(tuto.dir)
    assert tuto.exists()
    shutil.rmtree("topics")


def test_tutorial_has_workflow() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.has_workflow`."""
    # default
    tuto = Tutorial(training=training, topic=topic)
    assert not tuto.has_workflow()
    # with wf filepath
    tuto.init_wf_fp = WF_FP
    assert tuto.has_workflow()
    # with no wf filepah nor wf id
    tuto.init_wf_fp = None
    tuto.init_wf_id = ""
    assert not tuto.has_workflow()
    # with wf id
    tuto.init_wf_id = "ID"
    assert tuto.has_workflow()


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_tutorial_export_workflow_file() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.export_workflow_file`."""
    tuto = Tutorial(training=training, topic=topic)
    os.makedirs(tuto.wf_dir)
    # with worflow fp
    tuto.init_wf_fp = WF_FP
    tuto.export_workflow_file()
    assert os.path.exists(tuto.wf_fp)
    # with workflow id
    tuto.init_wf_fp = None
    os.remove(tuto.wf_fp)
    galaxy_url = f"http://{KWDS['host']}:{KWDS['port']}"
    with engine_context(CTX, **KWDS) as galaxy_engine:
        assert isinstance(galaxy_engine, LocalManagedGalaxyEngine)
        with galaxy_engine.ensure_runnables_served([RUNNABLE]) as config:
            tuto.init_wf_id = config.workflow_id(WF_FP)
            tuto.training.galaxy_url = galaxy_url
            tuto.training.galaxy_api_key = config.user_api_key
            tuto.export_workflow_file()
    assert os.path.exists(tuto.wf_fp)
    shutil.rmtree("topics")


def test_tutorial_get_files_from_zenodo() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.get_files_from_zenodo`."""
    tuto = Tutorial(training=training, topic=topic, zenodo_link=zenodo_link)
    tuto.datatype_fp = datatype_fp
    files, z_record = tuto.get_files_from_zenodo()
    assert z_record == "1321885"
    # test links
    file_link_prefix = "https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f"
    assert file_link_prefix in tuto.zenodo_file_links[0]
    # test files dict
    assert file_link_prefix in files[0]["url"]
    assert files[0]["src"] == "url"
    assert files[0]["info"] == zenodo_link
    assert "# Please add" in files[0]["ext"]
    assert files[1]["ext"] == "csv"


def test_tutorial_prepare_data_library_from_zenodo() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.prepare_data_library_from_zenodo`."""
    # without zenodo link
    tuto = Tutorial(training=training, topic=topic)
    tuto.datatype_fp = datatype_fp
    os.makedirs(tuto.wf_dir)
    tuto.prepare_data_library_from_zenodo()
    assert os.path.exists(tuto.data_lib_fp)
    with open(tuto.data_lib_fp) as fh:
        assert "DOI" not in fh.read()
    # with zenodo link
    tuto.zenodo_link = zenodo_link
    tuto.prepare_data_library_from_zenodo()
    with open(tuto.data_lib_fp) as fh:
        assert "DOI: 10.5281/zenodo" in fh.read()
    shutil.rmtree("topics")


def test_tutorial_write_hands_on_tutorial() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.write_hands_on_tutorial`."""
    tuto = Tutorial(training=training, topic=topic)
    os.makedirs(tuto.wf_dir)
    tuto.zenodo_file_links = ["URL1", "URL2"]
    tuto.write_hands_on_tutorial()
    assert os.path.exists(tuto.tuto_fp)
    with open(tuto.tuto_fp) as tuto_f:
        tuto_c = tuto_f.read()
        assert "layout: tutorial_hands_on" in tuto_c
        assert "# Introduction" in tuto_c
        assert "URL1" in tuto_c
        assert "# Conclusion" in tuto_c
    shutil.rmtree("topics")


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_tutorial_create_hands_on_tutorial() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.create_hands_on_tutorial`."""
    tuto = Tutorial(training=training, topic=topic)
    os.makedirs(tuto.wf_dir)
    # with init_wf_id and no Galaxy URL
    tuto.init_wf_id = "ID"
    tuto.training.galaxy_url = None
    with pytest.raises(Exception, match="No Galaxy URL given"):
        tuto.create_hands_on_tutorial(CTX)
    # with init_wf_id and no Galaxy API key
    tuto.init_wf_id = "ID"
    tuto.training.galaxy_url = f"http://{KWDS['host']}:{KWDS['port']}"
    tuto.training.galaxy_api_key = None
    with pytest.raises(Exception, match="No API key to access the given Galaxy instance"):
        tuto.create_hands_on_tutorial(CTX)
    # with init_wf_id
    with engine_context(CTX, **KWDS) as galaxy_engine:
        assert isinstance(galaxy_engine, LocalManagedGalaxyEngine)
        with galaxy_engine.ensure_runnables_served([RUNNABLE]) as config:
            tuto.init_wf_id = config.workflow_id(WF_FP)
            tuto.training.galaxy_api_key = config.user_api_key
            tuto.create_hands_on_tutorial(CTX)
    assert os.path.exists(tuto.tuto_fp)
    os.remove(tuto.tuto_fp)
    # with init_wf_fp
    tuto.init_wf_id = None
    tuto.init_wf_fp = WF_FP
    tuto.create_hands_on_tutorial(CTX)
    assert os.path.exists(tuto.tuto_fp)
    shutil.rmtree("topics")


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_tutorial_create_tutorial() -> None:
    """Test :func:`planemo.training.tutorial.tutorial.create_tutorial`."""
    tuto = Tutorial(training=training, topic=topic)
    tuto.init_from_kwds(
        {
            "tutorial_name": "my_tuto",
            "tutorial_title": "Title of tuto",
            "hands_on": True,
            "slides": True,
            "workflow": WF_FP,
            "workflow_id": None,
            "zenodo_link": zenodo_link,
            "datatypes": datatype_fp,
        }
    )
    tuto.create_tutorial(CTX)
    assert os.path.exists(tuto.dir)
    assert os.path.exists(tuto.tour_dir)
    assert os.path.exists(tuto.wf_dir)
    assert os.path.exists(tuto.data_lib_fp)
    assert os.path.exists(tuto.tuto_fp)
    assert os.path.exists(tuto.slide_fp)
    with open(tuto.slide_fp) as fh:
        assert "layout: tutorial_slides" in fh.read()
    shutil.rmtree("topics")


def assert_body_contains(body: str, contents: str) -> None:
    if contents not in body:
        message = f"Expected to find contents [{contents}] in body [{body}]"
        raise AssertionError(message)
