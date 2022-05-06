import json
import os
from tempfile import NamedTemporaryFile

from planemo.engine.test import test_runnables as t_runnables
from planemo.runnable import for_paths
from .test_utils import (
    create_test_context,
    mark,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    skip_unless_module,
    target_galaxy_branch,
    TEST_DATA_DIR,
)


@skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
@skip_unless_module("toil")
def test_toil_tests():
    ctx = create_test_context()
    test_artifact = os.path.join(TEST_DATA_DIR, "int_tool.cwl")
    runnables = for_paths([test_artifact])
    exit_code = t_runnables(
        ctx,
        runnables,
        engine="toil",
    )
    assert exit_code == 0


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
@mark.tests_galaxy_branch
def test_galaxy_wf_tests():
    ctx = create_test_context()
    random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
    cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
    test_artifact = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
    runnables = for_paths([test_artifact])
    kwds = {
        "engine": "galaxy",
        "no_dependency_resolution": True,
        "paste_test_data_paths": False,
        "galaxy_branch": target_galaxy_branch(),
        "extra_tools": [random_lines, cat],
    }
    exit_code = t_runnables(ctx, runnables, **kwds)
    assert exit_code == 0


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
@mark.tests_galaxy_branch
def test_galaxy_workflow_collection_output():
    ctx = create_test_context()
    test_artifact = os.path.join(TEST_DATA_DIR, "wf7-collection-output.gxwf.yml")
    collection_creates_pair = os.path.join(TEST_DATA_DIR, "collection_creates_pair_2.xml")
    runnables = for_paths([test_artifact])
    kwds = {
        "engine": "galaxy",
        "no_dependency_resolution": True,
        "paste_test_data_paths": False,
        "galaxy_branch": target_galaxy_branch(),
        "extra_tools": [collection_creates_pair],
    }
    exit_code = t_runnables(ctx, runnables, **kwds)
    assert exit_code == 0


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
@mark.tests_galaxy_branch
def test_galaxy_workflow_collection_output_fail():
    ctx = create_test_context()
    test_artifact = os.path.join(TEST_DATA_DIR, "wf7-collection-output-fail.gxwf.yml")
    collection_creates_pair = os.path.join(TEST_DATA_DIR, "collection_creates_pair_2.xml")
    runnables = for_paths([test_artifact])
    kwds = {
        "engine": "galaxy",
        "no_dependency_resolution": True,
        "paste_test_data_paths": False,
        "galaxy_branch": target_galaxy_branch(),
        "extra_tools": [collection_creates_pair],
    }
    exit_code = t_runnables(ctx, runnables, **kwds)
    assert exit_code == 1


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
# @mark.tests_galaxy_branch only works >= 20.09 or newer
def test_galaxy_workflow_tags():
    with NamedTemporaryFile(prefix="data_manager_test_json") as json_out:
        ctx = create_test_context()
        test_artifact = os.path.join(TEST_DATA_DIR, "wf10-tags-and-rules.gxwf.yml")
        collection_cat_list = os.path.join(TEST_DATA_DIR, "cat_list.xml")
        runnables = for_paths([test_artifact])
        kwds = {
            "engine": "galaxy",
            "no_dependency_resolution": True,
            "paste_test_data_paths": False,
            "galaxy_branch": "dev",
            "extra_tools": [collection_cat_list],
            "test_output_json": json_out.name,
        }
        try:
            exit_code = t_runnables(ctx, runnables, **kwds)
            assert exit_code == 0
        except Exception:
            with open(json_out.name) as f:
                print(f.read())
            raise


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
@mark.tests_galaxy_branch
def test_galaxy_workflow_nested_collection_inputs():
    ctx = create_test_context()
    test_artifact = os.path.join(TEST_DATA_DIR, "wf8-collection-nested-input.gxwf.yml")
    collection_cat_pair = os.path.join(TEST_DATA_DIR, "cat_pair.xml")
    collection_cat_list = os.path.join(TEST_DATA_DIR, "cat_list.xml")
    runnables = for_paths([test_artifact])
    kwds = {
        "engine": "galaxy",
        "no_dependency_resolution": True,
        "paste_test_data_paths": False,
        "galaxy_branch": target_galaxy_branch(),
        "extra_tools": [collection_cat_pair, collection_cat_list],
    }
    exit_code = t_runnables(ctx, runnables, **kwds)
    assert exit_code == 0


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
@mark.tests_galaxy_branch
def test_galaxy_workflow_non_data_inputs():
    ctx = create_test_context()
    test_artifact = os.path.join(TEST_DATA_DIR, "wf9-int-input.gxwf.yml")
    runnables = for_paths([test_artifact])
    kwds = {
        "engine": "galaxy",
        "no_dependency_resolution": True,
        "paste_test_data_paths": False,
        "galaxy_branch": target_galaxy_branch(),
    }
    exit_code = t_runnables(ctx, runnables, **kwds)
    assert exit_code == 0


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
@mark.tests_galaxy_branch
def test_galaxy_workflow_step_failed():
    ctx = create_test_context()
    test_artifact = os.path.join(TEST_DATA_DIR, "wf_failed_step.ga")
    runnables = for_paths([test_artifact])
    with NamedTemporaryFile(prefix="result_json") as json_out:
        kwds = {
            "engine": "galaxy",
            "no_dependency_resolution": True,
            "paste_test_data_paths": False,
            "extra_tools": ["$GALAXY_FUNCTIONAL_TEST_TOOLS"],
            "test_output_json": json_out.name,
            "galaxy_branch": target_galaxy_branch(),
        }
        exit_code = t_runnables(ctx, runnables, **kwds)
        assert exit_code == 1
        report = json.load(json_out)
    data = report["tests"][0]["data"]
    assert data["status"] == "error"
    assert data["execution_problem"]
    invocation_steps = data["invocation_details"]["steps"]
    assert len(invocation_steps) == 2
    first_step, second_step = invocation_steps.values()
    assert first_step["state"] == "scheduled"
    job = first_step["jobs"][0]
    assert job["exit_code"] == 127
    assert job["state"] == "error"
    assert second_step["state"] == "scheduled"
    assert second_step["jobs"][0]["state"] == "paused"
