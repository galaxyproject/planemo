import os

from planemo.engine.test import (
    test_runnables as t_runnables
)
from planemo.runnable import (
    for_paths,
)
from .test_utils import (
    mark,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    skip_unless_module,
    target_galaxy_branch,
    test_context,
    TEST_DATA_DIR,
)


@skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
@skip_unless_module("toil")
def test_toil_tests():
    ctx = test_context()
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
    ctx = test_context()
    random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
    cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
    test_artifact = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
    runnables = for_paths([test_artifact])
    kwds = {
        "engine": "galaxy",
        "no_dependency_resolution": True,
        "paste_test_data_paths": False,
        "galaxy_branch": target_galaxy_branch(),
        "extra_tools": [random_lines, cat]
    }
    exit_code = t_runnables(
        ctx,
        runnables,
        **kwds
    )
    assert exit_code == 0
