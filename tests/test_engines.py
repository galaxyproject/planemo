"""Unit tests for engines and runnables."""

import os

from planemo.engine import engine_context
from planemo.runnable import for_path
from planemo.runnable import get_outputs
from .test_utils import CWL_DRAFT3_DIR, test_context, TEST_DATA_DIR

A_CWL_TOOL = os.path.join(CWL_DRAFT3_DIR, "cat3-tool.cwl")
A_CWL_WORKFLOW = os.path.join(CWL_DRAFT3_DIR, "count-lines1-wf.cwl")

A_GALAXY_TOOL = os.path.join(TEST_DATA_DIR, "tools", "ok_select_param.xml")
A_GALAXY_GA_WORKFLOW = os.path.join(TEST_DATA_DIR, "test_workflow_1.ga")
A_GALAXY_YAML_WORKFLOW = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")

CAN_HANDLE = {
    'galaxy': {
        A_CWL_TOOL: True,
        A_CWL_WORKFLOW: True,
        A_GALAXY_TOOL: True,
        A_GALAXY_GA_WORKFLOW: True,
        A_GALAXY_YAML_WORKFLOW: True,
    },
    'cwltool': {
        A_CWL_TOOL: True,
        A_CWL_WORKFLOW: True,
        A_GALAXY_TOOL: False,
        A_GALAXY_GA_WORKFLOW: False,
        A_GALAXY_YAML_WORKFLOW: False,
    }
}


def test_can_handle():
    ctx = test_context()
    for engine_type in ["galaxy", "cwltool"]:
        with engine_context(ctx, engine=engine_type) as e:
            for key, value in CAN_HANDLE[engine_type].items():
                assert bool(e.can_run(for_path(key))) is value


def test_outputs():
    outputs = get_outputs(for_path(A_CWL_TOOL))
    assert len(outputs) == 1
    output_id = outputs[0].get_id()
    assert output_id == "output_file"
