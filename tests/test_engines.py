"""Unit tests for engines and runnables."""

import copy
import json
import os

import pytest

from planemo.engine import engine_context
from planemo.engine.galaxy import log_service_logs_on_failure
from planemo.engine.interface import (
    _absolutize_job_paths,
    materialized_job_paths,
)
from planemo.runnable import (
    cases,
    for_path,
    get_outputs,
    RunnableType,
)
from .test_utils import (
    create_test_context,
    TEST_DATA_DIR,
)

A_CWL_TOOL = os.path.join(TEST_DATA_DIR, "tools", "ok-cat1-tool.cwl")
A_CWL_WORKFLOW = os.path.join(TEST_DATA_DIR, "count-lines2-wf.cwl")

A_GALAXY_TOOL = os.path.join(TEST_DATA_DIR, "tools", "ok_select_param.xml")
A_GALAXY_GA_WORKFLOW = os.path.join(TEST_DATA_DIR, "test_workflow_1.ga")
A_GALAXY_YAML_WORKFLOW = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")

# workflows whose test definitions embed the job rather than pointing at a job file
A_COLLECTION_INPUT_WORKFLOW = os.path.join(TEST_DATA_DIR, "wf5-collection-input.gxwf.yml")
A_COMPOSITE_INPUT_WORKFLOW = os.path.join(TEST_DATA_DIR, "wf6-composite-inputs.gxwf.yml")
A_NESTED_COLLECTION_WORKFLOW = os.path.join(TEST_DATA_DIR, "wf8-collection-nested-input.gxwf.yml")
A_REMOTE_INPUT_WORKFLOW = os.path.join(TEST_DATA_DIR, "wf13_tool_shed_repository_gxformat2.yml")
# ... and one that points at tests/data/wf2-job.yml
A_FILE_JOB_WORKFLOW = os.path.join(TEST_DATA_DIR, "wf2.ga")

CAN_HANDLE = {
    "galaxy": {
        A_CWL_TOOL: True,
        A_CWL_WORKFLOW: True,
        A_GALAXY_TOOL: True,
        A_GALAXY_GA_WORKFLOW: True,
        A_GALAXY_YAML_WORKFLOW: True,
    },
    "cwltool": {
        A_CWL_TOOL: True,
        A_CWL_WORKFLOW: True,
        A_GALAXY_TOOL: False,
        A_GALAXY_GA_WORKFLOW: False,
        A_GALAXY_YAML_WORKFLOW: False,
    },
}


def test_can_handle():
    ctx = create_test_context()
    for engine_type in ["galaxy", "cwltool"]:
        with engine_context(ctx, engine=engine_type) as e:
            for key, value in CAN_HANDLE[engine_type].items():
                assert bool(e.can_run(for_path(key))) is value


def test_outputs():
    outputs = get_outputs(for_path(A_CWL_WORKFLOW))
    assert len(outputs) == 1
    output_id = outputs[0].get_id()
    assert output_id == "count_output"


def test_runnable_types():
    assert RunnableType.galaxy_tool.is_galaxy_artifact
    assert RunnableType.galaxy_workflow.is_galaxy_artifact
    assert not RunnableType.cwl_tool.is_galaxy_artifact
    assert not RunnableType.cwl_workflow.is_galaxy_artifact


class _RecordingContext:
    def __init__(self):
        self.messages = []

    def log(self, msg, *args):
        self.messages.append(msg)


class _ConfigWithServiceLogs:
    service_log_contents = {"celery.log": "Task galaxy.fetch_data raised unexpected"}


def test_service_logs_not_logged_when_tests_pass():
    ctx = _RecordingContext()
    log_service_logs_on_failure(ctx, _ConfigWithServiceLogs(), [{"data": {"status": "success"}}])
    assert ctx.messages == []


def test_service_logs_logged_when_test_fails():
    ctx = _RecordingContext()
    log_service_logs_on_failure(ctx, _ConfigWithServiceLogs(), [{"data": {"status": "error"}}])
    assert len(ctx.messages) == 1
    assert "celery.log" in ctx.messages[0]
    assert "raised unexpected" in ctx.messages[0]


def test_service_logs_logged_when_no_result_registered():
    """verify_tool blew up before registering anything - still want the logs."""
    ctx = _RecordingContext()
    log_service_logs_on_failure(ctx, _ConfigWithServiceLogs(), [])
    assert len(ctx.messages) == 1


def _materialized_jobs(workflow_path):
    """Run a real workflow's test cases through job materialization, return the job dicts."""
    jobs = []
    with materialized_job_paths(cases(for_path(workflow_path))) as job_paths:
        for job_path in job_paths:
            with open(job_path) as f:
                jobs.append(json.load(f))
    return jobs


def test_inline_jobs_are_not_written_to_the_test_directory():
    test_cases = cases(for_path(A_COLLECTION_INPUT_WORKFLOW))
    before = sorted(os.listdir(TEST_DATA_DIR))
    with materialized_job_paths(test_cases) as job_paths:
        assert len(job_paths) == len(test_cases)
        for job_path in job_paths:
            assert os.path.exists(job_path)
            assert os.path.dirname(job_path) != TEST_DATA_DIR
        assert sorted(os.listdir(TEST_DATA_DIR)) == before
    for job_path in job_paths:
        assert not os.path.exists(job_path)


def test_inline_job_directory_removed_when_the_engine_raises():
    with pytest.raises(RuntimeError, match="engine failed"):
        with materialized_job_paths(cases(for_path(A_COLLECTION_INPUT_WORKFLOW))) as job_paths:
            recorded = list(job_paths)
            raise RuntimeError("engine failed")
    for job_path in recorded:
        assert not os.path.exists(job_path)


def test_inline_job_relative_paths_resolve_to_real_test_data():
    collection_job, cwl_style_job = _materialized_jobs(A_COLLECTION_INPUT_WORKFLOW)
    hello = os.path.join(TEST_DATA_DIR, "hello.txt")
    assert collection_job["input1"]["elements"][0]["path"] == hello
    assert cwl_style_job["input1"][0]["path"] == hello
    assert os.path.exists(hello)


def test_nested_collection_paths_resolve_to_real_test_data():
    (job,) = _materialized_jobs(A_NESTED_COLLECTION_WORKFLOW)
    pair = job["input1"]["elements"][0]["elements"]
    assert [element["path"] for element in pair] == [os.path.join(TEST_DATA_DIR, "hello.txt")] * 2


def test_composite_data_paths_resolve_to_real_test_data():
    (job,) = _materialized_jobs(A_COMPOSITE_INPUT_WORKFLOW)
    paths = [item["path"] for item in job["input1"]["composite_data"]]
    assert paths == [
        os.path.join(TEST_DATA_DIR, "Example_Continuous.imzML"),
        os.path.join(TEST_DATA_DIR, "Example_Continuous.ibd"),
    ]
    assert all(os.path.exists(path) for path in paths)


def test_remote_inputs_are_left_alone():
    (job,) = _materialized_jobs(A_REMOTE_INPUT_WORKFLOW)
    pair = job["pe-fastq"]["elements"][0]["elements"]
    assert [element["location"] for element in pair] == [
        "https://github.com/GoekeLab/bioinformatics-workflows/raw/master/test_data/reads_1.fq.gz",
        "https://github.com/GoekeLab/bioinformatics-workflows/raw/master/test_data/reads_2.fq.gz",
    ]


def test_test_case_job_is_left_unmodified():
    test_cases = cases(for_path(A_COLLECTION_INPUT_WORKFLOW))
    jobs_before = copy.deepcopy([test_case.job for test_case in test_cases])
    with materialized_job_paths(test_cases):
        pass
    assert [test_case.job for test_case in test_cases] == jobs_before


def test_job_files_are_used_in_place_and_survive():
    test_cases = cases(for_path(A_FILE_JOB_WORKFLOW))
    job_file = os.path.join(TEST_DATA_DIR, "wf2-job.yml")
    with materialized_job_paths(test_cases) as job_paths:
        assert job_paths == [job_file]
    assert os.path.exists(job_file)


def test_paths_that_are_not_local_test_data_are_left_alone():
    job = {
        "absolute": {"class": "File", "path": "/already/absolute.txt"},
        # planemo test definitions do use URLs under "path" - see tests/data/cat_tool_url_job.json
        "url_under_path": {"class": "File", "path": "https://example.org/input.txt"},
        "cwl_reference": {"class": "File", "path": "#main/step/out"},
        "not_a_file": {"path": "this-is-just-a-parameter-value"},
    }
    assert _absolutize_job_paths(job, TEST_DATA_DIR) == job
