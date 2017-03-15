"""Module provides generic interface to running Galaxy tools and workflows."""

import json
import os
import tempfile

from bioblend.galaxy.client import Client
from galaxy.tools.parser import get_tool_source
from six import iteritems

from planemo.io import wait_on
from planemo.runnable import (
    ErrorRunResponse,
    get_outputs,
    RunnableType,
    SuccessfulRunResponse,
)

DEFAULT_HISTORY_NAME = "CWL Target History"
ERR_NO_SUCH_TOOL = ("Failed to find tool with ID [%s] in Galaxy - cannot execute job. "
                    "You may need to enable verbose logging and determine why the tool did not load. [%s]")


def execute(config, runnable, job_path, **kwds):
    """Execute a Galaxy activity."""
    try:
        return _execute(config, runnable, job_path, **kwds)
    except Exception as e:
        return ErrorRunResponse(str(e))


def _execute(config, runnable, job_path, **kwds):
    user_gi = config.user_gi
    admin_gi = config.gi

    history_id = _history_id(user_gi, **kwds)

    galaxy_paths, job_dict, datasets = stage_in(config, user_gi, history_id, job_path, **kwds)

    if runnable.type in [RunnableType.galaxy_tool, RunnableType.cwl_tool]:
        response_class = GalaxyToolRunResponse
        tool_id = _tool_id(runnable.path)
        if runnable.type == RunnableType.cwl_tool:
            inputs_representation = "cwl"
        else:
            inputs_representation = "galaxy"
        try:
            user_gi.tools.show_tool(tool_id)
        except Exception as e:
            raise Exception(ERR_NO_SUCH_TOOL % (tool_id, e))
        run_tool_payload = dict(
            history_id=history_id,
            tool_id=tool_id,
            inputs=job_dict,
            inputs_representation=inputs_representation,
        )
        tool_run_response = user_gi.tools._tool_post(run_tool_payload)

        job = tool_run_response["jobs"][0]
        job_id = job["id"]
        final_state = _wait_for_job(user_gi, job_id)
        if final_state != "ok":
            msg = "Failed to run CWL job final job state is [%s]." % final_state
            with open("errored_galaxy.log", "w") as f:
                f.write(config.log_contents)
            raise Exception(msg)

        job_info = admin_gi.jobs.show_job(job_id)
        response_kwds = {
            'job_info': job_info,
            'api_run_response': tool_run_response,
        }

    elif runnable.type in [RunnableType.galaxy_workflow]:
        response_class = GalaxyWorkflowRunResponse
        workflow_id = config.workflow_id(runnable.path)
        # TODO: update bioblend to allow inputs_by.
        # invocation = user_gi.worklfows.invoke_workflow(
        #    workflow_id,
        #    history_id=history_id,
        #    inputs=job_dict,
        # )
        payload = dict(
            workflow_id=workflow_id,
            history_id=history_id,
            inputs=job_dict,
            inputs_by="name",
            allow_tool_state_corrections=True,
        )
        invocations_url = "%s/%s/invocations" % (
            user_gi._make_url(user_gi.workflows),
            workflow_id,
        )
        invocation = Client._post(user_gi.workflows, payload, url=invocations_url)
        invocation_id = invocation["id"]
        _wait_for_invocation(user_gi, workflow_id, invocation_id)
        final_state = _wait_for_history(user_gi, history_id)
        if final_state != "ok":
            msg = "Failed to run CWL job final job state is [%s]." % final_state
            with open("errored_galaxy.log", "w") as f:
                f.write(config.log_contents)
            raise Exception(msg)
        response_kwds = {
            'workflow_id': workflow_id,
            'invocation_id': invocation_id,
        }
    else:
        raise NotImplementedError()

    run_response = response_class(
        runnable=runnable,
        user_gi=user_gi,
        history_id=history_id,
        galaxy_paths=galaxy_paths,
        log=config.log_contents,
        **response_kwds
    )
    output_directory = kwds.get("output_directory", None)
    run_response.collect_outputs(output_directory)
    return run_response


def stage_in(config, user_gi, history_id, job_path, **kwds):

    # Figure out what "." should be here instead.
    def upload(file_path):
        return user_gi.tools.upload_file(file_path, history_id)

    with open(job_path, "r") as f:
        job = json.load(f)

    job_dir = os.path.dirname(job_path)
    job_dict, datasets = galactic_job_json(job, job_dir, upload)

    if datasets:
        final_state = _wait_for_history(user_gi, history_id)
    else:
        # Mark uploads as ok because nothing to do.
        final_state = "ok"

    if final_state != "ok":
        msg = "Failed to run CWL job final job state is [%s]." % final_state
        with open("errored_galaxy.log", "w") as f:
            f.write(config.log_contents)
        raise Exception(msg)

    galaxy_paths = []
    for (dataset, local_path) in datasets:
        dataset_full = user_gi.datasets.show_dataset(dataset["id"])
        galaxy_path = dataset_full["file_name"]
        job_path = os.path.join(job_dir, local_path)
        galaxy_paths.append((job_path, galaxy_path))

    return galaxy_paths, job_dict, datasets


class GalaxyBaseRunResponse(SuccessfulRunResponse):

    def __init__(
        self,
        runnable,
        user_gi,
        history_id,
        galaxy_paths,
        log,
    ):
        self._runnable = runnable
        self._user_gi = user_gi
        self._history_id = history_id
        self._log = log

        self._job_info = None

        self.galaxy_paths = galaxy_paths

        self._outputs_dict = None

    def collect_outputs(self, output_directory):
        assert self._outputs_dict is None

        outputs_dict = {}
        if not output_directory:
            # TODO: rather than creating a directory just use
            # Galaxy paths if they are available in this
            # configuration.
            output_directory = tempfile.mkdtemp()

        for output in get_outputs(self._runnable):
            output_id = output.get_id()
            dataset = self.get_dataset_metadata(output)
            destination = self.download_output_to(output, output_directory)

            if dataset["file_ext"] == "expression.json":
                with open(destination, "r") as f:
                    dict_value = json.load(f)
            else:
                dict_value = {
                    "path": destination,
                    "class": "File",
                }
            outputs_dict[output_id] = dict_value
        self._outputs_dict = outputs_dict

    @property
    def log(self):
        return self._log

    @property
    def job_info(self):
        if self._job_info is not None:
            return dict(
                stdout=self._job_info["stdout"],
                stderr=self._job_info["stderr"],
                command_line=self._job_info["command_line"],
            )
        return None

    @property
    def outputs_dict(self):
        return self._outputs_dict

    def get_dataset_metadata(self, output):
        output_dataset_id = self.output_dataset_id(output)
        return self._user_gi.histories.show_dataset(
            self._history_id,
            output_dataset_id,
        )

    def download_output_to(self, output, output_directory):
        output_id = output.get_id()
        output_dataset_id = self.output_dataset_id(output)
        destination = os.path.join(output_directory, output_id)
        self._user_gi.histories.download_dataset(
            self._history_id,
            output_dataset_id,
            file_path=destination,
            use_default_filename=False,
        )
        return destination


class GalaxyToolRunResponse(GalaxyBaseRunResponse):

    def __init__(
        self,
        runnable,
        user_gi,
        history_id,
        galaxy_paths,
        log,
        job_info,
        api_run_response,
    ):
        super(GalaxyToolRunResponse, self).__init__(
            runnable=runnable,
            user_gi=user_gi,
            history_id=history_id,
            galaxy_paths=galaxy_paths,
            log=log,
        )
        self._job_info = job_info
        self.api_run_response = api_run_response

    def output_dataset_id(self, output):
        outputs = self.api_run_response["outputs"]
        output_id = output.get_id()
        output_dataset_id = None
        for output in outputs:
            if output["output_name"] == output_id:
                output_dataset_id = output["id"]

        return output_dataset_id


class GalaxyWorkflowRunResponse(GalaxyBaseRunResponse):

    def __init__(
        self,
        runnable,
        user_gi,
        history_id,
        galaxy_paths,
        log,
        workflow_id,
        invocation_id,
    ):
        super(GalaxyWorkflowRunResponse, self).__init__(
            runnable=runnable,
            user_gi=user_gi,
            history_id=history_id,
            galaxy_paths=galaxy_paths,
            log=log,
        )
        self._workflow_id = workflow_id
        self._invocation_id = invocation_id

    def output_dataset_id(self, output):
        invocation = self._user_gi.workflows.show_invocation(
            self._workflow_id,
            self._invocation_id,
        )
        workflow_output = output.workflow_output
        order_index = workflow_output.order_index

        invocation_steps = invocation["steps"]
        output_steps = [s for s in invocation_steps if s["order_index"] == order_index]
        assert len(output_steps) == 1, "More than one step matching outputs, behavior undefined."
        output_step = output_steps[0]
        job_id = output_step["job_id"]
        assert job_id, "Output doesn't define a job_id, behavior undefined."
        job_info = self._user_gi.jobs.show_job(job_id, full_details=True)
        job_outputs = job_info["outputs"]
        output_name = workflow_output.output_name
        assert output_name in job_outputs, "No output [%s] found for output job."
        job_output = job_outputs[output_name]
        assert "id" in job_output, "Job output [%s] does not contain 'id'." % job_output
        return job_output["id"]


def _tool_id(tool_path):
    tool_source = get_tool_source(tool_path)
    return tool_source.parse_id()


def _history_id(gi, **kwds):
    history_id = kwds.get("history_id", None)
    if history_id is None:
        history_name = kwds.get("history_name", DEFAULT_HISTORY_NAME)
        history_id = gi.histories.create_history(history_name)["id"]
    return history_id


def _wait_for_invocation(gi, workflow_id, invocation_id):
    def state_func():
        return gi.workflows.show_invocation(workflow_id, invocation_id)

    return _wait_on_state(state_func)


def _wait_for_history(gi, history_id):
    def state_func():
        return gi.histories.show_history(history_id)

    return _wait_on_state(state_func)


def _wait_for_job(gi, job_id):
    def state_func():
        return gi.jobs.show_job(job_id, full_details=True)

    return _wait_on_state(state_func)


def _wait_on_state(state_func):

    def get_state():
        response = state_func()
        state = response["state"]
        if str(state) not in ["running", "queued", "new", "ready"]:
            return state
        else:
            return None

    final_state = wait_on(get_state, "state", timeout=100)
    return final_state


# Now this is the newer version of this function as of 4/29
def galactic_job_json(job, test_data_directory, upload_func):
    datasets = []

    def upload(file_path):
        if not os.path.isabs(file_path):
            file_path = os.path.join(test_data_directory, file_path)
        _ensure_file_exists(file_path)
        return upload_func(file_path)

    def replacement_item(value):
        if not isinstance(value, dict):
            return value

        type_class = value.get("class", None)
        if type_class != "File":
            return value

        file_path = value.get("location", None)
        if file_path is None:
            return value

        upload_response = upload(file_path)
        dataset = upload_response["outputs"][0]
        datasets.append((dataset, file_path))
        dataset_id = dataset["id"]
        return {"src": "hda", "id": dataset_id}

    replace_keys = {}
    for key, value in iteritems(job):
        if isinstance(value, dict):
            replace_keys[key] = replacement_item(value)
        elif isinstance(value, list):
            new_list = []
            for item in value:
                new_list.append(replacement_item(item))
            replace_keys[key] = new_list

    job.update(replace_keys)
    return job, datasets


def _ensure_file_exists(file_path):
    if not os.path.exists(file_path):
        template = "File [%s] does not exist - parent directory [%s] does %sexist, cwd is [%s]"
        parent_directory = os.path.dirname(file_path)
        message = template % (
            file_path,
            parent_directory,
            "" if os.path.exists(parent_directory) else "not ",
            os.getcwd(),
        )
        raise Exception(message)


__all__ = (
    "execute",
)
