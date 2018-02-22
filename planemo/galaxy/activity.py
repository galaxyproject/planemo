"""Module provides generic interface to running Galaxy tools and workflows."""

import json
import os
import tempfile
import time

import bioblend
import requests
import yaml
from bioblend.galaxy.client import Client
from galaxy.tools.cwl.util import (
    DirectoryUploadTarget,
    FileUploadTarget,
    galactic_job_json,
    invocation_to_output,
    output_properties,
    output_to_cwl_json,
    tool_response_to_output,
)
from galaxy.tools.parser import get_tool_source
from galaxy.util import safe_makedirs

from planemo.galaxy.api import summarize_history
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


def execute(ctx, config, runnable, job_path, **kwds):
    """Execute a Galaxy activity."""
    try:
        return _execute(ctx, config, runnable, job_path, **kwds)
    except Exception as e:
        return ErrorRunResponse(str(e))


def _verified_tool_id(runnable, user_gi):
    tool_id = _tool_id(runnable.path)
    try:
        user_gi.tools.show_tool(tool_id)
    except Exception as e:
        raise Exception(ERR_NO_SUCH_TOOL % (tool_id, e))
    return tool_id


def _inputs_representation(runnable):
    if runnable.type == RunnableType.cwl_tool:
        inputs_representation = "cwl"
    else:
        inputs_representation = "galaxy"
    return inputs_representation


def _execute(ctx, config, runnable, job_path, **kwds):
    user_gi = config.user_gi
    admin_gi = config.gi

    history_id = _history_id(user_gi, **kwds)

    galaxy_paths, job_dict, _ = stage_in(ctx, runnable, config, user_gi, history_id, job_path, **kwds)

    if runnable.type in [RunnableType.galaxy_tool, RunnableType.cwl_tool]:
        response_class = GalaxyToolRunResponse
        tool_id = _verified_tool_id(runnable, user_gi)
        inputs_representation = _inputs_representation(runnable)
        run_tool_payload = dict(
            history_id=history_id,
            tool_id=tool_id,
            inputs=job_dict,
            inputs_representation=inputs_representation,
        )
        ctx.vlog("Post to Galaxy tool API with payload [%s]" % run_tool_payload)
        tool_run_response = user_gi.tools._tool_post(run_tool_payload)

        job = tool_run_response["jobs"][0]
        job_id = job["id"]
        try:
            final_state = _wait_for_job(user_gi, job_id)
        except Exception:
            summarize_history(ctx, user_gi, history_id)
            raise
        if final_state != "ok":
            msg = "Failed to run CWL tool job final job state is [%s]." % final_state
            summarize_history(ctx, user_gi, history_id)
            with open("errored_galaxy.log", "w") as f:
                f.write(config.log_contents)
            raise Exception(msg)

        ctx.vlog("Final job state was ok, fetching details for job [%s]" % job_id)
        job_info = admin_gi.jobs.show_job(job_id)
        response_kwds = {
            'job_info': job_info,
            'api_run_response': tool_run_response,
        }
        if ctx.verbose:
            summarize_history(ctx, user_gi, history_id)
    elif runnable.type in [RunnableType.galaxy_workflow, RunnableType.cwl_workflow]:
        response_class = GalaxyWorkflowRunResponse
        workflow_id = config.workflow_id(runnable.path)
        ctx.vlog("Found Galaxy workflow ID [%s] for path [%s]" % (workflow_id, runnable.path))
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
        ctx.vlog("Waiting for invocation [%s]" % invocation_id)
        try:
            final_invocation_state = _wait_for_invocation(ctx, user_gi, history_id, workflow_id, invocation_id)
        except Exception:
            ctx.vlog("Problem waiting on invocation...")
            summarize_history(ctx, user_gi, history_id)
            raise
        ctx.vlog("Final invocation state is [%s]" % final_invocation_state)
        final_state = _wait_for_history(ctx, user_gi, history_id)
        if final_state != "ok":
            msg = "Failed to run workflow final history state is [%s]." % final_state
            summarize_history(ctx, user_gi, history_id)
            with open("errored_galaxy.log", "w") as f:
                f.write(config.log_contents)
            raise Exception(msg)
        ctx.vlog("Final history state is 'ok'")
        response_kwds = {
            'workflow_id': workflow_id,
            'invocation_id': invocation_id,
        }
    else:
        raise NotImplementedError()

    run_response = response_class(
        ctx=ctx,
        runnable=runnable,
        user_gi=user_gi,
        history_id=history_id,
        galaxy_paths=galaxy_paths,
        log=config.log_contents,
        **response_kwds
    )
    output_directory = kwds.get("output_directory", None)
    ctx.vlog("collecting outputs from run...")
    run_response.collect_outputs(ctx, output_directory)
    ctx.vlog("collecting outputs complete")
    return run_response


def stage_in(ctx, runnable, config, user_gi, history_id, job_path, **kwds):

    def upload_func(upload_target):
        if isinstance(upload_target, FileUploadTarget):
            file_path = upload_target.path
            upload_payload = user_gi.tools._upload_payload(
                history_id,
                file_type="auto",
            )
            name = os.path.basename(file_path)
            upload_payload["inputs"]["files_0|auto_decompress"] = False
            upload_payload["inputs"]["auto_decompress"] = False
            upload_payload["inputs"]["files_0|url_paste"] = "file://%s" % os.path.abspath(file_path)
            upload_payload["inputs"]["files_0|NAME"] = name
            if upload_target.secondary_files:
                upload_payload["inputs"]["files_1|url_paste"] = "file://%s" % os.path.abspath(upload_target.secondary_files)
                upload_payload["inputs"]["files_1|type"] = "upload_dataset"
                upload_payload["inputs"]["files_1|auto_decompress"] = True
                upload_payload["inputs"]["file_count"] = "2"
                upload_payload["inputs"]["force_composite"] = "True"

            ctx.vlog("upload_payload is %s" % upload_payload)
            return user_gi.tools._tool_post(upload_payload, files_attached=False)
        elif isinstance(upload_target, DirectoryUploadTarget):
            tar_path = upload_target.tar_path

            upload_payload = user_gi.tools._upload_payload(
                history_id,
                file_type="tar",
            )
            upload_payload["inputs"]["files_0|auto_decompress"] = False
            upload_payload["inputs"]["files_0|url_paste"] = "file://%s" % tar_path
            tar_upload_response = user_gi.tools._tool_post(upload_payload, files_attached=False)
            convert_response = user_gi.tools.run_tool(
                tool_id="CONVERTER_tar_to_directory",
                tool_inputs={"input1": {"src": "hda", "id": tar_upload_response["outputs"][0]["id"]}},
                history_id=history_id,
            )
            assert "outputs" in convert_response, convert_response
            return convert_response
        else:
            content = json.dumps(upload_target.object)
            return user_gi.tools.paste_content(
                content,
                history_id,
                file_type="expression.json",
            )

    def create_collection_func(element_identifiers, collection_type):
        payload = {
            "name": "dataset collection",
            "instance_type": "history",
            "history_id": history_id,
            "element_identifiers": element_identifiers,
            "collection_type": collection_type,
            "fields": None if collection_type != "record" else "auto",
        }
        dataset_collections_url = user_gi.url + "/dataset_collections"
        dataset_collection = Client._post(user_gi.histories, payload, url=dataset_collections_url)
        return dataset_collection

    with open(job_path, "r") as f:
        job = yaml.load(f)

    # Figure out what "." should be here instead.
    job_dir = os.path.dirname(job_path)
    job_dict, datasets = galactic_job_json(
        job,
        job_dir,
        upload_func,
        create_collection_func,
        tool_or_workflow="tool" if runnable.type in [RunnableType.cwl_tool, RunnableType.galaxy_tool] else "workflow",
    )

    if datasets:
        final_state = _wait_for_history(ctx, user_gi, history_id)

        for (dataset, path) in datasets:
            dataset_details = user_gi.histories.show_dataset(
                history_id,
                dataset["id"],
            )
            ctx.vlog("Uploaded dataset for path [%s] with metadata [%s]" % (path, dataset_details))
    else:
        # Mark uploads as ok because nothing to do.
        final_state = "ok"

    ctx.vlog("final state is %s" % final_state)
    if final_state != "ok":
        msg = "Failed to run job final job state is [%s]." % final_state
        summarize_history(ctx, user_gi, history_id)
        with open("errored_galaxy.log", "w") as f:
            f.write(config.log_contents)
        raise Exception(msg)

    galaxy_paths = []
    for (dataset, upload_target) in datasets:
        if isinstance(upload_target, FileUploadTarget):
            local_path = upload_target.path
            ctx.vlog("fetching full dataset for %s, %s" % (dataset, local_path))
            dataset_full = user_gi.datasets.show_dataset(dataset["id"])
            galaxy_path = dataset_full["file_name"]
            ctx.vlog("galaxy_path is %s" % galaxy_path)
            job_path = os.path.join(job_dir, local_path)
            galaxy_paths.append((job_path, galaxy_path))

    ctx.vlog("galaxy_paths are %s" % galaxy_paths)
    return galaxy_paths, job_dict, datasets


class GalaxyBaseRunResponse(SuccessfulRunResponse):

    def __init__(
        self,
        ctx,
        runnable,
        user_gi,
        history_id,
        galaxy_paths,
        log,
    ):
        self._ctx = ctx
        self._runnable = runnable
        self._user_gi = user_gi
        self._history_id = history_id
        self._log = log

        self._job_info = None

        self.galaxy_paths = galaxy_paths

        self._outputs_dict = None

    def to_galaxy_output(self, output):
        """Convert runnable output to a GalaxyOutput object.

        Subclasses for workflow and tool execution override this.
        """
        raise NotImplementedError()

    def _get_extra_files(self, dataset_details):
        extra_files_url = "%s/%s/contents/%s/extra_files" % (
            self._user_gi._make_url(self._user_gi.histories), self._history_id, dataset_details["id"]
        )
        extra_files = Client._get(self._user_gi.jobs, url=extra_files_url)
        return extra_files

    def _get_metadata(self, history_content_type, content_id):
        if history_content_type == "dataset":
            return self._user_gi.histories.show_dataset(
                self._history_id,
                content_id,
            )
        elif history_content_type == "dataset_collection":
            return self._user_gi.histories.show_dataset_collection(
                self._history_id,
                content_id,
            )
        else:
            raise Exception("Unknown history content type encountered [%s]" % history_content_type)

    def collect_outputs(self, ctx, output_directory):
        assert self._outputs_dict is None, "collect_outputs pre-condition violated"

        outputs_dict = {}
        if not output_directory:
            # TODO: rather than creating a directory just use
            # Galaxy paths if they are available in this
            # configuration.
            output_directory = tempfile.mkdtemp()

        def get_dataset(dataset_details, filename=None):
            parent_basename = dataset_details.get("cwl_file_name")
            if not parent_basename:
                parent_basename = dataset_details.get("name")
            file_ext = dataset_details["file_ext"]
            if file_ext == "directory":
                # TODO: rename output_directory to outputs_directory because we can have output directories
                # and this is confusing...
                the_output_directory = os.path.join(output_directory, parent_basename)
                safe_makedirs(the_output_directory)
                destination = self.download_output_to(dataset_details, the_output_directory, filename=filename)
            else:
                destination = self.download_output_to(dataset_details, output_directory, filename=filename)
            if filename is None:
                basename = parent_basename
            else:
                basename = os.path.basename(filename)

            return {"path": destination, "basename": basename}

        ctx.vlog("collecting outputs to directory %s" % output_directory)
        for runnable_output in get_outputs(self._runnable):
            output_id = runnable_output.get_id()
            output_dict_value = None
            if self._runnable.type in [RunnableType.cwl_workflow, RunnableType.cwl_tool]:
                galaxy_output = self.to_galaxy_output(runnable_output)
                cwl_output = output_to_cwl_json(
                    galaxy_output,
                    self._get_metadata,
                    get_dataset,
                    self._get_extra_files,
                    pseduo_location=True,
                )
                output_dict_value = cwl_output
            else:
                # TODO: deprecate this route for finding workflow outputs,
                # it is a brittle and bad approach...
                output_dataset_id = self.output_dataset_id(runnable_output)
                dataset = self._get_metadata("dataset", output_dataset_id)
                dataset_dict = get_dataset(dataset)
                ctx.vlog("populated destination [%s]" % dataset_dict["path"])

                if dataset["file_ext"] == "expression.json":
                    with open(dataset_dict["path"], "r") as f:
                        output_dict_value = json.load(f)
                else:
                    output_dict_value = output_properties(**dataset_dict)

            outputs_dict[output_id] = output_dict_value

        self._outputs_dict = outputs_dict
        ctx.vlog("collected outputs [%s]" % self._outputs_dict)

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

    def download_output_to(self, dataset_details, output_directory, filename=None):
        if filename is None:
            local_filename = dataset_details.get("cwl_file_name") or dataset_details.get("name")
        else:
            local_filename = filename
        destination = os.path.join(output_directory, local_filename)
        self._history_content_download(
            self._history_id,
            dataset_details["id"],
            to_path=destination,
            filename=filename,
        )
        return destination

    def _history_content_download(self, history_id, dataset_id, to_path, filename=None):
        user_gi = self._user_gi
        url = user_gi.url + "/histories/%s/contents/%s/display" % (history_id, dataset_id)

        data = {}
        if filename:
            data["filename"] = filename

        r = requests.get(url, params=data, verify=user_gi.verify, stream=True, timeout=user_gi.timeout)
        r.raise_for_status()

        with open(to_path, 'wb') as fp:
            for chunk in r.iter_content(chunk_size=bioblend.CHUNK_SIZE):
                if chunk:
                    fp.write(chunk)


class GalaxyToolRunResponse(GalaxyBaseRunResponse):

    def __init__(
        self,
        ctx,
        runnable,
        user_gi,
        history_id,
        galaxy_paths,
        log,
        job_info,
        api_run_response,
    ):
        super(GalaxyToolRunResponse, self).__init__(
            ctx=ctx,
            runnable=runnable,
            user_gi=user_gi,
            history_id=history_id,
            galaxy_paths=galaxy_paths,
            log=log,
        )
        self._job_info = job_info
        self.api_run_response = api_run_response

    def is_collection(self, output):
        # TODO: Make this more rigorous - search both output and output
        # collections - throw an exception if not found in either place instead
        # of just assuming all non-datasets are collections.
        return self.output_dataset_id(output) is None

    def to_galaxy_output(self, runnable_output):
        output_id = runnable_output.get_id()
        return tool_response_to_output(self.api_run_response, self._history_id, output_id)

    def output_dataset_id(self, output):
        outputs = self.api_run_response["outputs"]
        output_id = output.get_id()
        output_dataset_id = None
        self._ctx.vlog("Looking for id [%s] in outputs [%s]" % (output_id, outputs))
        for output in outputs:
            if output["output_name"] == output_id:
                output_dataset_id = output["id"]

        return output_dataset_id


class GalaxyWorkflowRunResponse(GalaxyBaseRunResponse):

    def __init__(
        self,
        ctx,
        runnable,
        user_gi,
        history_id,
        galaxy_paths,
        log,
        workflow_id,
        invocation_id,
    ):
        super(GalaxyWorkflowRunResponse, self).__init__(
            ctx=ctx,
            runnable=runnable,
            user_gi=user_gi,
            history_id=history_id,
            galaxy_paths=galaxy_paths,
            log=log,
        )
        self._workflow_id = workflow_id
        self._invocation_id = invocation_id

    def to_galaxy_output(self, runnable_output):
        output_id = runnable_output.get_id()
        self._ctx.vlog("checking for output in invocation [%s]" % self._invocation)
        return invocation_to_output(self._invocation, self._history_id, output_id)

    def output_dataset_id(self, output):
        invocation = self._invocation
        if "outputs" in invocation:
            # Use newer workflow outputs API.

            output_name = output.get_id()
            if output_name in invocation["outputs"]:
                return invocation["outputs"][output.get_id()]["id"]
            else:
                raise Exception("Failed to find output [%s] in invocation outputs [%s]" % (output_name, invocation["outputs"]))
        else:
            # Assume the output knows its order_index and such - older line of
            # development not worth persuing.
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

    @property
    def _invocation(self):
        invocation = self._user_gi.workflows.show_invocation(
            self._workflow_id,
            self._invocation_id,
        )
        return invocation


def _tool_id(tool_path):
    tool_source = get_tool_source(tool_path)
    return tool_source.parse_id()


def _history_id(gi, **kwds):
    history_id = kwds.get("history_id", None)
    if history_id is None:
        history_name = kwds.get("history_name", DEFAULT_HISTORY_NAME)
        history_id = gi.histories.create_history(history_name)["id"]
    return history_id


def _wait_for_invocation(ctx, gi, history_id, workflow_id, invocation_id):

    def state_func():
        if _retry_on_timeouts(ctx, gi, lambda gi: has_jobs_in_states(gi, history_id, ["error", "deleted", "deleted_new"])):
            raise Exception("Problem running workflow, one or more jobs failed.")

        return _retry_on_timeouts(ctx, gi, lambda gi: gi.workflows.show_invocation(workflow_id, invocation_id))

    return _wait_on_state(state_func)


def _retry_on_timeouts(ctx, gi, f):
    gi.timeout = 60
    try_count = 5
    try:
        for try_num in range(try_count):
            start_time = time.time()
            try:
                return f(gi)
            except Exception:
                end_time = time.time()
                if end_time - start_time > 45 and (try_num + 1) < try_count:
                    ctx.vlog("Galaxy seems to have timedout, retrying to fetch status.")
                    continue
                else:
                    raise
    finally:
        gi.timeout = None


def has_jobs_in_states(gi, history_id, states):
    params = {"history_id": history_id}
    jobs_url = gi._make_url(gi.jobs)
    jobs = Client._get(gi.jobs, params=params, url=jobs_url)

    target_jobs = [j for j in jobs if j["state"] in states]

    return len(target_jobs) > 0


def _wait_for_history(ctx, gi, history_id):

    def has_active_jobs(gi):
        if has_jobs_in_states(gi, history_id, ["new", "upload", "waiting", "queued", "running"]):
            return True
        else:
            return None

    wait_on(lambda: _retry_on_timeouts(ctx, gi, has_active_jobs), "active jobs", timeout=60 * 60 * 24)

    def state_func():
        return _retry_on_timeouts(ctx, gi, lambda gi: gi.histories.show_history(history_id))

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

    final_state = wait_on(get_state, "state", timeout=60 * 60 * 24)
    return final_state


__all__ = (
    "execute",
)
