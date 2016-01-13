""" High-level client sitting on top of bioblend for running CWL
stuff in Galaxy.
"""
from __future__ import print_function

import json
import os

from six import iteritems

from planemo.io import wait_on


DEFAULT_HISTORY_NAME = "CWL Target History"


def run_cwl_tool(tool_path, job_path, config, **kwds):
    user_gi = config.user_gi
    admin_gi = config.gi

    tool_id = _tool_id(tool_path)
    history_id = _history_id(user_gi, **kwds)
    job_dict = _galactic_job_json(job_path, user_gi, history_id)
    final_state = _wait_for_history(user_gi, history_id)
    if final_state != "ok":
        msg = "Failed to run CWL job final job state is [%s]." % final_state
        with open("errored_galaxy.log", "w") as f:
            f.write(config.log_contents)
        raise Exception(msg)
    run_tool_payload = dict(
        history_id=history_id,
        tool_id=tool_id,
        tool_inputs=job_dict,
        inputs_representation="cwl",
    )
    run_response = user_gi.tools._tool_post(run_tool_payload)
    job = run_response["jobs"][0]
    job_id = job["id"]
    final_state = _wait_for_job(user_gi, job_id)
    if final_state != "ok":
        msg = "Failed to run CWL job final job state is [%s]." % final_state
        with open("errored_galaxy.log", "w") as f:
            f.write(config.log_contents)
        raise Exception(msg)

    job_info = admin_gi.jobs.show_job(job_id)
    cwl_command_state = job_info["cwl_command_state"]

    return CwlRunResponse(cwl_command_state, run_response)


class CwlRunResponse(object):

    def __init__(self, cwl_command_state, api_run_response):
        self.cwl_command_state = cwl_command_state
        self.api_run_response = api_run_response


def _history_id(gi, **kwds):
    history_id = kwds.get("history_id", None)
    if history_id is None:
        history_name = kwds.get("history_name", DEFAULT_HISTORY_NAME)
        history_id = gi.histories.create_history(history_name)["id"]
    return history_id


def _tool_id(tool_path):
    tool_id, _ = os.path.splitext(os.path.basename(tool_path))
    return tool_id


def _galactic_job_json(job_path, gi, history_id):
    with open(job_path, "r") as f:
        job_as_dict = json.load(f)

    replace_keys = {}
    for key, value in iteritems(job_as_dict):
        if isinstance(value, dict):
            type_class = value.get("class", None)
            if type_class != "File":
                continue

            file_path = value.get("path", None)
            if file_path is None:
                continue

            if not os.path.isabs(file_path):
                directory = os.path.dirname(job_path)
                file_path = os.path.join(directory, file_path)

            upload_response = gi.tools.upload_file(file_path, history_id)
            dataset_id = upload_response["outputs"][0]["id"]

            replace_keys[key] = {"src": "hda", "id": dataset_id}

    job_as_dict.update(replace_keys)
    return job_as_dict


def _wait_for_history(gi, history_id):
    def state_func():
        return gi.histories.show_history(history_id)

    return _wait_on_state(state_func)


def _wait_for_job(gi, job_id):
    def state_func():
        return gi.jobs.show_job(job_id)

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
