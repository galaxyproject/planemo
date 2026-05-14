"""Module describing the planemo ``invocation_download`` command."""

import json
import os

import click
from bioblend.galaxy import GalaxyInstance

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.galaxy.activity import invocation_to_run_response
from planemo.io import info
from planemo.output_models import PlanemoInvocationDownloadManifest
from planemo.runnable import get_outputs
from planemo.runnable_resolve import for_runnable_identifier


def invocation_download_manifest(
    invocation_id,
    output_directory,
    outputs_dict,
    output_json=None,
    path_type="relative",
    missing_output_reasons=None,
):
    """Build the machine-readable manifest for ``invocation_download``."""
    output_directory = os.path.abspath(output_directory)
    missing_output_reasons = missing_output_reasons or {}
    outputs = {}
    missing_outputs = []
    for output_id, output_value in outputs_dict.items():
        if output_value is None:
            missing_outputs.append({"id": output_id, "reason": missing_output_reasons.get(output_id, "missing")})
        else:
            outputs[output_id] = _manifest_output_value(output_value, output_directory, path_type)

    manifest = {
        "invocation_id": invocation_id,
        "output_directory": _manifest_path(output_directory, output_directory, path_type),
        "path_type": path_type,
        "outputs": outputs,
        "missing_outputs": missing_outputs,
    }
    if output_json is not None:
        manifest["output_json"] = _manifest_path(output_json, output_directory, path_type)
    return manifest


def _manifest_output_value(value, output_directory, path_type):
    if isinstance(value, dict):
        output_value = {}
        for key, item in value.items():
            if key == "path":
                output_value[key] = _manifest_path(item, output_directory, path_type)
            else:
                output_value[key] = _manifest_output_value(item, output_directory, path_type)
        return output_value
    elif isinstance(value, list):
        return [_manifest_output_value(item, output_directory, path_type) for item in value]
    else:
        return value


def _manifest_path(path, output_directory, path_type):
    if path_type == "absolute":
        return os.path.abspath(path)
    else:
        return os.path.relpath(path, output_directory)


def _missing_output_reasons(runnable, gi):
    reasons = {}
    for runnable_output in get_outputs(runnable, gi=gi):
        output_id = runnable_output.get_id()
        if output_id:
            reasons[output_id] = "skipped" if runnable_output.is_optional() else "missing"
    return reasons


@click.command("invocation_download")
@options.run_output_directory_option()
@options.run_output_json_option()
@click.option(
    "--output_json_path_type",
    type=click.Choice(["relative", "absolute"]),
    default="relative",
    show_default=True,
    help="Path style to use in the output JSON manifest.",
)
@options.profile_option()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@click.argument("invocation_id", type=click.STRING)
@click.option("--ignore_missing_output/--no_ignore_missing_output", default=True, help="Ignore missing output files")
@command_function
def cli(ctx, invocation_id, output_directory, output_json_path_type, ignore_missing_output, **kwds):
    """Download output files from a completed Galaxy workflow invocation.

    This command allows downloading output files after a workflow has been executed
    through Galaxy's web interface or through planemo.

    \b
        % planemo invocation_download INVOCATION_ID
    """

    profile = kwds.get("profile")
    if profile is not None:
        profile = profiles.ensure_profile(ctx, profile)
        key = profile["galaxy_admin_key"] or profile["galaxy_user_key"]
        url = profile["galaxy_url"]
    else:
        url = kwds.get("galaxy_url")
        key = kwds.get("galaxy_admin_key") or kwds.get("galaxy_user_key")

    gi = GalaxyInstance(url=url, key=key)

    invocation_data = gi.invocations.show_invocation(invocation_id)
    workflow_id = gi.workflows.show_workflow(workflow_id=invocation_data["workflow_id"], instance=True)["id"]
    runnable = for_runnable_identifier(ctx, workflow_id, kwds)
    run_response = invocation_to_run_response(ctx, gi, runnable, invocation_data)
    if output_directory is None:
        output_directory = f"output_{invocation_id}"
        os.makedirs(output_directory, exist_ok=True)
        info(f"No output_directory provided, saving to {output_directory}")
    run_response.collect_outputs(output_directory, ignore_missing_output=ignore_missing_output)
    output_json = kwds.get("output_json")
    if output_json:
        manifest = invocation_download_manifest(
            invocation_id,
            output_directory,
            run_response.outputs_dict,
            output_json=output_json,
            path_type=output_json_path_type,
            missing_output_reasons=_missing_output_reasons(runnable, gi),
        )
        manifest = PlanemoInvocationDownloadManifest.model_validate(manifest).model_dump(mode="json")
        with open(output_json, "w") as f:
            json.dump(manifest, f, ensure_ascii=False)
