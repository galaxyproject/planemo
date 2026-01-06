"""Utilities for Galaxy workflows."""

import json
import os
import tempfile
from collections import namedtuple
from functools import lru_cache
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Tuple,
    TYPE_CHECKING,
)
from urllib.parse import (
    quote,
    urlparse,
)

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

import requests
import yaml
from ephemeris import (
    generate_tool_list_from_ga_workflow_files,
    shed_tools,
)
from gxformat2.converter import python_to_workflow
from gxformat2.interface import (
    BioBlendImporterGalaxyInterface,
    ImporterGalaxyInterface,
)
from gxformat2.normalize import (
    inputs_normalized,
    outputs_normalized,
)

from planemo.galaxy.api import (
    get_dict_from_workflow,
    gi,
)
from planemo.io import warn

FAILED_REPOSITORIES_MESSAGE = "Failed to install one or more repositories."
GALAXY_WORKFLOWS_PREFIX = "gxid://workflows/"
GALAXY_WORKFLOW_INSTANCE_PREFIX = "gxid://workflow-instance/"
TRS_WORKFLOWS_PREFIX = "trs://"
MAIN_TOOLSHED_URL = "https://toolshed.g2.bx.psu.edu"


def parse_trs_id(trs_id: str) -> Optional[Dict[str, str]]:
    """Parse a TRS ID into a full TRS URL.

    Args:
        trs_id: TRS ID in format: [#]workflow/github.com/org/repo/workflow_name[/version]
                Examples:
                - workflow/github.com/org/repo/main
                - #workflow/github.com/org/repo/main/v0.1.14
                - workflow/github.com/iwc-workflows/parallel-accession-download/main

    Returns:
        Dict with key 'trs_url' containing the full TRS API URL,
        or None if invalid
    """
    # Remove leading # if present
    if trs_id.startswith("#"):
        trs_id = trs_id[1:]

    # Expected format: workflow/github.com/org/repo[/workflow_name][/version]
    # Some workflows use the repo name as the workflow name (4 parts for tool ID)
    # Others have a separate workflow name (5 parts for tool ID)
    parts = trs_id.split("/")
    if len(parts) < 4:
        return None

    artifact_type = parts[0]  # workflow or tool
    service = parts[1]  # github.com
    owner = parts[2]
    repo = parts[3]

    # Determine if we have a workflow name and/or version
    # Format could be:
    #   workflow/github.com/org/repo (4 parts) - no workflow name, no version
    #   workflow/github.com/org/repo/workflow_name (5 parts) - with workflow name, no version
    #   workflow/github.com/org/repo/workflow_name/version (6 parts) - with both
    #   workflow/github.com/org/repo/workflow_name/versions/version (7 parts) - full URL format
    if len(parts) == 4:
        workflow_name = None
        version = None
    elif len(parts) == 5:
        # 5th part is the workflow name (e.g., "main" in most cases)
        workflow_name = parts[4]
        version = None
    elif len(parts) >= 6:
        # 6+ parts: 5th is workflow name, 6th might be "versions" keyword or version
        workflow_name = parts[4]
        # Check if this is full URL format with "versions" keyword
        if len(parts) >= 7 and parts[5] == "versions":
            # Full URL format: .../workflow_name/versions/version
            version = parts[6]
        else:
            # Short format: .../workflow_name/version
            version = parts[5]
    else:
        workflow_name = None
        version = None

    # Build the TRS tool ID
    # Format: #workflow/github.com/org/repo[/workflow_name]
    if workflow_name:
        trs_tool_id = f"#{artifact_type}/{service}/{owner}/{repo}/{workflow_name}"
    else:
        trs_tool_id = f"#{artifact_type}/{service}/{owner}/{repo}"
    # URL-encode the tool ID for API calls
    encoded_tool_id = quote(trs_tool_id, safe="")

    # Build the full TRS URL
    # Dockstore is the primary TRS server for GitHub workflows
    trs_base_url = "https://dockstore.org/api/ga4gh/trs/v2/tools/"

    if version:
        # Specific version requested
        trs_url = f"{trs_base_url}{trs_tool_id}/versions/{version}"
    else:
        # No version specified - fetch latest version from Dockstore
        try:
            # Query Dockstore API to get available versions (using encoded URL)
            versions_url = f"{trs_base_url}{encoded_tool_id}/versions"
            response = requests.get(versions_url, timeout=10)
            response.raise_for_status()
            versions = response.json()

            if versions and len(versions) > 0:
                # Get the first version (usually the latest/default)
                latest_version = versions[0].get("name") or versions[0].get("id")
                if latest_version:
                    trs_url = f"{trs_base_url}{trs_tool_id}/versions/{latest_version}"
                else:
                    # Fallback to just the tool ID without version
                    trs_url = f"{trs_base_url}{trs_tool_id}"
            else:
                # No versions found, use tool ID without version
                trs_url = f"{trs_base_url}{trs_tool_id}"
        except Exception:
            # If we can't fetch versions, just use the tool ID without version
            # Galaxy might handle this gracefully
            trs_url = f"{trs_base_url}{trs_tool_id}"

    return {"trs_url": trs_url}


def parse_trs_uri(trs_uri: str) -> Optional[Dict[str, str]]:
    """Parse a TRS URI into a full TRS URL.

    Args:
        trs_uri: TRS URI in format: trs://[#]workflow/github.com/org/repo/workflow_name[/version]
                 or trs://<full_dockstore_url>

    Returns:
        Dict with key 'trs_url' containing the full TRS API URL,
        or None if invalid
    """
    if not trs_uri.startswith(TRS_WORKFLOWS_PREFIX):
        return None

    # Remove trs:// prefix
    trs_content = trs_uri[len(TRS_WORKFLOWS_PREFIX) :]

    # Parse as a TRS ID (workflow/... or #workflow/...) to resolve versions
    return parse_trs_id(trs_content)


def import_workflow_from_trs(trs_uri: str, user_gi: "GalaxyInstance") -> Dict[str, Any]:
    """Import a workflow from a TRS endpoint using Galaxy's TRS import API.

    Args:
        trs_uri: TRS URI in format: trs://[#]workflow/github.com/org/repo/workflow_name[/version]
            Example: trs://workflow/github.com/iwc-workflows/parallel-accession-download/main

        user_gi: BioBlend GalaxyInstance for user API

    Returns:
        Workflow dict with 'id' and other metadata
    """
    trs_info = parse_trs_uri(trs_uri)
    if not trs_info:
        raise ValueError(f"Invalid TRS URI: {trs_uri}")

    # Create TRS import payload with full TRS URL
    # Example TRS URL: https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14
    trs_payload = {"archive_source": "trs_tool", "trs_url": trs_info["trs_url"]}

    # Use bioblend's _post method to import from TRS
    url = user_gi.workflows._make_url()
    workflow = user_gi.workflows._post(url=url, payload=trs_payload)

    return workflow


DOCKSTORE_TRS_BASE = "https://dockstore.org/api/ga4gh/trs/v2/tools/"


def _resolve_trs_url(trs_id: str) -> str:
    """Resolve a TRS identifier to a full TRS URL."""
    if trs_id.startswith(DOCKSTORE_TRS_BASE):
        return trs_id
    if trs_id.startswith(TRS_WORKFLOWS_PREFIX):
        trs_info = parse_trs_uri(trs_id)
        if not trs_info:
            raise ValueError(f"Invalid TRS URI: {trs_id}")
        return trs_info["trs_url"]
    # It's a short TRS ID
    trs_info = parse_trs_id(trs_id)
    if not trs_info:
        raise ValueError(f"Invalid TRS ID: {trs_id}")
    return trs_info["trs_url"]


def _encode_trs_url(trs_url: str) -> str:
    """URL-encode the tool ID portion of a TRS URL for API calls."""
    if not trs_url.startswith(DOCKSTORE_TRS_BASE):
        return trs_url
    tool_id_and_version = trs_url[len(DOCKSTORE_TRS_BASE) :]
    if "/versions/" in tool_id_and_version:
        tool_id, version = tool_id_and_version.split("/versions/", 1)
        return f"{DOCKSTORE_TRS_BASE}{quote(tool_id, safe='')}/versions/{version}"
    return f"{DOCKSTORE_TRS_BASE}{quote(tool_id_and_version, safe='')}"


def _find_primary_descriptor(files: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Find the primary descriptor file from a list of TRS files."""
    for f in files:
        if f.get("file_type") == "PRIMARY_DESCRIPTOR":
            return f
    for f in files:
        if f.get("path", "").endswith((".ga", ".gxwf.yml", ".gxwf.yaml")):
            return f
    return files[0] if files else None


def fetch_workflow_from_trs(trs_id: str) -> Dict[str, Any]:
    """Fetch a workflow definition directly from a TRS endpoint.

    Args:
        trs_id: TRS ID in format: [#]workflow/github.com/org/repo/workflow_name[/version]
                or a full TRS URL

    Returns:
        Workflow definition dict (gxformat2 or GA format)

    Raises:
        ValueError: If the TRS ID is invalid or workflow cannot be fetched
    """
    trs_url = _encode_trs_url(_resolve_trs_url(trs_id))
    files_url = f"{trs_url}/GALAXY/files"

    try:
        response = requests.get(files_url, timeout=30)
        response.raise_for_status()
        files = response.json()

        primary_file = _find_primary_descriptor(files)
        if not primary_file:
            raise ValueError(f"No workflow file found at TRS endpoint: {trs_url}")

        descriptor_url = f"{trs_url}/GALAXY/descriptor/{primary_file.get('path', '')}"
        file_response = requests.get(descriptor_url, timeout=30)
        file_response.raise_for_status()
        content = file_response.json().get("content", "")

        if not content:
            raise ValueError(f"Empty workflow content from TRS endpoint: {trs_url}")

        try:
            return json.loads(content)
        except json.JSONDecodeError:
            return yaml.safe_load(content)

    except requests.RequestException as e:
        raise ValueError(f"Failed to fetch workflow from TRS endpoint {trs_url}: {e}")


@lru_cache(maxsize=None)
def guess_tool_shed_url(tool_shed_fqdn: str) -> Optional[str]:
    if tool_shed_fqdn in MAIN_TOOLSHED_URL:
        return MAIN_TOOLSHED_URL
    else:
        # guess if tool shed is served over https or http
        https_tool_shed_url = f"https://{tool_shed_fqdn}"
        r = requests.get(https_tool_shed_url)
        if r.status_code == 200:
            return https_tool_shed_url
        else:
            http_tool_shed_url = f"http://{tool_shed_fqdn}"
            r = requests.get(http_tool_shed_url)
            if r.status_code == 200:
                return http_tool_shed_url
            else:
                warn(f"Could not connect to {tool_shed_fqdn}")
    return None


def get_toolshed_url_for_tool_id(tool_id: str) -> Optional[str]:
    components = tool_id.split("/repos")
    if len(components) > 1:
        tool_shed_fqdn = components[0]
        return guess_tool_shed_url(tool_shed_fqdn=tool_shed_fqdn)
    return None


def load_shed_repos(runnable):
    if runnable.type.name != "galaxy_workflow":
        return []
    path = runnable.path
    if path.endswith(".ga"):
        with tempfile.NamedTemporaryFile() as out:
            generate_tool_list_from_ga_workflow_files.generate_tool_list_from_workflow(
                [path], "Tools from workflows", out.name
            )
            with open(out.name) as f:
                tools = yaml.safe_load(f)["tools"]

    else:
        # It'd be better to just infer this from the tool shed ID somehow than
        # require explicit annotation like this... I think?
        with open(path) as f:
            workflow = yaml.safe_load(f)
        steps = workflow["steps"]
        if isinstance(steps, dict):
            steps = steps.values()
        tools = []
        for step in steps:
            repository = step.get("tool_shed_repository")
            if repository:
                repository["tool_panel_section_label"] = "Tools from workflows"
                tools.append(repository)
    for repo in tools:
        tool_shed = repo.get("tool_shed")
        if tool_shed:
            tool_shed_url = guess_tool_shed_url(tool_shed)
            if tool_shed_url:
                repo["tool_shed_url"] = tool_shed_url
    return tools


def _install_shed_repos_from_tools_info(
    tools_info: List[Dict[str, Any]],
    admin_gi: "GalaxyInstance",
    ignore_dependency_problems: bool,
    install_tool_dependencies: bool = False,
    install_resolver_dependencies: bool = True,
    install_repository_dependencies: bool = True,
    install_most_recent_revision: bool = False,
) -> Tuple[Optional[List[Any]], Optional[List[Any]]]:
    """Common logic for installing tool shed repositories from a tools_info list."""
    if not tools_info:
        return None, None

    install_tool_manager = shed_tools.InstallRepositoryManager(admin_gi)
    install_results = install_tool_manager.install_repositories(
        tools_info,
        default_install_tool_dependencies=install_tool_dependencies,
        default_install_resolver_dependencies=install_resolver_dependencies,
        default_install_repository_dependencies=install_repository_dependencies,
    )
    if install_most_recent_revision:  # for workflow autoupdates we also need the most recent tool versions
        update_results = install_tool_manager.update_repositories(
            tools_info,
            default_install_tool_dependencies=install_tool_dependencies,
            default_install_resolver_dependencies=install_resolver_dependencies,
            default_install_repository_dependencies=install_repository_dependencies,
        )
        install_results.errored_repositories.extend(update_results.errored_repositories)
        updated_repos = update_results.installed_repositories
    else:
        updated_repos = None

    if install_results.errored_repositories:
        if ignore_dependency_problems:
            warn(FAILED_REPOSITORIES_MESSAGE)
        else:
            raise Exception(FAILED_REPOSITORIES_MESSAGE)
    return install_results.installed_repositories, updated_repos


def install_shed_repos(
    runnable,
    admin_gi,
    ignore_dependency_problems,
    install_tool_dependencies=False,
    install_resolver_dependencies=True,
    install_repository_dependencies=True,
    install_most_recent_revision=False,
):
    tools_info = load_shed_repos(runnable)
    return _install_shed_repos_from_tools_info(
        tools_info,
        admin_gi,
        ignore_dependency_problems,
        install_tool_dependencies,
        install_resolver_dependencies,
        install_repository_dependencies,
        install_most_recent_revision,
    )


def install_shed_repos_for_workflow_id(
    workflow_id: str,
    user_gi: "GalaxyInstance",
    admin_gi: "GalaxyInstance",
    ignore_dependency_problems: bool,
    install_tool_dependencies: bool = False,
    install_resolver_dependencies: bool = True,
    install_repository_dependencies: bool = True,
    install_most_recent_revision: bool = False,
) -> Tuple[Optional[List[Any]], Optional[List[Any]]]:
    """Install tool shed repositories for a workflow that's already in Galaxy.

    This is used for TRS workflows that are imported via Galaxy's TRS API.
    We fetch the workflow definition from Galaxy and extract tool requirements.
    """
    # Fetch the workflow from Galaxy to get the GA format
    workflow_dict = user_gi.workflows.export_workflow_dict(workflow_id)

    # Use ephemeris to generate the tool list from the workflow
    with tempfile.NamedTemporaryFile(mode="w", suffix=".ga", delete=False) as wf_file:
        json.dump(workflow_dict, wf_file)
        wf_file.flush()
        wf_path = wf_file.name

    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as tool_file:
            tool_path = tool_file.name

        try:
            # Generate tool list from the GA workflow
            generate_tool_list_from_ga_workflow_files.generate_tool_list_from_workflow(
                [wf_path], "Tools from TRS workflow", tool_path
            )

            # Load the generated tool list
            with open(tool_path) as f:
                tools_data = yaml.safe_load(f)
                tools_info = tools_data.get("tools", []) if tools_data else []

            # Add tool shed URLs
            for repo in tools_info:
                tool_shed = repo.get("tool_shed")
                if tool_shed:
                    tool_shed_url = guess_tool_shed_url(tool_shed)
                    if tool_shed_url:
                        repo["tool_shed_url"] = tool_shed_url

            # Use common installation logic
            return _install_shed_repos_from_tools_info(
                tools_info,
                admin_gi,
                ignore_dependency_problems,
                install_tool_dependencies,
                install_resolver_dependencies,
                install_repository_dependencies,
                install_most_recent_revision,
            )
        finally:
            # Clean up tool list file
            if os.path.exists(tool_path):
                os.unlink(tool_path)
    finally:
        # Clean up workflow file
        if os.path.exists(wf_path):
            os.unlink(wf_path)


def import_workflow(path, admin_gi, user_gi, from_path=False):
    """Import a workflow path to specified Galaxy instance."""
    if not from_path:
        importer = BioBlendImporterGalaxyInterface(admin_gi=admin_gi, user_gi=user_gi)
        workflow = _raw_dict(path, importer)
        return user_gi.workflows.import_workflow_dict(workflow)
    else:
        path = os.path.abspath(path)
        workflow = user_gi.workflows.import_workflow_from_local_path(path)
        return workflow


def _raw_dict(path, importer=None):
    if path.endswith(".ga"):
        with open(path) as f:
            workflow = json.load(f)
    else:
        if importer is None:
            importer = DummyImporterGalaxyInterface()

        workflow_directory = os.path.dirname(path)
        workflow_directory = os.path.abspath(workflow_directory)
        with open(path) as f:
            workflow = yaml.safe_load(f)
            workflow = python_to_workflow(workflow, importer, workflow_directory)

    return workflow


def get_tool_ids_for_workflow(wf_dict: Dict[str, Any], tool_ids: Optional[List[str]] = None) -> List[str]:
    tool_ids = [] if tool_ids is None else tool_ids
    steps = wf_dict["steps"].values() if isinstance(wf_dict["steps"], dict) else wf_dict["steps"]
    for step in steps:
        if step.get("type", "tool") == "tool" and not step.get("run", {}).get("class") == "GalaxyWorkflow":
            tool_id = step["tool_id"]
            tool_ids.append(tool_id)
        elif step.get("type") == "subworkflow":  # GA SWF
            get_tool_ids_for_workflow(step["subworkflow"], tool_ids=tool_ids)
        elif step.get("run", {}).get("class") == "GalaxyWorkflow":  # gxformat2 SWF
            get_tool_ids_for_workflow(step["run"], tool_ids=tool_ids)
        else:
            continue
    return list(dict.fromkeys(tool_ids))


def find_tool_ids(path):
    workflow = _raw_dict(path)
    return get_tool_ids_for_workflow(workflow)


WorkflowOutput = namedtuple("WorkflowOutput", ["order_index", "output_name", "label", "optional"])


def remote_runnable_to_workflow_id(runnable):
    assert runnable.is_remote_workflow_uri
    parse_result = urlparse(runnable.uri)
    return parse_result.path[1:]


def describe_outputs(runnable, gi=None):
    """Return a list of :class:`WorkflowOutput` objects for target workflow."""
    if runnable.uri.startswith((GALAXY_WORKFLOWS_PREFIX, GALAXY_WORKFLOW_INSTANCE_PREFIX)):
        workflow_id = remote_runnable_to_workflow_id(runnable)
        assert gi is not None
        instance = runnable.uri.startswith(GALAXY_WORKFLOW_INSTANCE_PREFIX)
        workflow = get_dict_from_workflow(gi, workflow_id, instance)
    else:
        workflow = _raw_dict(runnable.path)

    outputs = []
    for order_index, step in workflow["steps"].items():
        optional = False
        if not step.get("tool_id"):
            # One of the parameter types ... need eliminate this guesswork on the Galaxy side
            tool_state = json.loads(step.get("tool_state", "{}"))
            optional = tool_state.get("optional", False)
        step_outputs = step.get("workflow_outputs", [])
        for step_output in step_outputs:
            output = WorkflowOutput(
                int(order_index),
                step_output["output_name"],
                step_output["label"],
                optional,
            )
            outputs.append(output)
    return outputs


class DummyImporterGalaxyInterface(ImporterGalaxyInterface):
    def import_workflow(self, workflow, **kwds):
        return None


def input_labels(workflow_path):
    """Get normalized labels for workflow artifact regardless of format."""
    steps = inputs_normalized(workflow_path=workflow_path)
    labels = []
    for step in steps:
        step_id = input_label(step)
        if step_id:
            labels.append(step_id)
    return labels


def required_input_steps(workflow_path):
    try:
        steps = inputs_normalized(workflow_path=workflow_path)
    except Exception:
        raise Exception("Input workflow could not be successfully normalized - try linting with planemo workflow_lint.")
    required_steps = []
    for input_step in steps:
        if input_step.get("optional", False) or input_step.get("default"):
            continue
        required_steps.append(input_step)
    return required_steps


def required_input_labels(workflow_path):
    return map(input_label, required_input_steps(workflow_path))


def input_label(input_step):
    """Get the normalized label of a step returned from inputs_normalized."""
    step_id = input_step.get("id") or input_step.get("label")
    return step_id


def output_labels(workflow_path):
    outputs = outputs_normalized(workflow_path=workflow_path)
    return [o["id"] for o in outputs]


def output_stubs_for_workflow(workflow_path, **kwds):
    """
    Return output labels and class.
    """
    if kwds.get("from_invocation"):
        return _job_outputs_template_from_invocation(workflow_path, kwds["galaxy_url"], kwds["galaxy_user_key"])
    outputs = {}
    for label in output_labels(workflow_path):
        if not label.startswith("_anonymous_"):
            outputs[label] = {"class": ""}
    return outputs


def job_template(workflow_path, **kwds):
    """Return a job template for specified workflow.

    A dictionary describing non-optional inputs that must be specified to
    run the workflow.
    """
    if kwds.get("from_invocation"):
        return _job_inputs_template_from_invocation(workflow_path, kwds["galaxy_url"], kwds["galaxy_user_key"])

    template = {}
    for required_input_step in required_input_steps(workflow_path):
        i_label = input_label(required_input_step)
        input_type = required_input_step["type"]
        if input_type == "data":
            template[i_label] = {
                "class": "File",
                "path": "todo_test_data_path.ext",
            }
        elif input_type == "collection":
            template[i_label] = {
                "class": "Collection",
                "collection_type": "list",
                "elements": [
                    {
                        "class": "File",
                        "identifier": "todo_element_name",
                        "path": "todo_test_data_path.ext",
                    }
                ],
            }
        elif input_type in ["string", "int", "float", "boolean", "color"]:
            template[i_label] = "todo_param_value"
        else:
            template[i_label] = {
                "TODO",  # Does this work yet?
            }
    return template


def _collection_elements_for_type(collection_type):
    """Generate appropriate sample elements for a collection type."""
    if collection_type == "paired":
        return [
            {
                "class": "File",
                "identifier": "forward",
                "path": "todo_test_data_path_forward.ext",
            },
            {
                "class": "File",
                "identifier": "reverse",
                "path": "todo_test_data_path_reverse.ext",
            },
        ]
    elif collection_type == "list:paired":
        return [
            {
                "class": "Collection",
                "type": "paired",
                "identifier": "todo_element_name",
                "elements": [
                    {
                        "class": "File",
                        "identifier": "forward",
                        "path": "todo_test_data_path_forward.ext",
                    },
                    {
                        "class": "File",
                        "identifier": "reverse",
                        "path": "todo_test_data_path_reverse.ext",
                    },
                ],
            }
        ]
    else:
        # Default to list
        return [
            {
                "class": "File",
                "identifier": "todo_element_name",
                "path": "todo_test_data_path.ext",
            }
        ]


def _build_template_and_metadata_from_inputs(
    all_inputs: List[Dict[str, Any]],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Build job template and metadata from normalized workflow inputs.

    Args:
        all_inputs: List of normalized input step definitions from gxformat2

    Returns:
        Tuple of (template_dict, metadata_dict)
    """
    template: Dict[str, Any] = {}
    metadata: Dict[str, Any] = {}
    for input_step in all_inputs:
        i_label = input_label(input_step)
        input_type = input_step["type"]
        input_doc = input_step.get("doc", "")
        is_optional = input_step.get("optional", False)
        default_value = input_step.get("default")
        has_default = default_value is not None
        input_format = input_step.get("format", "")
        collection_type = input_step.get("collection_type", "")

        # Store metadata for this input
        metadata[i_label] = {
            "type": input_type,
            "doc": input_doc,
            "optional": is_optional or has_default,
            "default": default_value,
            "format": input_format,
            "collection_type": collection_type,
        }

        if input_type == "data":
            template[i_label] = {
                "class": "File",
                "path": "todo_test_data_path.ext",
            }
        elif input_type == "collection":
            coll_type = collection_type or "list"
            template[i_label] = {
                "class": "Collection",
                "collection_type": coll_type,
                "elements": _collection_elements_for_type(coll_type),
            }
        elif input_type in ["string", "int", "float", "boolean", "color"]:
            # Use default value if available, otherwise use placeholder or false for booleans
            if has_default:
                template[i_label] = default_value
            elif input_type == "boolean":
                template[i_label] = False
            else:
                template[i_label] = "todo_param_value"
        else:
            template[i_label] = {
                "TODO",  # Does this work yet?
            }
    return template, metadata


def _job_template_with_metadata_from_dict(workflow_dict: Dict[str, Any]):
    """Generate job template and metadata from a workflow dictionary.

    This handles both GA format (.ga) and gxformat2 format workflows.

    Args:
        workflow_dict: Workflow definition dict

    Returns:
        Tuple of (template_dict, metadata_dict)
    """
    # Write workflow to temp file and use inputs_normalized
    # This handles both GA and gxformat2 formats consistently
    is_ga_format = (
        workflow_dict.get("a_galaxy_workflow", False)
        or "steps" in workflow_dict
        and isinstance(workflow_dict.get("steps"), dict)
    )

    suffix = ".ga" if is_ga_format else ".gxwf.yml"

    with tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False) as wf_file:
        if is_ga_format:
            json.dump(workflow_dict, wf_file)
        else:
            yaml.dump(workflow_dict, wf_file)
        wf_file.flush()
        wf_path = wf_file.name

    try:
        all_inputs = inputs_normalized(workflow_path=wf_path)
    except Exception:
        raise Exception("Input workflow could not be successfully normalized from TRS endpoint.")
    finally:
        if os.path.exists(wf_path):
            os.unlink(wf_path)

    return _build_template_and_metadata_from_inputs(all_inputs)


def is_trs_identifier(identifier: str) -> bool:
    """Check if the identifier is a TRS ID or TRS URI.

    Args:
        identifier: The workflow identifier to check

    Returns:
        True if it's a TRS ID or URI, False otherwise
    """
    # Full Dockstore TRS URL
    if identifier.startswith(DOCKSTORE_TRS_BASE):
        return True
    # TRS URI
    if identifier.startswith(TRS_WORKFLOWS_PREFIX):
        return True
    # Short TRS ID format
    if identifier.startswith(("workflow/", "tool/", "#workflow/", "#tool/")) and "/github.com/" in identifier:
        return True
    return False


def job_template_with_metadata(workflow_path, **kwds):
    """Return a job template with metadata for each input.

    Returns a tuple of (template_dict, metadata_dict) where metadata_dict
    contains type, doc (description), optional status, default value, and format for each input label.

    The workflow_path can be:
    - A local file path to a workflow file
    - A TRS ID (e.g., workflow/github.com/org/repo/workflow_name)
    - A TRS URI (e.g., trs://workflow/github.com/org/repo/workflow_name)
    - A full Dockstore TRS URL
    """
    if kwds.get("from_invocation"):
        # For invocation-based templates, we don't have metadata
        return _job_inputs_template_from_invocation(workflow_path, kwds["galaxy_url"], kwds["galaxy_user_key"]), {}

    # Check if this is a TRS identifier
    if is_trs_identifier(workflow_path):
        # Fetch workflow from TRS and write to temp file for processing
        workflow_dict = fetch_workflow_from_trs(workflow_path)
        return _job_template_with_metadata_from_dict(workflow_dict)

    try:
        all_inputs = inputs_normalized(workflow_path=workflow_path)
    except Exception:
        raise Exception("Input workflow could not be successfully normalized - try linting with planemo workflow_lint.")

    return _build_template_and_metadata_from_inputs(all_inputs)


def new_workflow_associated_path(workflow_path, suffix="tests"):
    """Generate path for test or job YAML file next to workflow."""
    base, input_ext = os.path.splitext(workflow_path)
    # prefer -tests.yml but if the author uses underscores or .yaml respect that.
    sep = "-"
    if "_" in base and "-" not in base:
        sep = "_"
    ext = "yml"
    if "yaml" in input_ext:
        ext = "yaml"
    return base + sep + suffix + "." + ext


def rewrite_job_file(input_file, output_file, job):
    """Rewrite a job file with galaxy_ids for upload_data subcommand"""
    with open(input_file) as f:
        job_contents = yaml.safe_load(f)
        for job_input, job_input_name in job_contents.items():
            if isinstance(job[job_input], dict):  # dataset or collection
                job_contents[job_input] = {"class": job_input_name["class"], "galaxy_id": job[job_input]["id"]}
            # else: presumably a parameter, no need to modify
    with open(output_file, "w") as f:
        yaml.dump(job_contents, f)


def get_workflow_from_invocation_id(invocation_id, galaxy_url, galaxy_api_key):
    user_gi = gi(url=galaxy_url, key=galaxy_api_key)
    workflow_id = user_gi.invocations.show_invocation(invocation_id)["workflow_id"]
    workflow = get_dict_from_workflow(user_gi, workflow_id, instance=True)
    workflow_name = "-".join(workflow["name"].split())
    with open(f"{workflow_name}.ga", "w") as workflow_out:
        json.dump(workflow, workflow_out, ensure_ascii=False, indent=4)
    return workflow_name


def _elements_to_test_def(
    elements: List[Dict[str, Any]],
    test_data_base_path: str,
    download_function: Callable,
    definition_style: str = "input",
):
    element_test_def = []
    output_element_test_def = {}
    if definition_style == "output":
        elements = elements[:1]
    for element in elements:
        if element["element_type"] == "dataset_collection":
            nested_elements = _elements_to_test_def(
                element["object"]["elements"],
                test_data_base_path,
                download_function,
                definition_style=definition_style,
            )
            test_def = {}
            if definition_style == "input":
                test_def["class"] = "Collection"
                test_def["type"] = element["object"]["collection_type"]
                test_def["identifier"] = element["element_identifier"]
                test_def["elements"] = nested_elements
                element_test_def.append(test_def)
            else:
                output_element_test_def[element["element_identifier"]] = {"elements": nested_elements}
        elif element["element_type"] == "hda":
            ext = element["object"]["file_ext"]
            path = f"{test_data_base_path}_{element['element_identifier']}.{ext}"
            download_function(
                element["object"]["id"],
                use_default_filename=False,
                file_path=path,
            )
            if definition_style == "input":
                test_def = {
                    "class": "File",
                    "identifier": element["element_identifier"],
                    "path": path,
                }
                element_test_def.append(test_def)
            else:
                output_element_test_def[element["element_identifier"]] = {"path": path}
    if definition_style == "input":
        return element_test_def
    else:
        return output_element_test_def


def _job_inputs_template_from_invocation(invocation_id, galaxy_url, galaxy_api_key):
    user_gi = gi(url=galaxy_url, key=galaxy_api_key)
    invocation = user_gi.invocations.show_invocation(invocation_id)
    template = {}
    for input_step in invocation["inputs"].values():
        if input_step["src"] == "hda":
            ext = user_gi.datasets.show_dataset(input_step["id"])["extension"]
            user_gi.datasets.download_dataset(
                input_step["id"], use_default_filename=False, file_path=f"test-data/{input_step['label']}.{ext}"
            )
            template[input_step["label"]] = {
                "class": "File",
                "path": f"test-data/{input_step['label']}.{ext}",
                "filetype": ext,
            }
        elif input_step["src"] == "hdca":
            collection = user_gi.dataset_collections.show_dataset_collection(input_step["id"])
            test_def = {
                "class": "Collection",
                "collection_type": collection["collection_type"],
                "elements": _elements_to_test_def(
                    collection["elements"],
                    test_data_base_path=f"test-data/{input_step['label']}",
                    download_function=user_gi.datasets.download_dataset,
                ),
            }
            template[input_step["label"]] = test_def
    for param, param_step in invocation["input_step_parameters"].items():
        template[param] = param_step["parameter_value"]

    return template


def _job_outputs_template_from_invocation(invocation_id, galaxy_url, galaxy_api_key):
    user_gi = gi(url=galaxy_url, key=galaxy_api_key)
    invocation = user_gi.invocations.show_invocation(invocation_id)
    outputs = {}
    for label, output in invocation["outputs"].items():
        ext = user_gi.datasets.show_dataset(output["id"])["extension"]
        user_gi.datasets.download_dataset(
            output["id"], use_default_filename=False, file_path=f"test-data/{label}.{ext}"
        )
        outputs[label] = {"path": f"test-data/{label}.{ext}"}
    for label, output in invocation["output_collections"].items():
        collection = user_gi.dataset_collections.show_dataset_collection(output["id"])
        element_tests = _elements_to_test_def(
            collection["elements"],
            test_data_base_path=f"test-data/{label}",
            download_function=user_gi.datasets.download_dataset,
            definition_style="outputs",
        )
        outputs[label] = {"element_tests": element_tests}
    return outputs


__all__ = (
    "describe_outputs",
    "DOCKSTORE_TRS_BASE",
    "fetch_workflow_from_trs",
    "import_workflow",
    "import_workflow_from_trs",
    "is_trs_identifier",
    "TRS_WORKFLOWS_PREFIX",
)
