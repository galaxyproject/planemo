"""Unit tests for the generated job and TPV configurations (``planemo.galaxy.job_conf``)."""

import os

import yaml

from planemo.galaxy.config import _handle_job_config_file
from planemo.galaxy.job_conf import (
    TPV_CONFIG_FILE_NAME,
    TPV_DEFAULT_TOOL,
    TPV_DESTINATION_NAME,
    tpv_install_command,
)

SHARED_TPV_DATABASE = "https://gxy.io/tpv/db.yml"


def _write(tmp_path, **kwds):
    """Write a job configuration for the given Planemo options and load the result."""
    kwds.setdefault("docker", True)
    _handle_job_config_file(str(tmp_path), "planemo", None, set(), kwds)
    with open(kwds["job_config_file"]) as fh:
        job_config = yaml.safe_load(fh)
    return job_config, _load_tpv_config(tmp_path)


def _user_config(directory, contents="tools: {}\n"):
    """Write a TPV configuration of the user and return its path."""
    os.makedirs(str(directory), exist_ok=True)
    path = os.path.join(str(directory), "user.tpv.yml")
    with open(path, "w") as fh:
        fh.write(contents)
    return path


def _load_tpv_config(tmp_path):
    tpv_config_path = os.path.join(str(tmp_path), TPV_CONFIG_FILE_NAME)
    if not os.path.exists(tpv_config_path):
        return None
    with open(tpv_config_path) as fh:
        return yaml.safe_load(fh)


def _tpv_environment(job_config):
    return job_config["execution"]["environments"]["tpv"]


def _tpv_references(job_config):
    tpv_environment = _tpv_environment(job_config)
    for key in ("tpv_configs", "tpv_config_files"):
        if key in tpv_environment:
            return key, tpv_environment[key]
    raise AssertionError(f"No TPV configuration referenced by {tpv_environment}")


def _tpv_params(tpv_config):
    return tpv_config["destinations"][TPV_DESTINATION_NAME].get("params") or {}


def test_no_tpv_config_by_default(tmp_path):
    """Without --tpv (or --tpv_config_file) no TPV configuration is generated."""
    job_config, tpv_config = _write(tmp_path)
    assert tpv_config is None
    assert job_config["execution"]["default"] == "local"


def test_tpv_singularity_destination(tmp_path):
    """Singularity is passed through to the generated TPV destination."""
    _, tpv_config = _write(tmp_path, docker=False, singularity=True, tpv_config_file=[_user_config(tmp_path)])
    params = _tpv_params(tpv_config)
    assert params["singularity_enabled"] is True
    assert "docker_enabled" not in params


def test_tpv_maps_all_tools_to_generated_destination(tmp_path):
    """--tpv_config_file replaces the tool based environment assignment by TPV."""
    job_config, tpv_config = _write(tmp_path, tpv_config_file=[_user_config(tmp_path)])

    assert job_config["execution"]["default"] == "tpv"
    assert "tools" not in job_config

    config_key, references = _tpv_references(job_config)
    assert config_key in ("tpv_configs", "tpv_config_files")
    assert references[0].endswith(TPV_CONFIG_FILE_NAME)

    destination = tpv_config["destinations"][TPV_DESTINATION_NAME]
    assert destination["runner"] == "local"
    assert destination["params"]["docker_enabled"] is True
    assert "docker" in destination["scheduling"]["accept"]


def test_tpv_without_config_file_is_unchanged(tmp_path):
    """--tpv alone keeps the generated configuration, including the shared database."""
    job_config, tpv_config = _write(tmp_path, tpv=True)

    assert tpv_config is None
    assert job_config["execution"]["default"] == "tpv"
    _, references = _tpv_references(job_config)
    assert SHARED_TPV_DATABASE in references
    # the tool based environment assignment of the generated configuration is kept
    assert job_config["tools"]


def test_tpv_config_is_self_contained(tmp_path):
    """A TPV configuration needs no shared database and keeps unlisted tools."""
    job_config, tpv_config = _write(tmp_path, tpv_config_file=[_user_config(tmp_path)])

    _, references = _tpv_references(job_config)
    assert SHARED_TPV_DATABASE not in references
    assert all("://" not in reference for reference in references)

    # Tools that no configuration mentions inherit the generated default destination.
    assert tpv_config["global"]["default_inherits"] == TPV_DEFAULT_TOOL
    assert TPV_DEFAULT_TOOL in tpv_config["tools"]


def test_tpv_keeps_container_settings(tmp_path):
    """The generated destination reproduces the container settings of the job config."""
    job_config, tpv_config = _write(
        tmp_path, tpv_config_file=[_user_config(tmp_path)], docker=True, docker_extra_volume=["/data"]
    )

    environment = job_config["execution"]["environments"]["local"]
    params = _tpv_params(tpv_config)
    assert params["docker_enabled"] == environment["docker_enabled"]
    assert params["docker_volumes"] == environment["docker_volumes"]


def test_tpv_config_file_implies_tpv(tmp_path):
    """--tpv_config_file enables TPV and is layered after the generated configuration."""
    user_config = _user_config(
        tmp_path, "tools:\n  tool1:\n    params:\n      docker_run_extra_arguments: --net none\n"
    )
    job_config, _ = _write(tmp_path, tpv_config_file=[user_config])

    assert job_config["execution"]["default"] == "tpv"
    _, references = _tpv_references(job_config)
    assert references[-1] == user_config
    assert references[0].endswith(TPV_CONFIG_FILE_NAME)
    assert SHARED_TPV_DATABASE not in references


def test_tpv_config_file_keeps_requested_order(tmp_path):
    """Several configurations are layered in the order they were passed in."""
    first = _user_config(tmp_path / "first")
    second = _user_config(tmp_path / "second")

    job_config, _ = _write(tmp_path, tpv_config_file=[first, second])
    _, references = _tpv_references(job_config)
    assert references[1:] == [first, second]


def test_tpv_galaxy_24_2_uses_config_files(tmp_path):
    """Galaxy 24.2 references TPV configurations with tpv_config_files."""
    job_config, _ = _write(tmp_path, tpv_config_file=[_user_config(tmp_path)], galaxy_version="24.2")

    assert _tpv_environment(job_config)["runner"] == "dynamic"
    config_key, references = _tpv_references(job_config)
    assert config_key == "tpv_config_files"
    assert references[0].endswith(TPV_CONFIG_FILE_NAME)


def test_tpv_destination_uses_requested_runner(tmp_path):
    """The generated destination runs on the configured job runner."""
    _, tpv_config = _write(tmp_path, runner="slurm", tpv_config_file=[_user_config(tmp_path)])
    destination = tpv_config["destinations"][TPV_DESTINATION_NAME]
    assert destination["runner"] == "slurm"
    assert destination["params"]["native_specification"] == ""


def test_tpv_install_command():
    """TPV is installed in the Galaxy virtualenv only when TPV is requested."""
    assert tpv_install_command({}) is None
    assert "total-perspective-vortex" in tpv_install_command({"tpv": True})
    assert "total-perspective-vortex" in tpv_install_command({"tpv_config_file": ["user.tpv.yml"]})
