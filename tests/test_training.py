"""Training training functions."""
import json
import os
import shutil

import pytest

from planemo import cli
from planemo.runnable import for_path
from planemo.training import Training
from .test_utils import (
    skip_if_environ,
    TEST_DATA_DIR,
)

datatype_fp = os.path.join(TEST_DATA_DIR, "training_datatypes.yaml")
tuto_fp = os.path.join(TEST_DATA_DIR, "training_tutorial.md")
tuto_wo_zenodo_fp = os.path.join(TEST_DATA_DIR, "training_tutorial_wo_zenodo.md")
zenodo_link = "https://zenodo.org/record/1321885"
# load a workflow generated from Galaxy
WF_FP = os.path.join(TEST_DATA_DIR, "training_workflow.ga")
with open(WF_FP) as wf_f:
    wf = json.load(wf_f)
# load wf_param_values (output of tutorial.get_wf_param_values on wf['steps']['4'])
with open(os.path.join(TEST_DATA_DIR, "training_wf_param_values.json")) as wf_param_values_f:
    wf_param_values = json.load(wf_param_values_f)
# configuration
RUNNABLE = for_path(WF_FP)
CTX = cli.PlanemoCliContext()
CTX.planemo_directory = "/tmp/planemo-test-workspace"
KWDS = {
    "topic_name": "my_new_topic",
    "topic_title": "New topic",
    "topic_target": "use",
    "topic_summary": "Topic summary",
    "tutorial_name": "new_tuto",
    "tutorial_title": "Title of tuto",
    "hands_on": True,
    "slides": True,
    "workflow": None,
    "workflow_id": None,
    "zenodo_link": None,
    "datatypes": os.path.join(TEST_DATA_DIR, "training_datatypes.yaml"),
    "templates": None,
    # planemo configuation
    "conda_auto_init": True,
    "conda_auto_install": True,
    "conda_copy_dependencies": False,
    "conda_debug": False,
    "conda_dependency_resolution": False,
    "conda_ensure_channels": "iuc,conda-forge,bioconda,defaults",
    "conda_exec": None,
    "conda_prefix": None,
    "conda_use_local": False,
    "brew_dependency_resolution": False,
    "daemon": False,
    "database_connection": None,
    "database_type": "auto",
    "dependency_resolvers_config_file": None,
    "docker": False,
    "docker_cmd": "docker",
    "docker_extra_volume": None,
    "docker_galaxy_image": "quay.io/bgruening/galaxy",
    "docker_host": None,
    "docker_sudo": False,
    "docker_sudo_cmd": "sudo",
    "engine": "galaxy",
    "extra_tools": (),
    "file_path": None,
    "galaxy_api_key": None,
    "galaxy_branch": None,
    "galaxy_email": "planemo@galaxyproject.org",
    "galaxy_root": None,
    "galaxy_single_user": True,
    "galaxy_source": None,
    "galaxy_url": None,
    "host": "127.0.0.1",
    "ignore_dependency_problems": False,
    "install_galaxy": False,
    "job_config_file": None,
    "mulled_containers": False,
    "no_cleanup": False,
    "no_cache_galaxy": False,
    "no_dependency_resolution": True,
    "non_strict_cwl": False,
    "pid_file": None,
    "port": "9090",
    "postgres_database_host": None,
    "postgres_database_port": None,
    "postgres_database_user": "postgres",
    "postgres_psql_path": "psql",
    "profile": None,
    "shed_dependency_resolution": False,
    "shed_install": True,
    "shed_tool_conf": None,
    "shed_tool_path": None,
    "skip_client_build": True,
    "skip_venv": False,
    "test_data": None,
    "tool_data_table": None,
    "tool_dependency_dir": None,
}


def test_training_init() -> None:
    """Test :func:`planemo.training.Training.init`."""
    train = Training(KWDS)
    assert train.topics_dir == "topics"
    assert train.topic is not None
    assert train.tuto is None


def test_training_init_training() -> None:
    """Test :func:`planemo.training.Training.init_training`."""
    train = Training(KWDS)
    # new topic, nothing else
    train.kwds["tutorial_name"] = None
    train.kwds["slides"] = None
    train.kwds["workflow"] = None
    train.kwds["workflow_id"] = None
    train.kwds["zenodo_link"] = None
    train.init_training(CTX)
    assert os.path.exists(train.topic.dir)
    assert not os.listdir(os.path.join(train.topic.dir, "tutorials"))
    # no new topic, no tutorial name but hands-on
    train.kwds["slides"] = True
    with pytest.raises(Exception, match="A tutorial name is needed to create the skeleton of a tutorial slide deck"):
        train.init_training(CTX)
    # no new topic, no tutorial name but workflow
    train.kwds["workflow"] = WF_FP
    train.kwds["slides"] = False
    with pytest.raises(
        Exception, match="A tutorial name is needed to create the skeleton of the tutorial from a workflow"
    ):
        train.init_training(CTX)
    # no new topic, no tutorial name but zenodo
    train.kwds["workflow"] = None
    train.kwds["zenodo_link"] = zenodo_link
    with pytest.raises(Exception, match="A tutorial name is needed to add Zenodo information"):
        train.init_training(CTX)
    # no new topic, new tutorial
    train.kwds["tutorial_name"] = "new_tuto"
    train.kwds["workflow"] = None
    train.kwds["zenodo_link"] = None
    train.init_training(CTX)
    assert os.path.exists(train.tuto.dir)
    assert os.path.exists(train.tuto.tuto_fp)
    with open(train.tuto.tuto_fp) as fh:
        assert train.kwds["tutorial_title"] in fh.read()
    # clean after
    shutil.rmtree(train.topics_dir)
    shutil.rmtree("metadata")


def create_existing_tutorial(exit_tuto_name, tuto_fp, topic) -> None:
    exist_tuto_dir = os.path.join(topic.dir, "tutorials", exit_tuto_name)
    os.makedirs(exist_tuto_dir)
    shutil.copyfile(tuto_fp, os.path.join(exist_tuto_dir, "tutorial.md"))


def test_training_check_topic_init_tuto() -> None:
    """Test :func:`planemo.training.Training.check_topic_init_tuto`."""
    train = Training(KWDS)
    # no topic
    with pytest.raises(Exception, match="The topic my_new_topic does not exists. It should be created"):
        train.check_topic_init_tuto()
    # add topic
    train.kwds["tutorial_name"] = None
    train.kwds["slides"] = None
    train.kwds["workflow"] = None
    train.kwds["workflow_id"] = None
    train.kwds["zenodo_link"] = None
    train.init_training(CTX)
    train.kwds["tutorial_name"] = "existing_tutorial"
    create_existing_tutorial("existing_tutorial", tuto_fp, train.topic)
    train.check_topic_init_tuto()
    assert train.tuto.name == train.kwds["tutorial_name"]
    assert train.tuto.datatype_fp
    # clean after
    shutil.rmtree(train.topics_dir)
    shutil.rmtree("metadata")


def test_fill_data_library() -> None:
    """Test :func:`planemo.training.fill_data_library`."""
    train = Training(KWDS)
    train.kwds["tutorial_name"] = None
    train.kwds["slides"] = False
    train.kwds["hands_on"] = False
    train.init_training(CTX)
    train.kwds["tutorial_name"] = "existing_tutorial"
    create_existing_tutorial("existing_tutorial", tuto_wo_zenodo_fp, train.topic)
    # no Zenodo link
    train.kwds["zenodo_link"] = None
    with pytest.raises(
        Exception, match="A Zenodo link should be provided either in the metadata file or as argument of the command"
    ):
        train.fill_data_library(CTX)
    # with a given Zenodo link and no Zenodo in metadata
    train.kwds["zenodo_link"] = zenodo_link
    train.fill_data_library(CTX)
    with open(train.tuto.data_lib_fp) as fh:
        assert "DOI: 10.5281/zenodo.1321885" in fh.read()
    with open(train.tuto.tuto_fp) as fh:
        assert "zenodo_link: %s" % zenodo_link in fh.read()
    # with a given Zenodo link and Zenodo in metadata
    new_z_link = "https://zenodo.org/record/1324204"
    train.kwds["zenodo_link"] = new_z_link
    train.tuto = None
    train.fill_data_library(CTX)
    assert train.tuto
    with open(train.tuto.data_lib_fp) as fh:
        assert "DOI: 10.5281/zenodo.1324204" in fh.read()
    with open(train.tuto.tuto_fp) as fh:
        assert "zenodo_link: %s" % new_z_link in fh.read()
    # with no given Zenodo link
    train.kwds["zenodo_link"] = None
    train.fill_data_library(CTX)
    with open(train.tuto.data_lib_fp) as fh:
        assert "DOI: 10.5281/zenodo.1324204" in fh.read()
    with open(train.tuto.tuto_fp) as fh:
        assert "zenodo_link: %s" % new_z_link in fh.read()
    # clean after
    shutil.rmtree(train.topics_dir)
    shutil.rmtree("metadata")


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_generate_tuto_from_wf() -> None:
    """Test :func:`planemo.training.generate_tuto_from_wf`."""
    train = Training(KWDS)
    train.kwds["tutorial_name"] = None
    train.kwds["slides"] = False
    train.init_training(CTX)
    train.kwds["tutorial_name"] = "existing_tutorial"
    create_existing_tutorial("existing_tutorial", tuto_fp, train.topic)
    # no workflow
    train.kwds["workflow"] = None
    with pytest.raises(
        Exception,
        match="A path to a local workflow or the id of a workflow on a running Galaxy instance should be provided",
    ):
        train.generate_tuto_from_wf(CTX)
    # with workflow
    train.kwds["workflow"] = WF_FP
    train.generate_tuto_from_wf(CTX)
    assert_file_contains(
        train.tuto.tuto_fp,
        "{% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.71) %} with the following parameters:",
    )
    assert os.path.exists(train.tuto.wf_fp)
    # clean after
    shutil.rmtree(train.topics_dir)
    shutil.rmtree("metadata")


def assert_file_contains(file: str, text: str) -> None:
    with open(file) as fh:
        contents = fh.read()
    if text not in contents:
        template = "Expected file [%s] to contain [%s], it did not - contents [%s]"
        message = template % (file, text, contents)
        raise AssertionError(message)
