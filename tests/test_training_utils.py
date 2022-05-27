"""Training:utils functions."""
import os

from planemo.training.utils import (
    load_yaml,
    Requirement,
    save_to_yaml,
)
from .test_utils import TEST_DATA_DIR

metadata_fp = os.path.join(TEST_DATA_DIR, "training_metadata.yaml")


def test_load_yaml():
    """Test :func:`planemo.training.utils.load_yaml`."""
    metadata = load_yaml(metadata_fp)
    # test if name there
    assert metadata["name"] == "test"
    # test if order of material is conserved
    assert metadata["maintainers"][0] == "maintainer1"


def test_save_to_yaml():
    """Test :func:`planemo.training.utils.save_to_yaml`."""
    metadata = load_yaml(metadata_fp)
    new_metadata_fp = "metadata.yaml"
    save_to_yaml(metadata, new_metadata_fp)
    assert os.path.exists(new_metadata_fp)
    os.remove(new_metadata_fp)


def test_requirement_init():
    """Test :func:`planemo.training.utils.Requirement.init`."""
    # test requirement with default parameter
    req = Requirement()
    assert req.title is None
    assert req.type == "internal"
    assert req.topic_name == "introduction"
    # test requirement with non default
    req = Requirement(title="Introduction", req_type="external", link="URL")
    assert req.title == "Introduction"
    assert req.type == "external"
    assert req.link == "URL"


def test_requirement_init_from_dict():
    """Test :func:`planemo.training.utils.Requirement.init_from_dict`."""
    req = Requirement()
    req.init_from_dict({"title": "The Requirement", "type": "external", "link": "http://URL"})
    assert req.title == "The Requirement"
    assert req.type == "external"
    assert req.link == "http://URL"


def test_requirement_export_to_ordered_dict():
    """Test :func:`planemo.training.utils.Requirement.export_to_ordered_dict`."""
    req = Requirement()
    exp_req = req.export_to_ordered_dict()
    assert "type" in exp_req
    assert exp_req["type"] == "internal"
    assert "topic_name" in exp_req
    assert "link" not in exp_req
