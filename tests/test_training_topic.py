"""Training:topic functions."""
import os
import shutil

from planemo.training.topic import Topic
from planemo.training.utils import load_yaml
from .test_utils import TEST_DATA_DIR


def test_topic_init():
    """Test :func:`planemo.training.topic.Topic.init`."""
    # test requirement with default parameter
    topic = Topic()
    assert topic.name == "new_topic"
    assert topic.type == "use"
    assert topic.title == "The new topic"
    assert topic.summary == "Summary"
    assert topic.docker_image == ""
    assert "maintainers" in topic.maintainers
    assert topic.parent_dir == "topics"
    assert topic.dir == "topics/new_topic"
    assert topic.requirements[0].topic_name == "introduction"
    # test requirement with non default
    topic = Topic(name="topic2", target="admin", title="The 2nd topic", summary="", parent_dir="dir")
    assert topic.name == "topic2"
    assert topic.type == "admin"
    assert topic.title == "The 2nd topic"
    assert topic.summary == ""
    assert topic.parent_dir == "dir"
    assert topic.dir == "dir/topic2"
    assert len(topic.requirements) == 0


def test_topic_init_from_kwds():
    """Test :func:`planemo.training.topic.Topic.init_from_kwds`."""
    topic = Topic()
    topic.init_from_kwds(
        {"topic_name": "topic", "topic_title": "New topic", "topic_target": "admin", "topic_summary": "Topic summary"}
    )
    assert topic.name == "topic"
    assert topic.type == "admin"
    assert topic.title == "New topic"
    assert topic.summary == "Topic summary"
    assert topic.dir == "topics/topic"
    assert len(topic.requirements) == 0


def test_topic_init_from_metadata():
    """Test :func:`planemo.training.topic.Topic.init_from_metadata`."""
    topic = Topic()
    os.makedirs(topic.dir)
    shutil.copy(os.path.join(TEST_DATA_DIR, "training_metadata.yaml"), topic.metadata_fp)
    topic.init_from_metadata()
    assert topic.name == "test"
    assert topic.title == "Test"
    assert topic.summary == "Summary"
    assert topic.requirements[0].topic_name == "introduction"
    assert topic.requirements[0].tutorials == ["peaks2genes"]
    assert "maintainer1" in topic.maintainers
    shutil.rmtree(topic.parent_dir)


def test_topic_get_requirements():
    """Test :func:`planemo.training.topic.Topic.get_requirements`."""
    topic = Topic()
    reqs = topic.get_requirements()
    assert len(reqs) == 1
    assert "topic_name" in reqs[0]


def test_topic_export_metadata_to_ordered_dict():
    """Test :func:`planemo.training.topic.Topic.export_metadata_to_ordered_dict`."""
    topic = Topic()
    metadata = topic.export_metadata_to_ordered_dict()
    assert "name" in metadata
    assert metadata["name"] == "new_topic"
    assert "type" in metadata
    assert "title" in metadata
    assert "summary" in metadata
    assert "requirements" in metadata
    assert "docker_image" in metadata
    assert "maintainers" in metadata


def test_topic_set_paths():
    """Test :func:`planemo.training.topic.Topic.set_paths`."""
    new_name = "the_new_name"
    topic = Topic()
    topic.name = new_name
    topic.set_paths()
    assert new_name in topic.dir
    assert new_name in topic.img_folder
    assert new_name in topic.tuto_folder
    assert new_name in topic.index_fp
    assert new_name in topic.readme_fp
    assert new_name in topic.metadata_fp
    assert new_name in topic.docker_folder
    assert new_name in topic.dockerfile_fp


def test_topic_exists():
    """Test :func:`planemo.training.topic.Topic.exists`."""
    topic = Topic()
    os.makedirs(topic.dir)
    assert topic.exists()
    shutil.rmtree(topic.parent_dir)


def test_topic_create_topic_structure():
    """Test :func:`planemo.training.topic.Topic.create_topic_structure`."""
    topic = Topic()
    topic.create_topic_structure()
    topic_name = "new_topic"
    topic_title = "The new topic"
    # check the folder and its structure
    assert topic.exists()
    assert os.path.exists(topic.img_folder)
    assert os.path.exists(topic.tuto_folder)
    # create the index.md and the topic name
    assert os.path.exists(topic.index_fp)
    with open(topic.index_fp) as fh:
        assert topic_name in fh.read()
    # create the README.md and the topic name
    assert os.path.exists(topic.readme_fp)
    with open(topic.readme_fp) as fh:
        assert topic_title in fh.read()
    # check metadata content
    assert os.path.exists(topic.metadata_fp)
    metadata = load_yaml(topic.metadata_fp)
    assert metadata["name"] == topic_name
    # check dockerfile
    assert os.path.exists(topic.dockerfile_fp)
    with open(topic.dockerfile_fp) as fh:
        assert topic_name in fh.read()
    with open(topic.dockerfile_fp) as fh:
        assert topic_title in fh.read()
    # check in metadata directory
    assert os.path.exists(os.path.join("metadata", "%s.yaml" % topic_name))
    # clean
    shutil.rmtree(topic.parent_dir)
    shutil.rmtree("metadata")
