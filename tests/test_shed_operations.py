"""Test some lower-level utilities in planemo.shed."""

import os

from planemo import shed
from .test_utils import (
    mock_shed_context,
    TEST_REPOS_DIR,
)


def test_find_repository_id():
    with mock_shed_context() as shed_context:
        repo_id = shed.find_repository_id(
            ctx=None,
            shed_context=shed_context,
            path=".",
            name="test_repo_1",
            owner="iuc",
        )
        assert repo_id == "r1"


def test_find_repository_id_missing():
    with mock_shed_context() as shed_context:
        repo_id = shed.find_repository_id(
            ctx=None, shed_context=shed_context, path=".", name="test_repo_absent", owner="iuc", allow_none=True
        )
        assert repo_id is None


def test_find_repository_id_missing_exception():
    with mock_shed_context() as shed_context:
        exception = None
        try:
            shed.find_repository_id(ctx=None, shed_context=shed_context, path=".", name="test_repo_absent", owner="iuc")
        except Exception as e:
            exception = e
        assert exception is not None


def test_find_category_ids():
    with mock_shed_context() as shed_context:
        category_ids = shed.find_category_ids(shed_context.tsi, ["Text Manipulation"])
        assert category_ids == ["c1"]


def test_create_simple():
    with mock_shed_context() as shed_context:
        path = os.path.join(TEST_REPOS_DIR, "single_tool")
        repo_config = shed.shed_repo_config(shed_context, path)
        create_response = shed.create_repository_for(
            None,
            shed_context.tsi,
            "single_tool",
            repo_config,
        )
        assert "id" in create_response
        repo_id = shed.find_repository_id(
            ctx=None, shed_context=shed_context, path=".", name="single_tool", owner="iuc"
        )
        assert repo_id == create_response["id"]


def test_suite_repositories_different_owners():
    with mock_shed_context() as shed_context:
        path = os.path.join(TEST_REPOS_DIR, "multi_repos_flat_configured_owners")
        repo_config = shed.shed_repo_config(shed_context, path)
        assert (
            '<repository owner="iuc" name="cs-cat1" />'
            in repo_config["repositories"]["suite_cat"]["_files"]["repository_dependencies.xml"]
        ), repo_config
        assert (
            '<repository owner="devteam" name="cs-cat2" />'
            in repo_config["repositories"]["suite_cat"]["_files"]["repository_dependencies.xml"]
        ), repo_config
