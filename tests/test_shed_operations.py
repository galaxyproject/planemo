""" Test some lower-level utilities in planemo.shed.
"""
import os

from .test_utils import (
    TEST_REPOS_DIR,
    mock_shed_client,
)

from planemo import shed


def test_find_repository_id():
    with mock_shed_client() as tsi:
        repo_id = shed.find_repository_id(
            ctx=None,
            tsi=tsi,
            path=".",
            name="test_repo_1",
            owner="iuc",
        )
        assert repo_id == "r1"


def test_find_repository_id_missing():
    with mock_shed_client() as tsi:
        repo_id = shed.find_repository_id(
            ctx=None,
            tsi=tsi,
            path=".",
            name="test_repo_absent",
            owner="iuc",
            allow_none=True
        )
        assert repo_id is None


def test_find_repository_id_missing_exception():
    with mock_shed_client() as tsi:
        exception = None
        try:
            shed.find_repository_id(
                ctx=None,
                tsi=tsi,
                path=".",
                name="test_repo_absent",
                owner="iuc"
            )
        except Exception as e:
            exception = e
        assert exception is not None


def test_find_category_ids():
    with mock_shed_client() as tsi:
        category_ids = shed.find_category_ids(
            tsi,
            ["Text Manipulation"]
        )
        assert category_ids == ["c1"]


def test_create_simple():
    with mock_shed_client() as tsi:
        path = os.path.join(TEST_REPOS_DIR, "single_tool")
        repo_config = shed.shed_repo_config(path)
        create_response = shed.create_repository_for(
            None,
            tsi,
            "single_tool",
            repo_config,
        )
        assert "id" in create_response
        repo_id = shed.find_repository_id(
            ctx=None,
            tsi=tsi,
            path=".",
            name="single_tool",
            owner="iuc"
        )
        assert repo_id == create_response["id"]
