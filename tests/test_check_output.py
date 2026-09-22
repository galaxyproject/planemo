"""Unit tests for :mod:`planemo.test._check_output` collection handling."""

from planemo.test._check_output import (
    check_output,
    for_collections,
)


def test_for_collections_element_tests():
    assert for_collections({"element_tests": {"el1": {}}})


def test_for_collections_elements_alias():
    assert for_collections({"class": "Collection", "elements": {"el1": {}}})


def test_for_collections_class_only():
    assert for_collections({"class": "Collection", "element_count": 2})


def test_for_collections_file_output():
    assert not for_collections({"file": "expected.txt"})


def test_collection_type_mismatch_detected():
    problems = check_output(
        None,
        {"collection_type": "list:paired", "elements": []},
        {"name": "out", "class": "Collection", "collection_type": "list", "element_tests": {}},
    )
    assert len(problems) == 1
    assert "expected to be of type [list]" in problems[0]


def test_collection_type_match_no_problems():
    problems = check_output(
        None,
        {"collection_type": "list", "elements": []},
        {"name": "out", "class": "Collection", "collection_type": "list", "element_tests": {}},
    )
    assert problems == []
