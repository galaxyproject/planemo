import os

from planemo.shed import diff
from .test_utils import TEST_DIR


def test_compare():
    local = os.path.join(TEST_DIR, "repository_dependencies.xml")
    shed = os.path.join(TEST_DIR, "repository_dependencies_shed.xml")
    assert not diff._shed_diff(local, shed)
