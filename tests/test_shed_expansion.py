import os

import yaml

from planemo import shed
from .test_utils import TempDirectoryTestCase


class ShedExpansionTestCase(TempDirectoryTestCase):
    def test_invalid_name_with_multiple_repos(self):
        self._make_shed_yml("repo1", name="repo1")
        self._make_shed_yml("repo2", name="repo2")

        with self.assertRaises(Exception):
            self._repos(name="repo1", recursive=True)

    def test_invalid_name_conflict(self):
        self._make_shed_yml("repo2", name="repo2")

        with self.assertRaises(Exception):
            self._repos(name="repo1", recursive=True)

    def _make_shed_yml(self, path, **kwds):
        shed_dir = os.path.join(self.temp_directory, path)
        os.makedirs(shed_dir)
        shed_yml = os.path.join(shed_dir, ".shed.yml")
        with open(shed_yml, "w") as f:
            yaml.dump(kwds, f)

    def _repos(self, **kwds):
        return shed.build_effective_repositories(self.temp_directory, **kwds)
