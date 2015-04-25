""" Integration tests for shed_diff command.
"""

import os
from .test_utils import CliShedTestCase

DIFF_LINES = [
    "diff -r _local_/related_file _custom_shed_/related_file",
    "< A related non-tool file (modified).",
    "> A related non-tool file.",
]


class ShedUploadTestCase(CliShedTestCase):

    def test_shed_diff(self):
        with self._isolate_repo("single_tool") as f:
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            with open(os.path.join(f, "related_file"), "w") as r_f:
                r_f.write("A related non-tool file (modified).\n")

            diff_command = ["shed_diff", "-o", "diff"]
            diff_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(diff_command)
            diff_path = os.path.join(f, "diff")
            with open(diff_path, "r") as diff_f:
                diff = diff_f.read()
            for diff_line in DIFF_LINES:
                assert diff_line in diff
