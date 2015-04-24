from .test_utils import CliShedTestCase


class ShedCreateTestCase(CliShedTestCase):

    def test_create_single(self):
        with self._isolate_repo("single_tool"):
            create_command = ["shed_create"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)

    def test_create_multiple(self):
        with self._isolate_repo("multi_repos_nested"):
            create_command = ["shed_create", "-r"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)
