from .test_utils import CliShedTestCase
from planemo import shed

CS1_DESCRIPTION = "The tool Cat 1 from the cat tool suite."
CS2_DESCRIPTION = "The tool Cat 2 from the cat tool suite."


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

    def test_create_expansion_configured(self):
        with self._isolate_repo("multi_repos_flat_configured"):
            self._multi_repo_create_and_verify()

    def test_create_expansion_flag(self):
        with self._isolate_repo("multi_repos_flat_flag"):
            self._multi_repo_create_and_verify()

    def _multi_repo_create_and_verify(self):
        create_command = ["shed_create", "-r"]
        create_command.extend(self._shed_args())
        self._check_exit_code(create_command)
        cat1_repo = self._get_repo_info("cs-cat1")
        cat2_repo = self._get_repo_info("cs-cat2")
        assert cat1_repo is not None
        assert cat1_repo["synopsis"] == CS1_DESCRIPTION
        assert cat2_repo is not None
        assert cat2_repo["synopsis"] == CS2_DESCRIPTION

    def _get_repo_info(self, name):
        self._tsi.repositories.get_repositories()
        return shed.find_repository(self._tsi, "iuc", name)
