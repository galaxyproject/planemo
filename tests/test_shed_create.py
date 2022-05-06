from planemo import shed
from .test_utils import (
    CliShedTestCase,
    skip,
)

CS1_DESCRIPTION = "The tool Cat 1 from the cat tool suite."
CS2_DESCRIPTION = "The tool Cat 2 from the cat tool suite."
SUITE_DESCRIPTION = "A suite of Cat tools."
SUITE_DESCRIPTION_LONG = "A longer description of all the cat tools."


class ShedCreateTestCase(CliShedTestCase):
    def test_create_single(self):
        with self._isolate_repo("single_tool"):
            create_command = ["shed_create", "--skip_upload"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)

    @skip("Broken, see https://github.com/galaxyproject/planemo/issues/437")
    def test_create_wrong_owner(self):
        with self._isolate_repo("single_tool_other_owner"):
            create_command = ["shed_create", "--skip_upload"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command, exit_code=1)

    def test_create_multiple(self):
        with self._isolate_repo("multi_repos_nested"):
            create_command = ["shed_create", "--skip_upload", "-r"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)

    def test_create_expansion_configured(self):
        with self._isolate_repo("multi_repos_flat_configured"):
            self._multi_repo_create_and_verify()

    def test_create_expansion_flag(self):
        with self._isolate_repo("multi_repos_flat_flag"):
            self._multi_repo_create_and_verify()

    def test_create_expansion_flag_suite(self):
        with self._isolate_repo("multi_repos_flat_flag_suite"):
            self._multi_repo_create_and_verify()
            suite_repo = self._get_repo_info("suite_cat")
            assert suite_repo["synopsis"] == SUITE_DESCRIPTION
            assert suite_repo["description"] == SUITE_DESCRIPTION_LONG
            assert suite_repo["type"] == "repository_suite_definition"

    def _multi_repo_create_and_verify(self):
        create_command = ["shed_create", "--skip_upload", "-r"]
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
