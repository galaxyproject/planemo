import os

import yaml

from .test_utils import (
    CliShedTestCase,
    TEST_REPOS_DIR,
)

TEST_WORKFLOW_PATH = os.path.join(TEST_REPOS_DIR, "workflow_1", "blast_top_hit_species.ga")


class ShedInitTestCase(CliShedTestCase):
    def test_minimal(self):
        with self._isolate() as f:
            self._check_exit_code(["shed_init", "--owner", "iuc", "--name", "samtools_filter"])
            shed_config_path = os.path.join(f, ".shed.yml")
            assert os.path.exists(shed_config_path)
            with open(shed_config_path) as fh:
                shed_config = yaml.safe_load(fh)
            assert shed_config["name"] == "samtools_filter"
            assert shed_config["owner"] == "iuc"
            assert shed_config["description"] == "samtools_filter"
            assert len(shed_config["categories"]) == 0

    def test_more_options(self):
        with self._isolate() as f:
            repo_url = "https://github.com/galaxyproject/tools-devteam"
            init_command = [
                "shed_init",
                "--owner",
                "devteam",
                "--name",
                "samtools_filter",
                "--description",
                "A samtools repo",
                "--long_description",
                "A longer description.",
                "--remote_repository_url",
                repo_url,
                "--homepage_url",
                "https://example.com/",
                "--category",
                "SAM",
                "--category",
                "Sequence Analysis",
                "--category",
                "Statistics",
            ]
            self._check_exit_code(init_command)
            shed_config_path = os.path.join(f, ".shed.yml")
            assert os.path.exists(shed_config_path)
            with open(shed_config_path) as fh:
                shed_config = yaml.safe_load(fh)
            assert shed_config["name"] == "samtools_filter"
            assert shed_config["owner"] == "devteam"
            assert shed_config["description"] == "A samtools repo"
            assert shed_config["long_description"] == "A longer description."
            assert shed_config["remote_repository_url"] == repo_url
            assert shed_config["homepage_url"] == "https://example.com/"

            categories = shed_config["categories"]
            assert len(categories) == 3
            assert "SAM" in categories
            assert "Statistics" in categories

    def test_from_workflow(self):
        with self._isolate() as f:
            init_command = ["shed_init", "--owner", "iuc"]
            init_command += ["--category", "Sequence Analysis"]
            init_command += ["--name", "blasthits"]
            init_command += ["--from_workflow", TEST_WORKFLOW_PATH]
            self._check_exit_code(init_command)

            repo_deps_path = os.path.join(f, "repository_dependencies.xml")
            wf_path = os.path.join(f, "blast_top_hit_species.ga")
            assert os.path.exists(repo_deps_path)
            assert os.path.exists(wf_path)

            # lint repository as a way of verifying repository_dependencies
            self._check_exit_code(["shed_lint"])

    def test_bad_name(self):
        with self._isolate():
            # has an invalid -
            self._check_exit_code(["shed_init", "--owner", "iuc", "--name", "samtools-filter"], exit_code=1)

    def test_bad_owner(self):
        with self._isolate():
            self._check_exit_code(["shed_init", "--owner", "IuC", "--name", "samtools_filter"], exit_code=1)
