import os

import yaml

from .test_utils import CliShedTestCase


class ShedInitTestCase(CliShedTestCase):

    def test_minimal(self):
        with self._isolate() as f:
            name = os.path.basename(os.path.abspath(f))
            self._check_exit_code(["shed_init", "--owner", "iuc"])
            shed_config_path = os.path.join(f, ".shed.yml")
            assert os.path.exists(shed_config_path)
            shed_config = yaml.load(open(shed_config_path, "r"))
            assert shed_config["name"] == name
            assert shed_config["owner"] == "iuc"
            assert shed_config["description"] == name
            assert len(shed_config["categories"]) == 0

    def test_more_options(self):
        with self._isolate() as f:
            repo_url = "https://github.com/galaxyproject/tools-devteam"
            init_command = [
                "shed_init",
                "--owner", "devteam",
                "--name", "samtools-filter",
                "--description", "A samtools repo",
                "--long_description", "A longer description.",
                "--remote_repository_url",
                repo_url,
                "--homepage_url", "https://example.com/",
                "--category", "SAM",
                "--category", "Sequence Analysis",
                "--category", "Statistics",
            ]
            self._check_exit_code(init_command)
            shed_config_path = os.path.join(f, ".shed.yml")
            assert os.path.exists(shed_config_path)
            shed_config = yaml.load(open(shed_config_path, "r"))
            assert shed_config["name"] == "samtools-filter"
            assert shed_config["owner"] == "devteam"
            assert shed_config["description"] == "A samtools repo"
            assert shed_config["long_description"] == "A longer description."
            assert shed_config["remote_repository_url"] == repo_url
            assert shed_config["homepage_url"] == "https://example.com/"

            categories = shed_config["categories"]
            assert len(categories) == 3
            assert "SAM" in categories
            assert "Statistics" in categories
