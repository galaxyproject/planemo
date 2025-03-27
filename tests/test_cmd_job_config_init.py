"""Tests for the ``job_config_init`` command."""

import os

import yaml

from .test_utils import CliTestCase


class CmdJobConfigInitTestCase(CliTestCase):
    def test_job_config_init_simple(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native"):
            init_cmd = ["job_config_init", "--tpv"]
            self._check_exit_code(init_cmd)
            assert os.path.exists("job_conf.yml")
            with open("job_conf.yml") as fh:
                config = yaml.safe_load(fh)
            assert config["execution"]["default"] == "tpv"

    def test_job_config_singularity(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native"):
            init_cmd = ["job_config_init", "--tpv", "--singularity"]
            self._check_exit_code(init_cmd)
            assert os.path.exists("job_conf.yml")
            with open("job_conf.yml") as fh:
                config = yaml.safe_load(fh)
            assert config["execution"]["default"] == "tpv"
            assert config["execution"]["environments"]["local"]["singularity_enabled"] is True

    def test_job_config_init_slurm(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native"):
            init_cmd = ["job_config_init", "--runner", "slurm"]
            self._check_exit_code(init_cmd)
            assert os.path.exists("job_conf.yml")
            with open("job_conf.yml") as fh:
                config = yaml.safe_load(fh)
            assert config["execution"]["default"] == "slurm"

    def test_job_config_init_25_0(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native"):
            init_cmd = ["job_config_init", "--tpv", "--runner", "drmaa", "--galaxy_version", "25.0"]
            self._check_exit_code(init_cmd)
            assert os.path.exists("job_conf.yml")
            with open("job_conf.yml") as fh:
                config = yaml.safe_load(fh)
            assert config["execution"]["default"] == "tpv"
            assert config["execution"]["environments"]["tpv"]["runner"] == "dynamic_tpv"

    def test_job_config_init_24_2(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native"):
            init_cmd = ["job_config_init", "--tpv", "--runner", "drmaa", "--galaxy_version", "24.2"]
            self._check_exit_code(init_cmd)
            assert os.path.exists("job_conf.yml")
            with open("job_conf.yml") as fh:
                config = yaml.safe_load(fh)
            assert config["execution"]["default"] == "tpv"
            assert config["execution"]["environments"]["tpv"]["runner"] == "dynamic"
