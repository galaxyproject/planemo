import os

import yaml

from .test_utils import CliTestCase


class ConfigInitTestCase(CliTestCase):
    def test_config(self):
        assert not os.path.exists(self.planemo_yaml_path)
        with self._isolate():
            self._check_exit_code(["config_init"])
            assert os.path.exists(self.planemo_yaml_path)
            with open(self.planemo_yaml_path) as f:
                # Ensure it is valid YAML
                yaml.safe_load(f)

    def test_config_wont_overwrite(self):
        with self._isolate():
            # create the config
            self._check_exit_code(["config_init"])
            # make sure it isn't overwritten
            self._check_exit_code(["config_init"], exit_code=1)
