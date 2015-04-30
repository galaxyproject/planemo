import os
import yaml

PLANEMO_CONFIG_ENV_PROP = "PLANEMO_GLOBAL_CONFIG_PATH"
DEFAULT_CONFIG = {
}


def global_config_path(config_path=None):
    if not config_path:
        config_path = os.environ.get(
            PLANEMO_CONFIG_ENV_PROP,
            "~/.planemo.yml"
        )
        config_path = os.path.expanduser(config_path)
    return config_path


def read_global_config(config_path):
    config_path = global_config_path(config_path)
    if not os.path.exists(config_path):
        return DEFAULT_CONFIG

    with open(config_path) as f:
        return yaml.load(f)
