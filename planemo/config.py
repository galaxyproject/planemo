import os
import yaml

DEFAULT_CONFIG = {
}


def global_config_path():
    config_path = os.environ.get(
        "PLANEMO_GLOBAL_CONFIG_PATH",
        "~/.planemo.yml"
    )
    config_path = os.path.expanduser(config_path)
    return config_path


def read_global_config():
    config_path = global_config_path()
    if not os.path.exists(config_path):
        return DEFAULT_CONFIG

    with open(config_path) as f:
        return yaml.load(f)
