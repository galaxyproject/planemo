import os
import yaml

PLANEMO_CONFIG_ENV_PROP = "PLANEMO_GLOBAL_CONFIG_PATH"
DEFAULT_CONFIG = {
}


def get_default_callback(default, name=None, resolve_path=False):

    def callback(ctx, param, value):
        planemo_ctx = ctx.obj
        config_name = name
        if config_name is None:
            config_name = param.name

        result = _default_option(planemo_ctx, config_name, value, default)
        if resolve_path and result:
            result = os.path.abspath(result)
        return result

    return callback


def _default_option(ctx, kwd_key, kwd_value, default):
    value = kwd_value

    if value is None:
        global_config = ctx.global_config
        global_config_key = "default_%s" % kwd_key
        if global_config_key in global_config:
            default_value = global_config[global_config_key]
        else:
            default_value = default

        value = default_value

    return value


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


__all__ = [
    "global_config_path",
    "read_global_config",
    "get_default_callback",
]
