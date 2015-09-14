import os
import yaml

PLANEMO_CONFIG_ENV_PROP = "PLANEMO_GLOBAL_CONFIG_PATH"
DEFAULT_CONFIG = {
}


def get_default_callback(default, name=None):

    def callback(ctx, param, value):
        planemo_ctx = ctx.obj
        config_name = name
        if config_name is None:
            config_name = param.name

        return _default_option(planemo_ctx, config_name, value, default)

    return callback


# TODO: delete the following stuff.
def populate_kwds(ctx, options_hash, kwds):
    """ Iterate over default option hash and use value from
    ~/.planemo.yml configuration if set or else the specified
    default. Populate kwds with this configured default value.
    """
    for name, default in options_hash.items():
        _populate_default_output(ctx, name, kwds, default)


def _populate_default_output(ctx, kwd_key, kwds, default):
    kwd_value = kwds.get(kwd_key, None)
    if kwd_value is None:
        global_config = ctx.global_config
        global_config_key = "default_%s" % kwd_key
        if global_config_key in global_config:
            default_value = global_config[global_config_key]
        else:
            default_value = default

        if default_value:
            default_value = os.path.abspath(default_value)
        kwds[kwd_key] = default_value


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
