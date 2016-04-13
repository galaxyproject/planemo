"""Module defines abstractions for configuring Planemo."""
import os
import yaml

import aenum
import click

PLANEMO_CONFIG_ENV_PROP = "PLANEMO_GLOBAL_CONFIG_PATH"
DEFAULT_CONFIG = {
}

VALUE_UNSET = object()


OptionSource = aenum.Enum(
    "OptionSource", 'cli profile global_config default'
)


def _default_callback(
    default, use_global_config=False, resolve_path=False,
):

    def callback(ctx, param, value):
        planemo_ctx = ctx.obj
        param_name = param.name
        if value is not None:
            result = value
            option_source = OptionSource.cli
        else:
            result, option_source = _find_default(
                planemo_ctx,
                param,
                use_global_config=use_global_config,
            )

        if result is VALUE_UNSET:
            result = default
            option_source = OptionSource.default

        assert option_source is not None
        assert result is not VALUE_UNSET

        if resolve_path and result is not None:
            result = os.path.abspath(result)

        planemo_ctx.set_option_source(param_name, option_source)
        return result

    return callback


def _find_default(ctx, param, use_global_config):
    if use_global_config:
        global_config = ctx.global_config
        global_config_key = "default_%s" % param.name
        if global_config_key in global_config:
            default_value = global_config[global_config_key]
            return default_value, OptionSource.global_config

    return VALUE_UNSET, None


def planemo_option(*args, **kwargs):
    """Extend ``click.option`` with planemo-config aware configuration.

    This extends click.option to use a callback when assigning default
    values, add ``use_global_config`` keyword argument to allow reading
    defaults from ~/.planemo.yml, and tracks how parameters are specified
    using the Planemo Context object.
    """
    option_type = kwargs.get("type", None)
    use_global_config = kwargs.pop("use_global_config", False)

    if "default" in kwargs:
        default = kwargs.pop("default")
        outer_callback = kwargs.pop("callback", None)

        def callback(ctx, param, value):
            resolve_path = option_type and getattr(option_type, "resolve_path", False)
            result = _default_callback(
                default,
                use_global_config=use_global_config,
                resolve_path=resolve_path,
            )(ctx, param, value)

            if outer_callback is not None:
                result = outer_callback(ctx, param, result)

            return result

        kwargs["callback"] = callback
        kwargs["default"] = None

    return click.option(*args, **kwargs)


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
    "planemo_option",
]
