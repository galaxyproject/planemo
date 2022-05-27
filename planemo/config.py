"""Module defines abstractions for configuring Planemo."""

import os
from enum import Enum
from typing import (
    Any,
    Dict,
)

import click
import yaml

PLANEMO_CONFIG_ENV_PROP = "PLANEMO_GLOBAL_CONFIG_PATH"
DEFAULT_CONFIG = {}  # type: Dict[str, Any]

VALUE_UNSET = object()


OptionSource = Enum("OptionSource", "cli profile global_config default")


def _default_callback(
    default,
    use_global_config=False,
    resolve_path=False,
    extra_global_config_vars=[],
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
                extra_global_config_vars=extra_global_config_vars,
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


def _find_default(ctx, param, use_global_config, extra_global_config_vars):
    if use_global_config:
        global_config = ctx.global_config
        global_config_keys = ["default_%s" % param.name] + extra_global_config_vars
        for global_config_key in global_config_keys:
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
    use_env_var = kwargs.pop("use_env_var", False)
    extra_global_config_vars = kwargs.pop("extra_global_config_vars", [])

    default_specified = "default" in kwargs
    default = None
    if default_specified:
        default = kwargs.pop("default")

    if default_specified or use_global_config or use_env_var:
        outer_callback = kwargs.pop("callback", None)

        def callback(ctx, param, value):
            resolve_path = option_type and getattr(option_type, "resolve_path", False)
            result = _default_callback(
                default,
                use_global_config=use_global_config,
                extra_global_config_vars=extra_global_config_vars,
                resolve_path=resolve_path,
            )(ctx, param, value)

            if outer_callback is not None:
                result = outer_callback(ctx, param, result)

            return result

        kwargs["callback"] = callback

    if default_specified:
        kwargs["default"] = None

    if use_env_var:
        name = None
        for arg in args:
            if arg.startswith("--"):
                name = arg[len("--") :]
        assert name
        kwargs["envvar"] = "PLANEMO_%s" % name.upper()

    option = click.option(*args, **kwargs)
    return option


def global_config_path(config_path=None):
    if not config_path:
        config_path = os.environ.get(PLANEMO_CONFIG_ENV_PROP, "~/.planemo.yml")
        config_path = os.path.expanduser(config_path)
    return config_path


def read_global_config(config_path: str):
    config_path = global_config_path(config_path)
    if not os.path.exists(config_path):
        return DEFAULT_CONFIG

    with open(config_path) as f:
        return yaml.safe_load(f)


__all__ = (
    "global_config_path",
    "read_global_config",
    "planemo_option",
)
