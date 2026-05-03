"""Structured metadata for Planemo's Click command tree."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import click

from planemo import __version__
from planemo.cli import (
    COMMAND_ALIASES,
    INTERNAL_COMMANDS,
    list_cmds,
    name_to_command,
)

SCHEMA_VERSION = "0.1"


def iter_command_names(include_internal: bool = False) -> list[str]:
    """Return command names included in public CLI metadata."""
    commands = list_cmds()
    if not include_internal:
        commands = [command for command in commands if command not in INTERNAL_COMMANDS]
    return commands


def load_command_metadata(command_name: str) -> dict[str, Any]:
    """Load metadata for one Planemo command."""
    return serialize_click_command(command_name, name_to_command(command_name))


def load_planemo_metadata(include_internal: bool = False) -> dict[str, Any]:
    """Load metadata for the Planemo command line interface."""
    return {
        "schema_version": SCHEMA_VERSION,
        "program": "planemo",
        "planemo_version": __version__,
        "commands": [load_command_metadata(command_name) for command_name in iter_command_names(include_internal)],
        "aliases": dict(sorted(COMMAND_ALIASES.items())),
    }


def serialize_click_command(command_name: str, command: click.Command) -> dict[str, Any]:
    """Serialize one Click command into JSON-compatible metadata."""
    context = click.Context(command, info_name=command_name)
    return {
        "name": command_name,
        "module": f"planemo.commands.cmd_{command_name}",
        "help": command.help,
        "short_help": command.short_help,
        "usage": command.get_usage(context).removeprefix("Usage: "),
        "internal": command_name in INTERNAL_COMMANDS,
        "hidden": getattr(command, "hidden", False),
        "params": [serialize_click_param(param) for param in command.params],
    }


def serialize_click_param(param: click.Parameter) -> dict[str, Any]:
    """Serialize one Click parameter into JSON-compatible metadata."""
    planemo_config = getattr(param, "planemo_config", {})
    metadata = {
        "kind": "option" if isinstance(param, click.Option) else "argument",
        "name": param.name,
        "opts": list(getattr(param, "opts", [])),
        "secondary_opts": list(getattr(param, "secondary_opts", [])),
        "help": getattr(param, "help", None),
        "required": param.required,
        "human_readable_name": param.human_readable_name,
        "multiple": param.multiple,
        "nargs": param.nargs,
        "type": serialize_click_type(param.type),
        "default": _json_value(planemo_config.get("declared_default", getattr(param, "default", None))),
        "is_flag": getattr(param, "is_flag", False),
        "flag_value": _json_value(getattr(param, "flag_value", None)),
        "envvar": _json_value(getattr(param, "envvar", None)),
        "hidden": getattr(param, "hidden", False),
        "prompt": getattr(param, "prompt", None),
        "planemo_config": _json_value(planemo_config),
    }
    if isinstance(param, click.Option):
        metadata["is_bool_flag"] = param.is_bool_flag
    return metadata


def serialize_click_type(param_type: click.ParamType) -> dict[str, Any]:
    """Serialize Click parameter type details where Click exposes them."""
    metadata: dict[str, Any] = {"name": param_type.name}
    if isinstance(param_type, click.Choice):
        metadata["choices"] = list(param_type.choices)
        metadata["case_sensitive"] = param_type.case_sensitive
    elif isinstance(param_type, click.IntRange):
        metadata.update(
            {
                "min": param_type.min,
                "max": param_type.max,
                "min_open": param_type.min_open,
                "max_open": param_type.max_open,
            }
        )
    elif isinstance(param_type, click.FloatRange):
        metadata.update(
            {
                "min": param_type.min,
                "max": param_type.max,
                "min_open": param_type.min_open,
                "max_open": param_type.max_open,
            }
        )
    elif isinstance(param_type, click.Path):
        for attr in ("exists", "file_okay", "dir_okay", "writable", "readable", "resolve_path", "allow_dash"):
            if hasattr(param_type, attr):
                metadata[attr] = getattr(param_type, attr)
    return _json_value(metadata)


def _json_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_value(nested_value) for key, nested_value in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_json_value(nested_value) for nested_value in value]
    return repr(value)
