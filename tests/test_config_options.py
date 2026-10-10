"""Exercise configuration blending through real Click parsing and local profiles."""

import json
import os

import click
import pytest
import yaml
from click.testing import CliRunner
from gxjobconfinit.types import Runner

from planemo import options
from planemo.cli import (
    command_function,
    PlanemoCliContext,
)
from planemo.config import (
    OptionSource,
    planemo_option,
)


@pytest.fixture
def cli_context(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ctx = PlanemoCliContext()
    ctx.planemo_config = str(tmp_path / "config.yml")
    ctx.planemo_directory = str(tmp_path / "workspace")
    return ctx


def _invoke(ctx, decorator, args=(), env=None, config=None, profile=None):
    if config is not None:
        with open(ctx.planemo_config, "w") as config_file:
            yaml.safe_dump(config, config_file)
    if profile is not None:
        profile_directory = os.path.join(ctx.galaxy_profiles_directory, "example")
        os.makedirs(profile_directory)
        with open(os.path.join(profile_directory, "planemo_profile_options.json"), "w") as profile_file:
            json.dump({"engine": "docker_galaxy", **profile}, profile_file)
        args = ["--profile", "example", *args]

    @click.command()
    @options.profile_option()
    @decorator
    @command_function
    def command(ctx, **kwds):
        click.echo(json.dumps({"values": kwds, "sources": {k: v.name for k, v in ctx.option_source.items()}}))

    return CliRunner().invoke(command, args, obj=ctx, env=env)


@pytest.mark.parametrize("source", ["cli", "env", "global_config", "default", "profile"])
def test_path_conversion_for_all_sources(cli_context, tmp_path, source):
    (tmp_path / "target").mkdir()
    (tmp_path / "link").symlink_to(tmp_path / "target", target_is_directory=True)
    decorator = planemo_option(
        "--path",
        type=click.Path(resolve_path=True),
        use_global_config=True,
        use_env_var=True,
        default="link" if source == "default" else None,
    )
    result = _invoke(
        cli_context,
        decorator,
        args=["--path", "link"] if source == "cli" else [],
        env={"PLANEMO_PATH": "link"} if source == "env" else {},
        config={"default_path": "link"} if source == "global_config" else {},
        profile={"path": "link"} if source == "profile" else None,
    )
    assert result.exit_code == 0, result.output
    parsed = json.loads(result.output)
    assert parsed["values"]["path"] == str(tmp_path / "target")
    assert parsed["sources"]["path"] == ("cli" if source == "env" else source)


@pytest.mark.parametrize("source", ["global_config", "default", "profile"])
@pytest.mark.parametrize(
    "path_type,value",
    [
        (click.Path(exists=True), "missing"),
        (click.Path(dir_okay=False), "."),
        (click.Path(file_okay=False), "config.yml"),
    ],
)
def test_configured_paths_are_validated(cli_context, source, path_type, value):
    result = _invoke(
        cli_context,
        planemo_option(
            "--path", type=path_type, use_global_config=True, default=value if source == "default" else None
        ),
        config={"default_path": value} if source == "global_config" else {},
        profile={"path": value} if source == "profile" else None,
    )
    assert result.exit_code == 2, result.output
    assert "Invalid value for '--path'" in result.output


@pytest.mark.parametrize("value", ["first", ["first", "second"]])
@pytest.mark.parametrize("source", ["global_config", "default", "profile"])
def test_multiple_paths(cli_context, tmp_path, source, value):
    (tmp_path / "first").mkdir()
    (tmp_path / "second").mkdir()
    result = _invoke(
        cli_context,
        planemo_option(
            "--paths",
            type=click.Path(exists=True, resolve_path=True),
            multiple=True,
            use_global_config=True,
            default=value if source == "default" else None,
        ),
        config={"default_paths": value} if source == "global_config" else {},
        profile={"paths": value} if source == "profile" else None,
    )
    assert result.exit_code == 0, result.output
    expected = [value] if isinstance(value, str) else value
    parsed = json.loads(result.output)
    assert parsed["values"]["paths"] == [str(tmp_path / p) for p in expected]
    assert parsed["sources"]["paths"] == source


def test_extra_config_key_and_outer_callback(cli_context, tmp_path):
    def callback(ctx, param, value):
        assert value == str(tmp_path / "files")
        return value

    result = _invoke(
        cli_context,
        planemo_option(
            "--path",
            type=click.Path(resolve_path=True),
            use_global_config=True,
            extra_global_config_vars=["legacy_path"],
            callback=callback,
        ),
        config={"legacy_path": "files"},
    )
    assert result.exit_code == 0, result.output


@pytest.mark.parametrize("value,exit_code", [(["2", "3"], 0), (["2"], 2), (["2", "invalid"], 2)])
def test_nargs_conversion(cli_context, value, exit_code):
    result = _invoke(
        cli_context,
        planemo_option("--pair", type=int, nargs=2, use_global_config=True),
        config={"default_pair": value},
    )
    assert result.exit_code == exit_code, result.output
    if exit_code == 0:
        assert json.loads(result.output)["values"]["pair"] == [2, 3]
    else:
        assert "Invalid value for '--pair'" in result.output


def test_cli_and_environment_override_profile_and_config(cli_context, tmp_path):
    for args, env in [(["--path", "cli"], {}), ([], {"PLANEMO_PATH": "env"})]:
        ctx = PlanemoCliContext()
        ctx.planemo_config = cli_context.planemo_config
        ctx.planemo_directory = str(tmp_path / ("cli-workspace" if args else "env-workspace"))
        result = _invoke(
            ctx,
            planemo_option("--path", type=click.Path(resolve_path=True), use_global_config=True, use_env_var=True),
            args=args,
            env=env,
            config={"default_path": "global"},
            profile={"path": "profile"},
        )
        assert result.exit_code == 0, result.output
        parsed = json.loads(result.output)
        assert parsed["values"]["path"] == str(tmp_path / ("cli" if args else "env"))
        assert parsed["sources"]["path"] == OptionSource.cli.name


def test_runner_default_remains_an_enum(cli_context):
    @click.command()
    @options.runner_target_option()
    @command_function
    def command(ctx, runner):
        assert runner is Runner.LOCAL

    result = CliRunner().invoke(command, obj=cli_context)
    assert result.exit_code == 0, result.output


def test_unset_path_remains_none(cli_context):
    result = _invoke(cli_context, options.file_path_option(), config={})
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["values"]["file_path"] is None


def test_real_docker_extra_volume_option(cli_context, tmp_path):
    (tmp_path / "mount").mkdir()
    result = _invoke(cli_context, options.docker_extra_volume_option(), config={"default_docker_extra_volume": "mount"})
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["values"]["docker_extra_volume"] == [str(tmp_path / "mount")]


@pytest.mark.parametrize("source", ["cli", "env", "global_config", "default", "profile"])
def test_values_are_converted_only_once(cli_context, source):
    class OnceType(click.ParamType):
        name = "once"

        def convert(self, value, param, ctx):
            assert isinstance(value, str), "Value was converted twice"
            return int(value)

    result = _invoke(
        cli_context,
        planemo_option(
            "--value",
            type=OnceType(),
            use_global_config=True,
            use_env_var=True,
            default="2" if source == "default" else None,
        ),
        args=["--value", "2"] if source == "cli" else [],
        env={"PLANEMO_VALUE": "2"} if source == "env" else {},
        config={"default_value": "2"} if source == "global_config" else {},
        profile={"value": "2"} if source == "profile" else None,
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["values"]["value"] == 2


@pytest.mark.parametrize("value", [True, False])
def test_boolean_config_values(cli_context, value):
    result = _invoke(
        cli_context,
        planemo_option("--flag/--no_flag", default=True, use_global_config=True),
        config={"default_flag": value},
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["values"]["flag"] is value
