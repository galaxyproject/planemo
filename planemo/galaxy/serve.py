"""Abstractions for serving out development Galaxy servers."""
import contextlib
import os
import time

from .config import galaxy_config
from .run import (
    setup_venv,
    run_galaxy_command,
)
from planemo import io
from planemo import network_util


def serve(ctx, runnables=[], **kwds):
    """Serve a Galaxy instance with artifacts defined by paths."""
    daemon = kwds.get("daemon", False)
    if daemon:
        kwds["no_cleanup"] = True

    with galaxy_config(ctx, runnables, **kwds) as config:
        pid_file = config.pid_file
        # TODO: Allow running dockerized Galaxy here instead.
        setup_venv_command = setup_venv(ctx, kwds)
        run_script = os.path.join(config.galaxy_root, "run.sh")
        run_script += " $COMMON_STARTUP_ARGS"
        if daemon:
            run_script += " --pid-file '%s' --daemon" % pid_file
            config.env["GALAXY_RUN_ALL"] = "1"
        else:
            run_script += " --server-name '%s' --reload" % config.server_name
        server_ini = os.path.join(config.config_directory, "galaxy.ini")
        config.env["GALAXY_CONFIG_FILE"] = server_ini
        cd_to_galaxy_command = "cd %s" % config.galaxy_root
        cmd = io.shell_join(
            cd_to_galaxy_command,
            setup_venv_command,
            run_script,
        )
        action = "Starting galaxy"
        run_galaxy_command(
            ctx,
            cmd,
            config.env,
            action,
        )
        host = kwds.get("host", "127.0.0.1")
        port = kwds.get("port", None)
        if port is None:
            port = network_util.get_free_port()
        assert network_util.wait_net_service(host, port)
        time.sleep(.1)
        assert network_util.wait_net_service(host, port)
        time.sleep(1)
        assert network_util.wait_net_service(host, port)
        config.install_workflows()
        return config


@contextlib.contextmanager
def shed_serve(ctx, install_args_list, **kwds):
    """Serve a daemon instance of Galaxy with specified repositories installed."""
    with serve_daemon(ctx, **kwds) as config:
        install_deps = not kwds.get("skip_dependencies", False)
        io.info("Installing repositories - this may take some time...")
        for install_args in install_args_list:
            install_args["install_tool_dependencies"] = install_deps
            install_args["install_repository_dependencies"] = True
            install_args["new_tool_panel_section_label"] = "Shed Installs"
            config.install_repo(
                **install_args
            )
        config.wait_for_all_installed()
        yield config


@contextlib.contextmanager
def serve_daemon(ctx, runnables=[], **kwds):
    """Serve a daemonized Galaxy instance with artifacts defined by paths."""
    config = None
    try:
        kwds["daemon"] = True
        config = serve(ctx, runnables, **kwds)
        yield config
    finally:
        if config:
            config.kill()
            if not kwds.get("no_cleanup", False):
                config.cleanup()

__all__ = ["serve", "serve_daemon", "shed_serve"]
