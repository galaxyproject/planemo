"""Abstractions for serving out development Galaxy servers."""
from __future__ import print_function

import contextlib
import os
import time

from planemo import (
    io,
    network_util
)
from .config import galaxy_config
from .ephemeris_sleep import sleep
from .run import (
    run_galaxy_command,
)
INSTALLING_MESSAGE = "Installing repositories - this may take some time..."


def serve(ctx, runnables=None, **kwds):
    if runnables is None:
        runnables = []
    """Serve a Galaxy instance with artifacts defined by paths."""
    try:
        return _serve(ctx, runnables, **kwds)
    except Exception as e:
        ctx.vlog("Problem serving Galaxy", exception=e)
        raise


def _serve(ctx, runnables, **kwds):
    engine = kwds.get("engine", "galaxy")
    if engine == "docker_galaxy":
        kwds["dockerize"] = True

    daemon = kwds.get("daemon", False)
    if daemon:
        kwds["no_cleanup"] = True

    port = kwds.get("port", None)
    if port is None:
        port = network_util.get_free_port()
        kwds["port"] = port

    with galaxy_config(ctx, runnables, **kwds) as config:
        cmd = config.startup_command(ctx, **kwds)
        action = "Starting Galaxy"
        exit_code = run_galaxy_command(
            ctx,
            cmd,
            config.env,
            action,
        )
        if exit_code:
            message = "Problem running Galaxy command [%s]." % config.log_contents
            io.warn(message)
            raise Exception(message)
        host = kwds.get("host", "127.0.0.1")

        timeout = 500
        galaxy_url = "http://%s:%s" % (host, port)
        galaxy_alive = sleep(galaxy_url, verbose=ctx.verbose, timeout=timeout)
        if not galaxy_alive:
            raise Exception("Attempted to serve Galaxy at %s, but it failed to start in %d seconds." % (galaxy_url, timeout))
        config.install_workflows()
        if kwds.get("pid_file"):
            real_pid_file = config.pid_file
            if os.path.exists(config.pid_file):
                os.symlink(real_pid_file, kwds["pid_file"])
            else:
                io.warn("Can't find Galaxy pid file [%s] to link" % real_pid_file)
        return config


@contextlib.contextmanager
def shed_serve(ctx, install_args_list, **kwds):
    """Serve a daemon instance of Galaxy with specified repositories installed."""
    with serve_daemon(ctx, **kwds) as config:
        install_deps = not kwds.get("skip_dependencies", False)
        print(INSTALLING_MESSAGE)
        io.info(INSTALLING_MESSAGE)
        for install_args in install_args_list:
            install_args["install_tool_dependencies"] = install_deps
            install_args["install_repository_dependencies"] = True
            install_args["new_tool_panel_section_label"] = "Shed Installs"
            config.install_repo(
                **install_args
            )
        try:
            config.wait_for_all_installed()
        except Exception:
            if ctx.verbose:
                print("Failed to install tool repositories, Galaxy log:")
                print(config.log_contents)
                print("Galaxy root:")
                io.shell(['ls', config.galaxy_root])
            raise
        yield config


@contextlib.contextmanager
def serve_daemon(ctx, runnables=None, **kwds):
    """Serve a daemonized Galaxy instance with artifacts defined by paths."""
    if runnables is None:
        runnables = []
    config = None
    try:
        kwds["daemon"] = True
        config = serve(ctx, runnables, **kwds)
        yield config
    finally:
        if config:
            if ctx.verbose:
                print("Galaxy Log:")
                print(config.log_contents)
            config.kill()
            if not kwds.get("no_cleanup", False):
                config.cleanup()


def sleep_for_serve():
    # This is bad, do something better...
    time.sleep(1000000)


__all__ = (
    "serve",
    "serve_daemon",
    "shed_serve",
)
