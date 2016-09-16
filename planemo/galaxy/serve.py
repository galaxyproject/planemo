"""Abstractions for serving out development Galaxy servers."""
from __future__ import print_function

import contextlib
import time

from planemo import io
from planemo import network_util

from .config import galaxy_config
from .run import (
    run_galaxy_command,
)


def serve(ctx, runnables=[], **kwds):
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

    with galaxy_config(ctx, runnables, **kwds) as config:
        cmd = config.startup_command(ctx, **kwds)
        action = "Starting galaxy"
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
        port = kwds.get("port", None)
        if port is None:
            port = network_util.get_free_port()

        ctx.vlog("Waiting for service on (%s, %s)" % (host, port))
        assert network_util.wait_net_service(host, port)
        time.sleep(.1)
        ctx.vlog("Waiting for service on (%s, %s)" % (host, port))
        assert network_util.wait_net_service(host, port)
        time.sleep(5)
        ctx.vlog("Waiting for service on (%s, %s)" % (host, port))
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
            if ctx.verbose:
                print("Galaxy Log:")
                print(config.log_contents)
            config.kill()
            if not kwds.get("no_cleanup", False):
                config.cleanup()

__all__ = ["serve", "serve_daemon", "shed_serve"]
