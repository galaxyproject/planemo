import os

from planemo import galaxy_config
from planemo import galaxy_run


def serve(ctx, path, **kwds):
    # TODO: Preceate a user.
    # TODO: Setup an admin user.
    # TODO: Pass through more parameters.
    # TODO: Populate test-data directory as FTP directory.
    daemon = kwds.get("daemon", False)
    if daemon:
        kwds["no_cleanup"] = True

    with galaxy_config.galaxy_config(ctx, path, **kwds) as config:
        # TODO: Allow running dockerized Galaxy here instead.
        run_script = os.path.join(config.galaxy_root, "run.sh")
        if daemon:
            run_script += " --daemon --wait"
            config.env["GALAXY_RUN_ALL"] = "1"
        else:
            run_script += " --server-name '%s' --reload" % config.server_name
        server_ini = os.path.join(config.config_directory, "galaxy.ini")
        config.env["GALAXY_CONFIG_FILE"] = server_ini
        cmds = [
            run_script,
        ]
        cmd = "; ".join(cmds)
        action = "Starting galaxy"
        galaxy_run.run_galaxy_command(
            ctx,
            cmd,
            config.env,
            action,
        )
        return config
