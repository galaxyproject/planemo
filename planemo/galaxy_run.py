from planemo.io import info
from galaxy.tools.deps.commands import shell


# Activate galaxy's virtualenv if present (needed for tests say but not for
# server because run.sh does this).
ACTIVATE_COMMAND = "[ -e .venv ] && . .venv/bin/activate"

# TODO: Mac-y curl variant of this.
DOWNLOAD_GALAXY = (
    "wget https://codeload.github.com/jmchilton/galaxy-central/tar.gz/master"
)


def run_galaxy_command(ctx, command, env, action):
    message = "%s with command [%s]" % (action, command)
    info(message)
    ctx.vlog("With environment variables:")
    ctx.vlog("============================")
    for key, value in env.items():
        ctx.vlog('%s="%s"' % (key, value))
    ctx.vlog("============================")
    return shell(command, env=env)
