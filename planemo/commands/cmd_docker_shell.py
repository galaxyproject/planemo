"""The docker_shell command takes the path to a single tool file and returns
a new command to launch a new shell into the corresponding Docker container
for that tool.

Galaxy will treat the identifier in the tool XML file as a published Docker
image - during development such the image an image may not yet exist or may
be not updated. So during development the --from_tag argument to used to treat
the identifier as a locally cached tag. (Tip: Use the docker_build to populate
such a tag from a Dockerfile located in the tool's directory.
"""
from __future__ import print_function
import click

from planemo.cli import pass_context
from planemo import options

from galaxy.tools.loader import load_tool
from galaxy.tools.deps import docker_util
from galaxy.tools.deps import dockerfiles
from galaxy.tools.deps.requirements import parse_requirements_from_xml


@click.command('docker_shell')
@options.required_tool_arg()
@click.option(
    '--from_tag',
    is_flag=True,
    help=(
        "Treat the tool's Docker container identifier as a locally cached tag."
    )
)
@click.option(
    '--shell',
    default="/bin/bash",
    help="Shell to launch in container (defaults to /bin/bash)."
)
@options.docker_cmd_option()
@options.docker_sudo_option()
@options.docker_sudo_cmd_option()
@options.docker_host_option()
@pass_context
def cli(ctx, path, **kwds):
    """Launch a shell in the Docker container referenced by the specified
    tool. Prints a command to do this the way Galaxy would in job files it
    generates - so be sure to wrap this in $(...) to launch the subshell.::

        % $(planemo docker_shell bowtie2.xml)
        ...
        root@b8754062f875:/#

    """
    tool_xml = load_tool(path)
    requirements, containers = parse_requirements_from_xml(tool_xml)
    identifier = None
    for container in containers:
        if container.type == "docker":
            identifier = container.identifier

    if kwds["from_tag"]:
        identifier = "-t %s" % identifier

    script = docker_util.build_docker_run_command(
        "/bin/bash",
        identifier,
        interactive=True,
        **dockerfiles.docker_host_args(**kwds)
    )
    print(script)
