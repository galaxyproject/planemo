"""Module describing the planemo ``container_register`` command."""
import os

import click

from galaxy.tools.deps.mulled.util import image_name, quay_repository

from planemo import options
from planemo.cli import command_function
from planemo.git import add, branch, commit, push
from planemo.github_util import clone_fork_branch, get_repository_object, pull_request
from planemo.conda import collect_conda_target_lists, best_practice_search
from planemo.mulled import conda_to_mulled_targets

REGISTERY_TARGET_NAME = "multireqcontainers"
REGISTERY_REPOSITORY = "jmchilton/multireqcontainers"
DEFAULT_MESSAGE = "Add %s (generated with Planemo)."


@click.command('container_register')
@options.optional_tools_arg(multiple=True)
@options.recursive_option()
@options.mulled_namespace_option()
@click.option(
    "output_directory",
    "--output_directory",
    type=click.Path(
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
    ),
    default=None,
    help=("Container registration directory (defaults to ~/.planemo/multireqcontainers."),
)
@click.option(
    "-m",
    "--message",
    default=DEFAULT_MESSAGE,
    help="Commit and pull request message template for registration interactions."
)
@click.option(
    "--pull_request/--no_pull_request",
    is_flag=True,
    default=True,
    help="Fork and create a pull request against %s for these changes." % REGISTERY_REPOSITORY
)
@command_function
def cli(ctx, paths, **kwds):
    """Register multi-requirement containers as needed.

    BioContainers publishes all Bioconda packages automatically as individual
    container images. These however are not enough for tools with multiple
    best-practice requirements. Such requirements should be recorded and published
    so that a container can be created and registered for these tools.
    """
    message = kwds["message"]
    output_directory = kwds["output_directory"]
    do_pull_request = kwds.get("pull_request", True)
    if output_directory is None:
        output_directory = os.path.join(ctx.workspace, REGISTERY_TARGET_NAME)
        clone_fork_branch(
            ctx,
            "https://github.com/%s" % REGISTERY_REPOSITORY,
            output_directory,
            fork=do_pull_request,
        )
        pr_titles = [pr.title for pr in open_prs(ctx)]

    combinations_added = 0
    for conda_targets in collect_conda_target_lists(ctx, paths):
        mulled_targets = conda_to_mulled_targets(conda_targets)
        if len(mulled_targets) < 2:
            # Skip these for now, we will want to revisit this for conda-forge dependencies and such.
            continue

        best_practice_requirements = True
        for conda_target in conda_targets:
            best_hit, exact = best_practice_search(conda_target)
            if not best_hit or not exact:
                best_practice_requirements = False

        if not best_practice_requirements:
            continue

        name = image_name(mulled_targets)
        tag = "0"
        name_and_tag = "%s:%s" % (name, tag)
        target_filename = os.path.join(output_directory, "%s.tsv" % name_and_tag)
        if os.path.exists(target_filename):
            continue

        message = kwds["message"] % name
        namespace = kwds["mulled_namespace"]
        repo_data = quay_repository(namespace, name)
        if "tags" in repo_data:
            continue

        if do_pull_request:
            if any([name in t for t in pr_titles]):
                continue

        with open(target_filename, "w") as f:
            f.write(",".join(["%s=%s" % (t.package_name, t.version) for t in mulled_targets]))

        if do_pull_request:
            branch_name = name
            branch(ctx, output_directory, branch_name)
            add(ctx, output_directory, target_filename)
            commit(ctx, output_directory, message=message)
            push(ctx, output_directory, "origin", branch_name)
            pull_request(ctx, output_directory, message=message)
        combinations_added += 1


def open_prs(ctx):
    repo = get_repository_object(ctx, REGISTERY_REPOSITORY)
    prs = [pr for pr in repo.get_pulls()]
    return prs
