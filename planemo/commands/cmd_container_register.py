"""Module describing the planemo ``container_register`` command."""
import os
import string

import click

from galaxy.tools.deps.mulled.util import conda_build_target_str, quay_repository, v2_image_name

from planemo import options
from planemo.cli import command_function
from planemo.conda import best_practice_search, collect_conda_target_lists_and_tool_paths
from planemo.git import add, branch, commit, push
from planemo.github_util import clone_fork_branch, get_repository_object, pull_request
from planemo.mulled import conda_to_mulled_targets

REGISTERY_TARGET_NAME = "multi-package-containers"
REGISTERY_TARGET_PATH = "combinations"
REGISTERY_REPOSITORY = "jmchilton/multi-package-containers"
DEFAULT_MESSAGE = "Add container $hash.\n**Hash**: $hash\n\n**Packages**:\n$packages\n\n**For** :\n$tools\n\nGenerated with Planemo."


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
    help=("Container registration directory (defaults to ~/.planemo/multi-package-containers."),
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
@click.option(
    "--force_push/--no_force_push",
    is_flag=True,
    default=False,
    help="Force push branch for pull request in case it already exists.",
)
@command_function
def cli(ctx, paths, **kwds):
    """Register multi-requirement containers as needed.

    BioContainers publishes all Bioconda packages automatically as individual
    container images. These however are not enough for tools with multiple
    best-practice requirements. Such requirements should be recorded and published
    so that a container can be created and registered for these tools.
    """
    registry_target = RegistryTarget(ctx, **kwds)

    combinations_added = 0
    conda_targets_list, tool_paths_list = collect_conda_target_lists_and_tool_paths(ctx, paths, recursive=kwds["recursive"])
    for conda_targets, tool_paths in zip(conda_targets_list, tool_paths_list):
        ctx.vlog("Handling conda_targets [%s]" % conda_targets)
        mulled_targets = conda_to_mulled_targets(conda_targets)
        mulled_targets_str = "- " + "\n- ".join(map(conda_build_target_str, mulled_targets))
        if len(mulled_targets) < 2:
            ctx.vlog("Skipping registeration, fewer than 2 targets discovered.")
            # Skip these for now, we will want to revisit this for conda-forge dependencies and such.
            continue

        best_practice_requirements = True
        for conda_target in conda_targets:
            best_hit, exact = best_practice_search(conda_target)
            if not best_hit or not exact:
                ctx.vlog("Target [%s] is not available in best practice channels - skipping" % conda_target)
                best_practice_requirements = False

        if not best_practice_requirements:
            continue

        name = v2_image_name(mulled_targets)
        tag = "0"
        name_and_tag = "%s-%s" % (name, tag)
        target_filename = os.path.join(registry_target.output_directory, "%s.tsv" % name_and_tag)
        ctx.vlog("Target filename for registeration is [%s]" % target_filename)
        if os.path.exists(target_filename):
            ctx.vlog("Target file already exists, skipping")
            continue

        namespace = kwds["mulled_namespace"]
        repo_data = quay_repository(namespace, name)
        if "tags" in repo_data:
            ctx.vlog("quay repository already exists, skipping")
            continue

        if registry_target.has_pull_request_for(name):
            ctx.vlog("Found matching open pull request for [%s], skipping" % name)
            continue

        registry_target.write_targets(ctx, target_filename, mulled_targets)
        tools_str = "\n".join(map(lambda p: "- " + os.path.basename(p), tool_paths))
        registry_target.handle_pull_request(ctx, name, target_filename, mulled_targets_str, tools_str, **kwds)
        combinations_added += 1


class RegistryTarget(object):
    """Abstraction around mulled container registery (both directory and Github repo)."""

    def __init__(self, ctx, **kwds):
        output_directory = kwds["output_directory"]
        pr_titles = []
        target_repository = None
        do_pull_request = kwds.get("pull_request", True)
        if output_directory is None:
            target_repository = os.path.join(ctx.workspace, REGISTERY_TARGET_NAME)
            output_directory = os.path.join(target_repository, REGISTERY_TARGET_PATH)
            clone_fork_branch(
                ctx,
                "https://github.com/%s" % REGISTERY_REPOSITORY,
                target_repository,
                fork=do_pull_request,
            )
            pr_titles = [pr.title for pr in open_prs(ctx)]

        self.do_pull_request = do_pull_request
        self.pr_titles = pr_titles
        self.output_directory = output_directory
        self.target_repository = target_repository

    def has_pull_request_for(self, name):
        has_pr = False
        if self.do_pull_request:
            if any([name in t for t in self.pr_titles]):
                has_pr = True

        return has_pr

    def handle_pull_request(self, ctx, name, target_filename, packages_str, tools_str, **kwds):
        if self.do_pull_request:
            message = kwds["message"]
            message = string.Template(message).safe_substitute({
                "hash": name,
                "packages": packages_str,
                "tools": tools_str,
            })
            branch_name = name.replace(":", "-")
            branch(ctx, self.target_repository, branch_name, from_branch="master")
            add(ctx, self.target_repository, target_filename)
            commit(ctx, self.target_repository, message=message)
            force_push = kwds.get("force_push", False)
            push(ctx, self.target_repository, os.environ.get("GITHUB_USER"), branch_name, force=force_push)
            pull_request(ctx, self.target_repository, message=message)

    def write_targets(self, ctx, target_filename, mulled_targets):
        with open(target_filename, "w") as f:
            contents = ",".join(["%s=%s" % (t.package_name, t.version) for t in mulled_targets])
            f.write(contents)
            ctx.vlog("Wrote requirements [%s] to file [%s]" % (contents, target_filename))


def open_prs(ctx):
    repo = get_repository_object(ctx, REGISTERY_REPOSITORY)
    prs = [pr for pr in repo.get_pulls()]
    return prs
