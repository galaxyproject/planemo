""" Planemo specific utilities for dealing with conda, extending Galaxy's
features with planemo specific idioms.
"""

from __future__ import absolute_import

from galaxy.tools.deps import conda_util
from planemo.io import shell
from planemo.io import info
from planemo import git
import os
from planemo.bioconductor_skeleton import write_recipe


from galaxy.tools.deps.requirements import parse_requirements_from_xml
from galaxy.tools.loader_directory import load_tool_elements_from_path


def build_conda_context(**kwds):
    """ Build a Galaxy CondaContext tailored to planemo use
    and common command-line arguments.
    """
    conda_prefix = kwds.get("conda_prefix", None)
    use_planemo_shell = kwds.get("use_planemo_shell_exec", True)
    ensure_channels = kwds.get("conda_ensure_channels", "")
    shell_exec = shell if use_planemo_shell else None
    return conda_util.CondaContext(conda_prefix=conda_prefix,
                                   ensure_channels=ensure_channels,
                                   shell_exec=shell_exec)


def collect_conda_targets(path, found_tool_callback=None, conda_context=None):
    conda_targets = []
    for (tool_path, tool_xml) in load_tool_elements_from_path(path):
        if found_tool_callback:
            found_tool_callback(tool_path)
        requirements, containers = parse_requirements_from_xml(tool_xml)
        conda_targets.extend(conda_util.requirements_to_conda_targets(requirements))
    return conda_targets


# Bioconda helper functions


def clone_bioconda_repo(path):
    """Clone bioconda repository in given path."""
    bioconda_repo = "git@github.com:bioconda/bioconda-recipes.git"
    git.clone(None, bioconda_repo, path)
    return "git clone of bioconda repo worked"


def write_bioconda_recipe(package_name, clone, update, bioconda_dir_path):
    """Make a bioconda recipe given the package name."""
    bioconda_recipe_path = os.path.join(bioconda_dir_path, "biconda-recipes")
    if clone and (not os.path.exists(bioconda_recipe_path)):
        clone_bioconda_repo(bioconda_recipe_path)
        info("bioconda-recipes cloned and writing to %s" % bioconda_dir_path)
    else:
        info("bioconda repository not cloned or already exists")
    # Check if package_name is in recipes
    presence = any(package_name in r for r, d, f in os.walk(bioconda_recipe_path))
    print("Is the package present? :", presence)
    if presence:
        # click.echo("Package already exisits in bioconda")
        if update:
            # click.echo("Package was updated in bioconda-repo")
            recipe_dir = os.path.join(bioconda_recipe_path, "recipes")
            write_recipe(package_name, recipe_dir, True)
    elif not presence:
        info("Package found in bioconda recipes")
        recipe_dir = os.path.join(bioconda_recipe_path, "recipes")
        write_recipe(package_name, recipe_dir, True)
    return
