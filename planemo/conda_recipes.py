"""Planemo specific utilities for dealing with conda recipe generation.

The extend Galaxy/galaxy-lib's features with planemo specific idioms.
"""

from __future__ import absolute_import

import os

from planemo import git
from planemo.bioconda_scripts import bioconductor_skeleton
from planemo.io import info


# from galaxy.tools.deps import conda_util

# from planemo.io import shell
# from planemo.tools import yield_tool_sources_on_paths


# Bioconda helper functions


def clone_bioconda_repo(path):
    """Clone bioconda repository in given path."""
    bioconda_repo = "git@github.com:bioconda/bioconda-recipes.git"
    git.clone(None, bioconda_repo, path)
    return "git clone of bioconda repo worked"


def write_bioconda_recipe(package_name, clone, update, bioconda_dir_path=None):
    """Make a bioconda recipe given the package name.

    clone: y/N , clone the whole bioconda repository and create recipe inside
    repository.

    update: The update feature differs from the one in bioconda, as it
    updates the specific package, as opposed to the every package in the
    biocoda repository.
    """
    # set bioconda path
    if bioconda_dir_path is None:
        bioconda_recipe_path = os.path.join(os.path.expanduser("~"), "bioconda-recipes")
    else:
        bioconda_recipe_path = os.path.join(bioconda_dir_path, "bioconda-recipes")

    # Clone
    if clone and (not os.path.exists(bioconda_recipe_path)):
        clone_bioconda_repo(bioconda_recipe_path)
        info("bioconda-recipes cloned and writing to %s" % bioconda_dir_path)
    else:
        info("Bioconda repository not cloned or already exists")

    # Check if package_name is in recipes
    presence = any(package_name in r for r, d, f in os.walk(bioconda_recipe_path))
    if presence:
        info("Package already exists in bioconda")
        if update:
            info("Package will be updated")
            recipe_dir = os.path.join(bioconda_recipe_path, "recipes")
            bioconductor_skeleton.write_recipe(package_name, recipe_dir, True)
    elif not presence:
        info("Package found in bioconda recipes")
        recipe_dir = os.path.join(bioconda_recipe_path, "recipes")
        bioconductor_skeleton.write_recipe(package_name, recipe_dir, True)
    return


__all__ = (
    "clone_bioconda_repo",
    "write_bioconda_recipe",
)
