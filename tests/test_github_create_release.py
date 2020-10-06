import os
import uuid

import pytest

from planemo import git
from planemo.github_util import (
    add_dir_contents_to_repo,
    get_or_create_repository,
)
from planemo.io import temp_directory
from .test_utils import (
    test_context,
)


def test_get_or_create_repo_with_existing_repo():
    ctx = test_context()
    ctx._global_config = {"github": {"access_token": 'ABCDEFG'}}
    repository_path = get_or_create_repository(ctx, owner='galaxyproject', repo='planemo', dry_run=False)
    assert os.path.exists(os.path.join(repository_path, 'README.rst'))


def test_get_or_create_repo_with_new_repo():
    ctx = test_context()
    ctx._global_config = {"github": {"access_token": 'ABCDEFG'}}
    with pytest.raises(RuntimeError) as excinfo:
        # Token isn't valid, so this errors out while running gh create
        get_or_create_repository(ctx, owner=str(uuid.uuid4()), repo='some-repo', dry_run=False)
    assert "Problem executing commands" in str(excinfo.value)
    assert "gh repo create -y" in str(excinfo.value)


def test_add_dir_contents_to_repo():
    ctx = test_context()
    ctx._global_config = {"github": {"access_token": 'ABCDEFG'}}
    with temp_directory() as test_dir, temp_directory() as repo_dir:
        with open(os.path.join(test_dir, 'Readme.md'), 'w') as readme:
            readme.write('#Very important!')
        git.init(ctx, repo_path=repo_dir)
        with pytest.raises(RuntimeError) as excinfo:
            # Can't push without remote
            add_dir_contents_to_repo(
                ctx,
                from_dir=test_dir,
                target_dir="workflows/my-cool-workflow",
                target_repository_path=repo_dir,
                version=1.0,
                notes='The big release!',
                dry_run=False
            )
        assert "Problem executing commands git push" in str(excinfo.value)


def test_add_dir_contents_to_repo_dry_run():
    ctx = test_context()
    ctx._global_config = {"github": {"access_token": 'ABCDEFG'}}
    with temp_directory() as test_dir, temp_directory() as repo_dir:
        with open(os.path.join(test_dir, 'Readme.md'), 'w') as readme:
            readme.write('#Very important!')
        git.init(ctx, repo_path=repo_dir)
        add_dir_contents_to_repo(
            ctx,
            from_dir=test_dir,
            target_dir="workflows/my-cool-workflow",
            target_repository_path=repo_dir,
            version=1.0,
            notes='The big release!',
            dry_run=True
        )
