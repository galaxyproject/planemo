==================
Release Checklist
==================

This page describes the process of releasing new versions of Planemo.

* Review ``git status`` for missing files.
* Verify the latest github workflows pass.
* Make sure that the `.venv` exists and contains the latest development requirements ``make setup-venv``
* Ensure the target release is set correctly in ``planemo/__init__.py`` (
  ``version`` will be a ``devN`` variant of target release).
* ``make add-history`` adds contributions to .dev0 version in HISTORY.rst
* ``make open-docs`` and review changelog.
* ``make clean && make lint``
* Review and commit outstanding changes.
* Update version and history, commit, add tag, mint a new version and push
  everything upstream with ``make release``
* The new tag should automatically push the new release to PyPI via the
  ``deploy`` job of the GitHub Actions workflow defined in
  ``.github/workflows/deploy.yaml`` .
  If this didn't work, you can ``git checkout`` the tag and push to PyPI by
  executing ``make release-artifacts``

Dependency management and distributions
=======================================

Declare runtime dependencies in ``pyproject.toml`` under ``project.dependencies``.
Development requirements live in ``dependency-groups``; ``tox.ini`` installs
them through ``dependency_groups``.

Type check with both mypy and Pyrefly using ``tox -e py310-mypy,py310-pyrefly``.
Their versions are pinned in the ``typecheck`` dependency group so new checker
releases can be adopted deliberately. ``pyrefly.toml`` was migrated from
``mypy.ini`` and uses the same fixture exclusion and missing-import policy;
keep those settings in sync when editing either configuration.

The Python CI matrix tests against ``uv.lock``. It exports all dependency groups
as constraints and sets ``PLANEMO_TEST_CONSTRAINTS`` so tox applies those pins
to both its test tools and Planemo's runtime dependencies. Galaxy instances
started by integration tests manage their dependencies independently. One
additional quick-test job resolves the latest allowed dependencies to check
compatibility with the library's unpinned requirements.

To run tox with the same constraints locally::

    uv export --locked --all-groups --no-emit-project --no-hashes --output-file /tmp/planemo-constraints.txt
    PLANEMO_TEST_CONSTRAINTS=/tmp/planemo-constraints.txt uv run --locked tox -e py310-unit-quick

Install a locked development environment with ``uv sync --locked``. To refresh
dependencies, run ``make update-dependencies`` and commit ``uv.lock``.
``make check-dependencies`` verifies the lock is current without changing it.
A weekly workflow (``.github/workflows/dependencies.yaml``) runs
``make update-dependencies`` and opens a pull request with the refreshed lock.

``make dist`` builds both ``planemo`` and ``planemo-cli`` and requires uv. The
CLI project is staged from the normal source distribution, with its
distribution name changed and its dependencies pinned to the runtime closure
exported from ``uv.lock``. Both wheels have identical code, assets,
versions, and command entry points. Each wheel is built from its own source
distribution; users do not need uv to rebuild either archive.

Both distributions use the Planemo release version and publish on the same tag.
Dependency-only updates use a new Planemo patch release. Release builds consume
the committed ``uv.lock`` and never refresh dependencies automatically. Configure
a PyPI trusted publisher for ``planemo-cli`` using this repository's existing
``deploy.yaml`` workflow before the first tagged release. Artifact installation
tests must pass before the publishing job runs.
