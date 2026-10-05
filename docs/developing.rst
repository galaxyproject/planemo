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
Development requirements live in ``dependency-groups``. ``requirements.txt`` and
``dev-requirements.txt`` are generated compatibility files for existing scripts.

Install a locked development environment with ``uv sync --locked``. To refresh
dependencies, run ``make update-dependencies`` and review ``uv.lock`` and the
three generated requirements files. Commit these files together. The CLI pins
are exported from the runtime dependency closure, without development groups,
and preserve Python/platform markers. ``make check-dependencies`` verifies the
lock's dependency declarations and generated files without changing them.

``make dist`` builds both ``planemo`` and ``planemo-cli``. The CLI project is
staged from the normal source distribution, with its distribution name and
runtime dependency metadata changed. Both wheels have identical code, assets,
versions, and command entry points. Each wheel is built from its own source
distribution; users do not need uv to rebuild either archive.

Both distributions use the Planemo release version and publish on the same tag.
Dependency-only updates use a new Planemo patch release. Release builds consume
the committed snapshot and never refresh dependencies automatically. Configure
a PyPI trusted publisher for ``planemo-cli`` using this repository's existing
``deploy.yaml`` workflow before the first tagged release. Artifact installation
tests must pass before the publishing job runs.
