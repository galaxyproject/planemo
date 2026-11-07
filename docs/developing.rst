==================
Release Checklist
==================

This page describes the process of releasing new versions of Planemo.

**Agents:** the ``ready-release`` skill
(``.agents/skills/ready-release/SKILL.md``) walks through this whole process
interactively. Pause for maintainer confirmation before ``make release``
(pushes to upstream and publishes to PyPI) and before publishing the GitHub
release (it notifies contributors).

* Verify the latest github workflows pass.
* Make sure that the `.venv` exists and contains the latest dev-requirements ``make setup-venv``
* Ensure the target release is set correctly in ``planemo/__init__.py`` (
  ``version`` will be a ``devN`` variant of target release). ``make release``
  mints the next patch version, so for a minor release edit ``__version__``
  and the ``.devN`` header in ``HISTORY.rst`` first.
* ``make add-history`` adds an entry for each PR merged since the last release
  to the .dev0 version in HISTORY.rst. It makes one GitHub API call per PR, so
  set ``GH_TOKEN`` (e.g. ``GH_TOKEN=$(gh auth token) make add-history``) to
  avoid the unauthenticated rate limit.
* ``make open-docs`` and review changelog; merge generated entries that
  duplicate hand-written ones.
* Review and commit outstanding changes.
* ``make check-release`` verifies the working tree is clean, the ``UPSTREAM``
  remote points to ``UPSTREAM_REPO`` (defaults: ``galaxyproject`` /
  ``galaxyproject/planemo``), you are on ``master`` and it contains
  ``UPSTREAM`` master, ``HISTORY.rst`` has entries under the current ``.devN``
  section, ``make clean && make lint`` pass, and ``make check-dist`` builds
  packages whose README/HISTORY pass ``twine check``.
* If ``UPSTREAM`` master moved since ``make add-history``, rebase and add
  entries for the newly merged PRs.
* ``make release`` reruns ``make check-release``, then updates version and
  history, commits, tags, checks the release packages, mints a new version and
  pushes master and the release tag (only) upstream.
* If ``make release`` fails after tagging: when only the push failed, rerun it
  with ``make push-release VERSION=<released version>``, plus the same
  ``UPSTREAM=...`` you gave ``make release`` if any (plain
  ``make push-release`` would use the next ``.dev0`` version). Otherwise delete
  the tag (``git tag -d <released version>``), drop the ``Version
  <released version>`` commit and the ``Starting work on ...`` commit if it
  was made, fix the problem and rerun ``make release``.
* The new tag should automatically push the new release to PyPI via the
  workflow defined in ``.github/workflows/deploy.yaml``.
  If this didn't work, you can ``git checkout`` the tag and push to PyPI by
  executing ``make release-artifacts``
* Create the GitHub release for the tag; ``make release`` only pushes the
  tag. With `gh <https://cli.github.com/>`__::

      gh release create X.Y.Z --repo galaxyproject/planemo \
          --title X.Y.Z --verify-tag --generate-notes
