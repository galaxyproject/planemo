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
* ``make add-history`` adds contributions to .dev0 version in HISTORY.rst
* ``make open-docs`` and review changelog.
* Review and commit outstanding changes.
* ``make check-release`` verifies the working tree is clean, the ``UPSTREAM``
  remote points to ``UPSTREAM_REPO`` (defaults: ``galaxyproject`` /
  ``galaxyproject/planemo``), ``HISTORY.rst`` has entries under the current
  ``.devN`` section, and ``make clean && make lint`` pass.
* Update version and history, commit, add tag, mint a new version and push
  everything upstream with ``make release``
* The new tag should automatically push the new release to PyPI via the
  workflow defined in ``.github/workflows/deploy.yaml``.
  If this didn't work, you can ``git checkout`` the tag and push to PyPI by
  executing ``make release-artifacts``
* Create the GitHub release for the tag; ``make release`` only pushes the
  tag. With `gh <https://cli.github.com/>`__::

      gh release create X.Y.Z --repo galaxyproject/planemo \
          --title X.Y.Z --verify-tag --generate-notes
