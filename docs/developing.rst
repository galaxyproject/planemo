==================
Release Checklist
==================

This page describes the process of releasing new versions of Planemo.

* Review ``git status`` for missing files.
* Verify the latest github workflows pass.
* Make sure that the `.venv` exists and contains the latest dev-requirements ``make setup-venv``
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
  ``.github/workflows/ci.yaml`` .
  If this didn't work, you can ``git checkout`` the tag and push to PyPI by
  executing ``make release-artifacts``
