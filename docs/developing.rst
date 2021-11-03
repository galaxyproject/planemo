==================
Release Checklist
==================

This page describes the process of releasing new versions of Planemo.

* Update dependencies by running ``make dependencies-update`` and commit changes.
* Review ``git status`` for missing files.
* Open Pull Request to verify tests pass with new dependencies.
* Update ``HISTORY.rst`` with the help of ``scripts/bootstrap_history.py``
* ``make open-docs`` and review changelog.
* Ensure the target release is set correctly in ``planemo/__init__.py`` (
  ``version`` will be a ``devN`` variant of target release) and pyproject.toml.
* ``make clean && make lint && make test``
* Commit outstanding changes.
* Update version and history, commit, add tag, mint a new version and push
  everything upstream with ``make release``
* The new tag should automatically push the new release to PyPI via the
  ``deploy`` job of the GitHub Actions workflow defined in
  ``.github/workflows/ci.yaml`` .
  If this didn't work, you can ``git checkout`` the tag and push to PyPI by
  executing ``make release-artifacts``
