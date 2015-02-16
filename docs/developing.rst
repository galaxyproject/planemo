==========
Developing
==========

This section contains documentation for developers of planemo.

Release Checklist
-----------------

This release checklist is based on the `Pocoo Release Management Workflow
<http://www.pocoo.org/internal/release-management/>`_.

 * Review ``git status`` for missing files.
 * Verify the latest Travis CI builds pass.
 * ``make lint && make test``
 * Update version info in ``planemo/__init__.py`` (drop ``-dev`` suffix).
 * Update release date and description in ``HISTORY.rst``.
 * ``make docs`` and review changelog.
 * ``git commit``
 * ``python2.7 setup.py bdist_wheel bdist_egg sdist upload``
 * ``python2.6 setup.py bdist_egg upload`` (TODO)
 * Check PyPI release page for obvious errors (https://pypi.python.org/pypi/planemo)
 * ``git tag <release>``
 * Update version info in ``planemo/__init__.py`` (n+1-dev).
 * ``git commit``
 * ``git push origin``
 * ``git push --tags origin``
 * Update planemo homebrew recipe to new version.

This requires ``~/.pypirc`` file.
