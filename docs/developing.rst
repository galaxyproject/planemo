==========
Developing
==========

This section contains documentation for the maintainers of planemo.

Release Checklist
-----------------

This release checklist is based on the `Pocoo Release Management Workflow
<http://www.pocoo.org/internal/release-management/>`_.

This assumes ``~/.pypirc`` file exists with the following fields (variations)
are fine.

::

    [distutils]
    index-servers =
        pypi
        test
    
    [pypi]
    username:<username>
    password:<password>
    
    [test]
    repository:https://testpypi.python.org/pypi
    username:<username>
    password:<password>


* Review ``git status`` for missing files.
* Verify the latest Travis CI builds pass.
* ``make open-docs`` and review changelog.
* ``make clean && make lint && make test``
* ``python scripts/commit_version.py <new_version>``
* ``make release``
    * Review `Test PyPI site <https://testpypi.python.org/pypi/planemo>`_
      for errors.
    * Test intall ``pip install -i https://testpypi.python.org/pypi planemo``.
* Update version info in ``planemo/__init__.py`` (n+1.dev0) and create new entry in HISTORY.rst.
* ``git add HISTORY.rst planemo/__init__.py; git commit -m "Start work on new version"``
* ``git push galaxyproject master``
* ``git push --tags galaxyproject``
* Update planemo homebrew recipe to new version.
