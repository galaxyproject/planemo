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
* ``make clean && make lint && make test``
* Update version info in ``planemo/__init__.py`` (drop ``.dev0`` suffix).
* Update release date and description in ``HISTORY.rst``.
* ``make docs`` and review changelog.
* ``git add HISTORY.rst planemo/__init__.py; git commit -m "Version <version>"``
* ``make release``
    * Review `Test PyPI site <https://testpypi.python.org/pypi/planemo>`_
      for errors.
    * Test intall ``pip install -i https://testpypi.python.org/pypi planemo``.
* ``git tag <release>``
* Update version info in ``planemo/__init__.py`` (n+1.dev0) and create new entry in HISTORY.rst.
* ``git add HISTORY.rst planemo/__init__.py; git commit -m "Start work on new version"``
* ``git push origin``
* ``git push --tags origin``
* Update planemo homebrew recipe to new version.
