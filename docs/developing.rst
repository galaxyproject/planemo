==================
Release Checklist
==================

This page describes the process of releasing new versions of Planemo.

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
* ``make release VERSION=<old_version> NEW_VERSION=<new_version>``
* ``make push-release``
* Update planemo homebrew recipe to new version.
