============
Contributing
============

Please note that this project is released with a `Contributor Code of Conduct 
<https://planemo.readthedocs.org/en/latest/conduct.html>`__. By participating
in this project you agree to abide by its terms.

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/galaxyproject/planemo/issues.

If you are reporting a bug, please include:

* Your operating system name and version, versions of other relevant software 
  such as Galaxy or Docker.
* Links to relevant tools.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with
"enhancement" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Planemo could always use more documentation, whether as part of the
official Planemo docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/galaxyproject/planemo/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* This will hopefully become a community-driven project and contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `planemo` for local development.

1. Fork the `planemo` repo on GitHub.
2. Clone your fork locally::

    $ git clone --recursive git@github.com:your_name_here/planemo.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ make setup-venv

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass ``flake8``
   and the tests
   
   ::

       $ make lint
       $ make test
   
   If the modification doesn't affect code that configures and runs Galaxy - 
   skipping a couple tests that will cause Galaxy and its dependencies to be
   downloaded results in a significant speed up. This subset of tests can be
   run with ``make quick-test``.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring.
2. The pull request should work for Python >=3.6. Check
   https://travis-ci.org/galaxyproject/planemo/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

    % make tox ENV=py37-unit ARGS='--tests tests/test_shed_upload.py'

This will use Tox_ to run the specified tests using Python 3.7. ``ENV`` here
can be used to specify different Python version (e.g. ``py36`` or
``py37``).

Even more granularity is also possible by specifying specific test methods.::

    make tox ENV=py37-unit ARGS='--tests tests/test_shed_upload.py:ShedUploadTestCase.test_tar_from_git'


``tox`` can be used to run tests directly also (use ``. .venv/bin/activate``
to ensure it is on your ``PATH``).

::

    tox -e py37-unit -- --tests tests/test_shed_upload.py

Tox_ itself is configured to wrap nose_. One can skip Tox_ and run
``nosetests`` directly.

::

    pytest tests/test_shed_upload.py

Tox_
~~~~~~~~~~~

Tox_ is a tool to automate testing across different Python versions. The
``tox`` executable can be supplied with a ``-e`` argument to specify a
testing environment. Planemo defines the following environments:

``py37-lint``
    Lint the planemo code using Python 3.7.

``py37-lint_docs``
    Lint the docs reStructuredText.

``py37-lint_docstrings``
    Lint the project Python docstrings (doesn't pass currently).

``py36-unit-quick``
    Run the fastest unit tests (with least external dependencies) on Python 3.6.

``py36-unit-nonredundant-noclientbuild-gx-2005``
    Run tests that are marked as targeting a Galaxy branch and test against Galaxy 20.05.
    Skip tests that are marked as redundant or that require a Galaxy client build.

``py37-unit-gx-dev``
    Run tests that are marked as targeting a Galaxy branch and test against Galaxy's dev branch.

Pre-commit Hooks
~~~~~~~~~~~~~~~~~~~~~

Planemo pull requests are automatically linted and tested using `TravisCI
<https://travis-ci.org/galaxyproject/planemo>`__. A git pre-commit `hook
<http://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>`__ can be setup
to lint and/or test Planemo before committing to catch problems that would
be detected by TravisCI as early as possible.

The following command will install a pre-commit hook that lints the Planemo
code::

    make setup-git-hook-lint

To also run the faster planemo tests, the following command can be used to 
setup a more rigorous pre-commit hook::

    make setup-git-hook-lint-and-test

.. _Tox: https://tox.readthedocs.org/en/latest/
.. _nose: https://nose.readthedocs.org/en/latest/
