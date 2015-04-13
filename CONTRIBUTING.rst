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

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

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

    $ git clone git@github.com:your_name_here/planemo.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv planemo
    $ cd planemo/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass ``flake8``
and the tests::

    $ flake8 planemo tests
    $ python setup.py test

To get ``flake8``, just ``pip install`` them into your virtualenv.

.. including testing other Python versions with tox
..    $ python setup.py test
..    $ tox
..
..   To get flake8 and tox, just pip install them into your virtualenv.

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
2. The pull request should work for Python 2.6 and 2.7. Check
   https://travis-ci.org/galaxyproject/planemo/pull_requests
   and make sure that the tests pass for all supported Python versions.

.. Tips
.. ----
.. To run a subset of tests::
..    $ python -m unittest tests.test_planemo
