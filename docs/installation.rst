============
Installation
============

The recommended approach for installing Planemo is to use Conda_
(most easily obtained via the
`Miniconda Python distribution <http://conda.pydata.org/miniconda.html>`__).
Afterwards run the following commands.

::
    $ conda config --add channels r
    $ conda config --add channels bioconda
    $ conda install planemo

Another approach for installing Planemo is to use Homebrew_ or
linuxbrew_. To install Planemo this way use the ``brew`` command as
follows.

::
    $ brew tap galaxyproject/tap
    $ brew install planemo

This will install the latest Planemo release. ``brew`` provide a very nice
upgrade path as new versions of Planemo are released.

::
    $ brew update
    $ brew upgrade planemo

To install or upgrade to the latest development branch of Planemo add
the argument ``--HEAD`` to either ``install`` or ``upgrade``.

For a more traditional Python installation set up a virtualenv
for ``planemo`` (this example creates a new one in ``.venv``) and then
install with ``pip``.

::
    $ virtualenv .venv; . .venv/bin/activate
    $ pip install planemo

.. _Homebrew: http://brew.sh/
.. _linuxbrew: https://github.com/Homebrew/linuxbrew
.. _Conda: http://conda.pydata.org/docs/
