============
Installation
============

The recommended approach for installing Planemo is to use Homebrew_ or
linuxbrew_. To install planemo this way simply use the ``brew`` command as
follows.

::

    $ brew tap galaxyproject/tap
    $ brew install planemo

This will install the latest Planemo release. ``brew`` provide a very nice
upgrade path as new versions of Planemo are released.

::

    $ brew update
    $ brew upgrade planemo

To install or upgrade to the latest development branch of Planemo simply add
the argument ``--HEAD`` to either ``install`` or ``upgrade``.

For a more traditional Python installation simply setup a virtualenv
for ``planemo`` (this example creates a new one in ``.venv``) and then
install with ``pip``.

::

    $ virtualenv .venv; . .venv/bin/activate
    $ pip install planemo

.. _Homebrew: http://brew.sh/
.. _linuxbrew: https://github.com/Homebrew/linuxbrew
