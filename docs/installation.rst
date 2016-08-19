============
Installation
============

pip_
============

For a more traditional Python installation of Planemo set up a virtualenv
for ``planemo`` (this example creates a new one in ``.venv``) and then
install with ``pip``. Planemo requires pip 7.0 or newer.

::

    $ virtualenv .venv; . .venv/bin/activate
    $ pip install "pip>=7" # Upgrade pip if needed.
    $ pip install planemo

When installed this way, planemo can be upgraded as follows:

::

    $ . .venv/bin/activate
    $ pip install -U planemo

To install or update to the latest development branch of planemo with ``pip``, 
use the  following ``pip install`` idiom instead:

::

    $ pip install -U git+git://github.com/galaxyproject/planemo.git


Homebrew_
============

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

Conda_
============

Another approach for installing Planemo is to use Conda_
(most easily obtained via the
`Miniconda Python distribution <http://conda.pydata.org/miniconda.html>`__).
Afterwards run the following commands.

::

    $ conda config --add channels r
    $ conda config --add channels bioconda
    $ conda install planemo

.. _pip: https://pip.pypa.io/
.. _Homebrew: http://brew.sh/
.. _linuxbrew: https://github.com/Homebrew/linuxbrew
.. _Conda: http://conda.pydata.org/docs/
