============
Installation
============

pip_
============

For a traditional Python installation of Planemo, first set up a virtual environment
for ``planemo`` (this example creates a new one in ``planemo``) and then
install with ``pip``. Planemo requires pip 7.0 or newer.

::

    $ python3 -m venv planemo
    $ . planemo/bin/activate
    $ pip install planemo

When installed this way, planemo can be upgraded as follows:

::

    $ . planemo/bin/activate
    $ pip install -U planemo

To install or update to the latest development branch of planemo with ``pip``, 
use the  following ``pip install`` idiom instead:

::

    $ pip install -U git+git://github.com/galaxyproject/planemo.git

If your ``PATH`` contains a Python installed through Conda it should likely not be used to run Planemo,
consider using the ``virtualenv`` argument ``-p`` to point at a non-Conda Python 2 executable installed
natively on your system or using a tool such pyenv_. ``virtualenv`` can be installed via Conda, pyenv_,
or a package manager - it should make no difference.

Planemo runs on Python 3.6 or newer. Planemo can be used to run multiple versions of Galaxy,
but please note that the last Galaxy release that fully supports Python 2.7 is 19.09.

Conda_ (Experimental)
=======================

Another approach for installing Planemo is to use Conda_
(most easily obtained via the
`Miniconda Python distribution <http://conda.pydata.org/miniconda.html>`__).
Afterwards run the following commands.

::

    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda install planemo

Galaxy is known to have issues when running with a Conda Python so this approach
should be considered experimental for now. If you have problems with it or hacks to
make it work better - please report them.

.. _pip: https://pip.pypa.io/
.. _Conda: http://conda.pydata.org/docs/
.. _pyenv: https://github.com/pyenv/pyenv
