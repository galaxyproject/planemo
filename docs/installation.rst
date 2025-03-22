============
Installation
============

pip_
====

For a traditional Python installation of Planemo, first set up a virtual environment
for ``planemo`` (this example creates a new one in ``planemo``) and then
install with ``pip``. Planemo requires pip 7.0 or newer.

::

    $ python -m venv planemo
    $ . planemo/bin/activate
    $ pip install planemo

When installed this way, planemo can be upgraded as follows:

::

    $ . planemo/bin/activate
    $ pip install -U planemo

To install or update to the latest development branch of planemo with ``pip``, 
use the following ``pip install`` idiom instead:

::

    $ pip install -U git+git://github.com/galaxyproject/planemo.git

If your ``PATH`` contains a Python installed through Conda it should likely not be used to run Planemo,
please consider using a tool such uv_ to manage native Python installations and environments.

Planemo runs on Python 3.8 or newer. Planemo can be used to run multiple versions of Galaxy,
but please note that the last Galaxy release that fully supports Python 2.7 is 19.09.

Conda_ (Experimental)
=====================

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

uv_
===

`uv <https://github.com/astral-sh/uv>`__ is a modern Python package manager that can be used to install Planemo.
Information on installing uv_ can be found at https://docs.astral.sh/uv/getting-started/installation/.
Once uv_ has been installed and its environment sourced in your shell, you can install Planemo using the following commands:

::

    $ uv tool install planemo 

Installing Planemo as a tool will setup a global shim that manages the Python environment used for Planemo
without having to manage this yourself like the ``pip`` method described above. The Planemo  executable 
``planemo`` will just be on your path whenever the uv environment is active.

uv can also be used to upgrade Planemo to the latest version with the following command:

::

    $ uv tool upgrade planemo

rye_
====

`rye <https://github.com/astral-sh/rye>`__ is another modern Python package manager.
Information on installing rye_ can be found at https://rye.astral.sh/guide/installation/.
Once rye_ has been installed and its environment sourced in your shell, you can install Planemo using the following commands:

::

    $ rye tools install planemo

Installing Planemo as a tool will setup a global shim that manages the Python environment used for Planemo
without having to manage this yourself like the ``pip`` method described above. The Planemo executable 
``planemo`` will just be on your path whenever the rye environment is active.

.. note::
    The options ``--docker`` and ``--biocontainers`` for ``planemo serve`` and ``planemo test`` require that 
    docker is installed on the system and that it can be run without root privileges. 
    See `Docker installation <https://docs.docker.com/engine/install>`__ and 
    `run Docker without root privileges <https://docs.docker.com/engine/install/linux-postinstall>`__ for further instructions.

.. _pip: https://pip.pypa.io/
.. _Conda: http://conda.pydata.org/docs/
.. _uv: https://github.com/astral-sh/uv
.. _rye: https://github.com/astral-sh/rye
