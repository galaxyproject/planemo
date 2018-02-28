.. figure:: https://raw.githubusercontent.com/jmchilton/planemo/master/docs/planemo_logo.png
   :alt: Planemo Logo
   :align: center
   :figwidth: 100%
   :target: https://github.com/galaxyproject/planemo

Command-line utilities to assist in building and publishing Galaxy_ tools.

.. image:: https://readthedocs.org/projects/planemo/badge/?version=latest
   :target: http://planemo.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://badge.fury.io/py/planemo.svg
   :target: https://pypi.python.org/pypi/planemo/
   :alt: Planemo on the Python Package Index (PyPI)

.. image:: https://travis-ci.org/galaxyproject/planemo.png?branch=master
   :target: https://travis-ci.org/galaxyproject/planemo
   :alt: Build Status

.. image:: https://coveralls.io/repos/galaxyproject/planemo/badge.svg?branch=master
   :target: https://coveralls.io/r/galaxyproject/planemo?branch=master
   :alt: Coverage Status

* Free software: Academic Free License version 3.0
* Documentation: https://planemo.readthedocs.io.
* Code: https://github.com/galaxyproject/planemo


Quick Start
-----------

This quick start demonstrates using ``planemo`` commands to help
develop Galaxy tools.

-----------------
Obtaining
-----------------

For a traditional Python installation of Planemo, first set up a virtualenv
for ``planemo`` (this example creates a new one in ``.venv``) and then
install with ``pip``. Planemo requires pip 7.0 or newer.

::

    $ virtualenv .venv; . .venv/bin/activate
    $ pip install "pip>=7" # Upgrade pip if needed.
    $ pip install planemo

For information on updating Planemo, installing the latest development release,
or installing planemo via bioconda - checkout the `installation
<http://planemo.readthedocs.io/en/latest/installation.html>`__ 
documentation.

Planemo is also available as a `virtual appliance
<https://planemo.readthedocs.io/en/latest/appliance.html>`_ bundled
with a preconfigured Galaxy server and set up for Galaxy tool development.
You can choose from open virtualization format (OVA_, .ova), Docker_,
or Vagrant_ appliances.

--------------
Basics
--------------

This quick start will assume you have a directory with one or more Galaxy
tool XML files. If no such directory is available, one can be quickly created for
demonstrating ``planemo`` as follows ``mkdir mytools; cd mytools; planemo
project_init --template=demo``.

Planemo can check tools containing XML for common problems and best
practices using the ``lint`` `command <http://planemo.readthedocs.org/en/latest/commands.html#lint-command>`_
(also aliased as ``l``).

::

    $ planemo lint

Like many ``planemo`` commands - by default this will search the
current directory and use all tool files it finds. It can be explicitly
passed a path to tool files or a directory of tool files.

::

    $ planemo l randomlines.xml

The ``lint`` command takes in additional options related to
reporting levels, exit code, etc. These options are described
in the `docs <http://planemo.readthedocs.org/en/latest/commands.html#lint-command>`_
or (like with all commands) can be accessed by passing ``--help`` to it.

::

    $ planemo l --help
    Usage: planemo lint [OPTIONS] TOOL_PATH

Once tools are syntactically correct - it is time to test. The ``test``
`command <http://planemo.readthedocs.org/en/latest/commands.html#test-command>`__
can be used to test a tool or a directory of tools.

::

	$ planemo test --galaxy_root=../galaxy randomlines.xml

If no ``--galaxy_root`` is defined, ``planemo`` will check for a default in
`~/.planemo.yml
<http://planemo.readthedocs.org/en/latest/configuration.html>`_ and finally
search the tool's parent directories for a Galaxy root directory.
Planemo can also download and configure a disposable Galaxy instance for
testing. Pass ``--install_galaxy`` instead of ``--galaxy_root``.

::

	$ planemo t --install_galaxy

Planemo will create a HTML output report in the current directory named
``tool_test_output.html`` (override with ``--test_output``). See an
`example <http://galaxyproject.github.io/planemo/tool_test_viewer.html?test_data_url=https://gist.githubusercontent.com/jmchilton/9d4351c9545d34209904/raw/9ed285d3cf98e435fc4a743320363275949ad63c/index>`_
of such a report for Tophat.

Once tools have been linted and tested - the tools can be viewed in a
Galaxy interface using the ``serve`` (``s``) `command
<http://planemo.readthedocs.org/en/latest/commands.html#serve-command>`__.

::

	$ planemo serve

Like ``test``, ``serve`` requires a Galaxy root and one can be
explicitly specified with ``--galaxy_root`` or installed dynamically
with ``--install_galaxy``.

For more information on building Galaxy tools in general please check out
`Building Galaxy Tools Using Planemo`_.

---------
Tool Shed
---------

Planemo can help you publish tools to the Galaxy Tool Shed.
Check out `Publishing to the Tool Shed`_ for more information.

------
Conda
------

Planemo can help develop tools and Conda packages in unison.
Check out `Dependencies and Conda`_ for more information.

--------
Docker
--------

Planemo can help develop tools deployable via Docker.
Check out `Dependencies and Docker`_ for more information.

--------------------------
Common Workflow Language
--------------------------

Planemo includes experimental support for running a subset of valid
`Common Workflow Language`_ (CWL) tools using either the reference implementation cwltool or
a fork of Galaxy enhanced to run CWL tools.

::

    $ planemo project_init --template cwl_draft3_spec
    $ planemo serve --cwl cwl_draft3_spec/cat1-tool.cwl

Checko out `Building Common Workflow Language Tools`_ for more information.

.. _Galaxy: http://galaxyproject.org/
.. _GitHub: https://github.com/
.. _Conda: http://conda.pydata.org/
.. _Docker: https://www.docker.com/
.. _Vagrant: https://www.vagrantup.com/
.. _Travis CI: http://travis-ci.org/
.. _`tools-devteam`: https://github.com/galaxyproject/tools-devteam
.. _`tools-iuc`: https://github.com/galaxyproject/tools-iuc
.. _Building Galaxy Tools Using Planemo: http://planemo.readthedocs.io/en/latest/writing_standalone.html
.. _Publishing to the Tool Shed: http://planemo.readthedocs.org/en/latest/publishing.html
.. _Dependencies and Conda: http://planemo.readthedocs.io/en/latest/writing_advanced.html#dependencies-and-conda
.. _Dependencies and Docker: http://planemo.readthedocs.io/en/latest/writing_advanced.html#dependencies-and-docker
.. _Common Workflow Language: http://common-workflow-language.github.io
.. _Building Common Workflow Language Tools: http://planemo.readthedocs.io/en/latest/writing_cwl_standalone.html
.. _OVA: https://en.wikipedia.org/wiki/Open_Virtualization_Format
