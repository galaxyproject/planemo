.. figure:: https://raw.githubusercontent.com/jmchilton/planemo/master/docs/planemo_logo.png
   :alt: Planemo Logo
   :align: center
   :figwidth: 100%
   :target: https://github.com/galaxyproject/planemo

Command-line utilities to assist in developing Galaxy_ and `Common Workflow Language`_ artifacts -
including tools, workflows, and training materials.

.. image:: https://readthedocs.org/projects/planemo/badge/?version=latest
   :target: http://planemo.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://badge.fury.io/py/planemo.svg
   :target: https://pypi.python.org/pypi/planemo/
   :alt: Planemo on the Python Package Index (PyPI)

.. image:: https://github.com/galaxyproject/planemo/workflows/Python%20CI/badge.svg
   :target: https://github.com/galaxyproject/planemo/actions?query=workflow%3A%22Python+CI%22

.. image:: https://travis-ci.org/galaxyproject/planemo.png?branch=master
   :target: https://travis-ci.org/galaxyproject/planemo
   :alt: Build Status

* Free software: Academic Free License version 3.0
* Documentation: https://planemo.readthedocs.io.
* Code: https://github.com/galaxyproject/planemo

Quick Start
-----------

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
or installing Planemo via `Bioconda <https://github.com/bioconda/bioconda-recipes>`__
- checkout the `installation <http://planemo.readthedocs.io/en/latest/installation.html>`__
documentation.

Planemo is also available as a `virtual appliance
<https://planemo.readthedocs.io/en/latest/appliance.html>`_ bundled
with a preconfigured Galaxy server and set up for Galaxy_ and
`Common Workflow Language`_ tool development.
You can choose from open virtualization format (OVA_, .ova) or Docker_
appliances.

-----------------
Basics - Galaxy
-----------------

This quick start will assume you have a directory with one or more Galaxy
tool XML files. If no such directory is available, one can be quickly created for
demonstrating ``planemo`` as follows ``project_init --template=demo mytools; cd mytools``.

Planemo can check tool XML files for common problems and best
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

If no ``--galaxy_root`` is defined, Planemo will download and configure
a disposable Galaxy instance for testing.

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

For more information on developing Galaxy workflows with Planemo checkout
`best practices for Galaxy Workflows`_ and the description of Planemo's
`test format`_. For information on developing Galaxy training materials
checkout the `contributing documentation <https://training.galaxyproject.org/training-material/topics/contributing/>`__
on training.galaxyproject.org.


------------
ToolFactory
------------

Planemo can serve the `Galaxy ToolFactory tool <https://github.com/fubar2/toolfactory>`__ on `http://localhost:9090`

::

  $ planemo tool_factory --galaxy_root /galaxy_central --extra_tools ./mytools

The ToolFactory works like other Galaxy tools through the Galaxy GUI to fill in a form, but it can generate a new Galaxy tool and test
it with Planemo. Sample data file inputs should be uploaded to the working history before starting the ToolFactory. These are needed to create a
test as part of the new generated tool archive. Galaxy as an integrated development environment is possible but clumsy for anything other than relatively simple tools.
Conda packages can be used or a script can be provided for a scripting language interpreter. The command line for the package or script is made up from
the input files needed and the output files generated, plus any number of other command line parameters with argparse or positional style
formatting. The number of files and parameters is limited only by your patience with the GUI. A toolshed ready archive is produced in the history containing a test
based on the sample inputs and parameter settings provided on the form. Note that tools are normal Galaxy tools and are provided as archives ready to upload to
a toolshed.

The optional `--galaxy_root` pointing to a Galaxy source code directory will save time downloading a new copy
each time it is run. The optional `--extra_tools` can be a path containing any tools you have previously built if you
download the [toolname].toolshed.tar.gz file and unarchive it into that extra tools path. These tools will appear
in the Galaxy served by planemo along with the ToolFactory tool itself.

For dedicated ToolFactory access, the `ToolFactory Docker <https://github.com/fubar2/toolfactory-galaxy-docker>`__ container offers persistence
and has an inbuilt toolshed so new tools can be installed back to the Galaxy for additional testing once they pass the inbuilt test.
It is based on `the stable Galaxy docker <https://github.com/bgruening/docker-galaxy-stable>`__ and also uses Planemo under the
hood.

----------------------------------
Basics - Common Workflow Language
----------------------------------

This quick start will assume you have a directory with one or more `Common Workflow
Language`_ YAML files. If no such directory is available, one can be quickly created for
demonstrating ``planemo`` as follows ``planemo project_init --template=seqtk_complete_cwl mytools; cd mytools``.

Planemo can check tools YAML files for common problems and best
practices using the ``lint`` `command <http://planemo.readthedocs.org/en/latest/commands.html#lint-command>`_
(also aliased as ``l``).

::

    $ planemo lint

Like many ``planemo`` commands - by default this will search the
current directory and use all tool files it finds. It can be explicitly
passed a path to tool files or a directory of tool files.

::

    $ planemo l seqtk_seq.cwl

The ``lint`` command takes in additional options related to
reporting levels, exit code, etc. These options are described
in the `docs <http://planemo.readthedocs.org/en/latest/commands.html#lint-command>`_
or (like with all commands) can be accessed by passing ``--help`` to it.

::

    $ planemo l --help
    Usage: planemo lint [OPTIONS] TOOL_PATH

Once tools are syntactically correct - it is time to test. The ``test``
`command <http://planemo.readthedocs.org/en/latest/commands.html#test-command>`__
can be used to test a CWL tool, workflow, or a directories thereof.

::

  $ planemo test --engine cwltool seqtk_seq.cwl

Planemo will create a HTML output report in the current directory named
``tool_test_output.html``. Check out the file ``seqtk_seq_tests.yml`` for
an example of Planemo test for a CWL tool. A test consists of any number of
jobs (with input descriptions) and corresponding output assertions.

Checkout the `Commmon Workflow User Guide`_ for more information on developing
CWL tools in general and  `Building Common Workflow Language Tools`_ for more
information on using Planemo to develop CWL tools.

---------
Tool Shed
---------

Planemo can help you publish tools to the Galaxy Tool Shed.
Check out `Publishing to the Tool Shed`_ for more information.

------
Conda
------

Planemo can help develop tools and Conda packages in unison.
Check out the `Galaxy <http://planemo.readthedocs.io/en/latest/writing_advanced.html#dependencies-and-conda>`__ or `CWL
<http://planemo.readthedocs.io/en/latest/writing_advanced_cwl.html#dependencies-and-conda-cwl>`__ version of the "Dependencies and Conda" tutorial
for more information.

-----------------------
Docker and Containers
-----------------------

Planemo can help develop tools that run in "Best Practice" containers for
scientific workflows. Check out the `Galaxy <http://planemo.readthedocs.io/en/latest/writing_advanced.html#dependencies-and-containers>`__ or `CWL
<http://planemo.readthedocs.io/en/latest/writing_advanced_cwl.html#dependencies-and-containers-cwl>`__ version of the "Dependencies and Containers" tutorial for more information.

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
.. _Common Workflow Language: https://www.commonwl.org/
.. _Commmon Workflow User Guide: http://www.commonwl.org/user_guide/
.. _Building Common Workflow Language Tools: http://planemo.readthedocs.io/en/latest/writing_cwl_standalone.html
.. _OVA: https://en.wikipedia.org/wiki/Open_Virtualization_Format
.. _test format: https://planemo.readthedocs.io/en/latest/test_format.html
.. _best practices for Galaxy Workflows: https://planemo.readthedocs.io/en/latest/best_practices_workflows.html
