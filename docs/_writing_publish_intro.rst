.. _shed:

Publishing to the Tool Shed
====================================================

Now that the tool is working and useful - it is time to publish it to the tool
shed. The `Galaxy Tool Shed`_ (referred to colloquially in Planemo as the
"shed") can store Galaxy tools, dependency definitions, and workflows among
other Galaxy artifacts.

-------------------------------------------------
Configuring a Shed Account
-------------------------------------------------

The `planemo <http://planemo.readthedocs.org/en/latest/appliance.html>`__
appliance comes pre-configured with a local Tool Shed and planemo is
configured to talk to it via ``~/.planemo.yml``. Check out the `publishing docs
<http://planemo.readthedocs.org/en/latest/publishing.html>`__ for information
on setting up this ``~/.planemo.yml`` file on your development environment.
 
-------------------------------------------------
Creating a Repository
-------------------------------------------------

Planemo can be used to used to publish "repositories" to the Tool Shed. A
single GitHub repository or locally managed directory of tools may correspond
to any number of Tool Shed repositories. Planemo maps files to Tool Shed
repositories using a special file called ``.shed.yml``.

From a directory containing tools the ``shed_init``
`command <http://planemo.readthedocs.org/en/latest/commands.html#shed-init-command>`__
can be used to bootstrap a new ``.shed.yml`` file.

::

    planemo shed_init --name=seqtk_seq
                      --owner=planemo
                      --description=seqtk_seq
                      --long_description="Tool that converts FASTQ to FASTA files using seqtk"
                      --category=Fastq Manipulation

There is not a lot of magic happening here, this file could easily be created
directly with a text editor - but the command has a ``--help`` to assist you
and does some very basic validation.

More information on ``.shed.yml`` can be found as part of the IUC `best
practice documentation
<http://galaxy-iuc-standards.readthedocs.org/en/latest/best_practices/shed_yml.html>`__.

After reviewing ``.shed.yml``, this configuration file and relevant shed
artifacts can be quickly linted using the following command.

::

    planemo shed_lint --tools

Once the details the ``.shed.yml`` are set and it is time to create the remote
repository and upload artifacts to it - the following two commands can be used
- the first only needs to be run once and creates the repository based the
metadata in ``.shed.yml`` and the second uploads your actual artifacts to it.

::

    planemo shed_create --shed_target local

Now navigate to the local tool shed (likely at `http://localhost:9009/
<http://localhost:9009/>`__). You can login with login ``planemo@test.com``
and password ``planemo``.

-------------------------------------------------
Updating a Repository
-------------------------------------------------

::

    planemo shed_update --check_diff --shed_target local

Once tools and reqiured dependency files have been published to the tool shed,
the actual shed dependencies can be automatically and installed and tool 
tests ran using the command::

    planemo shed_serve --shed_target local

Once your artifacts are ready for publication to the main Tool Shed, the
following commands to create a repository there and populate it with your
repository contents.

::

    planemo shed_create

The planemo machine isn't preconfigured to allow publishing to the main tool
shed so this command will not work. See the more complete `publishing docs
<http://planemo.readthedocs.org/en/latest/publishing.html>`__ for full details
about how to setup Planemo to publish to the main and test tool shed - the
process is very similar.

.. _Galaxy Tool Shed: https://toolshed.g2.bx.psu.edu/
.. _Tool Shed Wiki: https://wiki.galaxyproject.org/ToolShed
.. _package definition: https://wiki.galaxyproject.org/PackageRecipes
.. _`tools-iuc`: https://github.com/galaxyproject/tools-iuc
