.. _shed:
=============================
Publishing to the Tool Shed
=============================

The `Galaxy Tool Shed`_ (referred to colloquially in Planemo as the "shed")
can store Galaxy tools, dependency definitions, and workflows among other
Galaxy artifacts. This guide will assume some basic familiarity with the shed
- please review the `Tool Shed Wiki`_ for and introduction.

Configuring a Shed Account
=============================

Before getting started, it is a good idea to have accounts on the Galaxy `test 
<https://testtoolshed.g2.bx.psu.edu/>`__ and (optionally) `main
<https://toolshed.g2.bx.psu.edu/>`__ Tool Sheds. Also, if you haven't initialized a
global Planemo configuration file (``~/.planemo.yml``) this can be done with.

::

    planemo config_init

This will populate a template ``~/.planemo.yml`` file and provide locations to
fill in shed credentials for the test and main Tool Sheds. For each shed, fill
in  either an API ``key`` or an ``email`` and ``password``. Also specify the
``shed_username`` created when registering shed accounts. All these options
can be specified and/or overridden on each planemo command invocation - but
that becomes tedious quickly.

Creating a Repository
=============================

Planemo can be used to used to publish "repositories" to the Tool Shed. A
single GitHub repository or locally managed directory of tools may correspond
to any number of Tool Shed repositories. Planemo maps files to Tool Shed
repositories using a special file called ``.shed.yml``.

From a directory containing tools, a `package definition`_. etc... the ``shed_init``
`command <http://planemo.readthedocs.org/en/latest/commands.html#shed-init-command>`__
can be used to bootstrap a new ``.shed.yml`` file.

::

    planemo shed_init --name=<name>
                      --owner=<shed_username>
                      --description=<short description>
                      [--remote_repository_url=<URL to .shed.yml on github>]
                      [--homepage_url=<Homepage for tool.>]
                      [--long_description=<long description>]
                      [--category=<category name>]*

There is not a lot of magic happening here, this file could easily be created
directly with a text editor - but the command has a ``--help`` to assist you
and does some very basic validation.

.. note:: Periods and hyphens are disallowed in repository names, it is
          recommended replacing periods in the version number with underscores.

          The following naming conventions are recommended and in some cases
          Planemo will determine the repository type based on adherence to these
          conventions (for packages and suites specifically).
 
          +-----------------------+-----------------------------+-----------------------------+
          | Repository Type       | Recommended Name            | Examples                    |
          +=======================+=============================+=============================+
          | Data Managers         | ``data_manager_$name``      | ``data_manager_bowtie2``    |
          +-----------------------+-----------------------------+-----------------------------+
          | Packages              | ``package_$name_$version``  | ``package_aragorn_1_2_36``  |
          +-----------------------+-----------------------------+-----------------------------+
          | Tool Suites           | ``suite_$name``             | ``suite_samtools``          |
          +-----------------------+-----------------------------+-----------------------------+
          | Tools                 | ``$name``                   | ``stringtie``, ``bedtools`` |
          +-----------------------+-----------------------------+-----------------------------+

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

    planemo shed_create --shed_target testtoolshed


Updating a Repository
=============================

Ensure the Galaxy Test Tool Shed is enabled in Galaxy's
``config/tool_sheds_conf.xml`` file and install and test the new repository.

If modifications are required these can be reviewed using the ``shed_diff``
command.

::

    planemo shed_diff --shed_target testtoolshed

.. note:: If you look at `tools-iuc`_ you will see it is common practice to leave
          details such as shed target and changeset_revision from
          ``tool_dependencies.xml`` and ``repository_dependencies.xml`` files. These 
          are required by the Tool Shed but it will populate them on upload and 
          leaving them blank allows uploading the same artifacts to the Test and
          Main sheds. The upshot of this is however is that ``shed_diff`` will always 
          print diffs on these artifacts.

Modified artifacts can be uploaded using the following command.

::

    planemo shed_update --check_diff --shed_target testtoolshed

The ``--check_diff`` option here will ensure there are significant differnces
before uploading new contents to the tool shed.

Once your artifacts are ready for publication to the main Tool Shed, the
following commands to create a repository there and populate it with your
repository contents.

::

    planemo shed_create

Advanced Usage
=============================

The above usage is relatively straight forward - it will map the current
directory to a single repository in the Tool Shed.

See `Pull Request 143`_ and linked examples for details on more advanced
options such as mapping each tool to its own repository automatically (`a best
practice <https://wiki.galaxyproject.org/AToolOrASuitePerRepository>`__) or
building enitrely custom repository definitions manually.

.. _Galaxy Tool Shed: https://toolshed.g2.bx.psu.edu/
.. _Tool Shed Wiki: https://wiki.galaxyproject.org/ToolShed
.. _package definition: https://wiki.galaxyproject.org/PackageRecipes
.. _`tools-devteam`: https://github.com/galaxyproject/tools-devteam
.. _`tools-iuc`: https://github.com/galaxyproject/tools-iuc
.. _Pull Request 143: https://github.com/galaxyproject/planemo/pull/143
