==================
Virtual Appliance
==================

You can also use planemo as part of a Galaxy tool development virtual
appliance pre-configured with Planemo, Galaxy_, Docker_, linuxbrew_, and a
Codebox_ web-based IDE.

The Galaxy instance that runs in these appliances has been optimized for tool
development - Galaxy will monitor your tool directory for changes and reload
the tools as they are modified, the server will directly log you into Galaxy
as an admin (no need to worry about user management or configuration), and
Galaxy is configured to use a `PostgreSQL
<http://www.postgresql.org/>`_ database backend and execute jobs via `SLURM
<https://computing.llnl.gov/linux/slurm/>`_ for robustness.

Launching the Appliance (Docker)
==================================

The `Docker`_ version of the planemo appliance can be launched using the
following command (which will pull the appliance down from `Docker Hub
<https://registry.hub.docker.com/u/planemo/box/>`_).

:: 

    docker run -d -p 8010:80 -v `pwd`:/opt/galaxy/tools -t planemo/box

This assumes your tools are in your current working directory (just replace
`\`pwd\`` with a path to your tools if this is not the case).

The above command will print a container ID which you can later use to kill
Docker

:: 

    docker kill <container_id>

This Docker environment will contain your tools and modifications made to them
will be made directly to your filesystem - so they are persistent. Data loaded
into the Galaxy instance (history data for instance) will be lost when the
Docker container is stopped. Check out the `docker-galaxy-stable
<https://github.com/bgruening/docker-galaxy-stable>`_ project for information
on running persistent Galaxy processes in Docker.

Launching the Appliance (Vagrant)
==================================

The latest version`Vagrant`_ version of the planemo appliance can be found
`here <https://images.galaxyproject.org/planemo/latest.box>`_. Once you have
installed `Vagrant`_ (`download now <http://www.vagrantup.com/downloads>`_),
the appliance can be enabled by first creating a `Vagrantfile` in your tool
directory - the following demonstrates an example of such file.

.. literalinclude:: Vagrantfile

This file must literally be named ``Vagrantfile``. Next you will need to
startup the appliance. This is as easy as

::

	vagrant up


Using the Appliance
====================

If you have started the appliance using one of the above commands, your Galaxy
instance should be running at http://localhost:8010/ and your tools should
appear in the left panel.

A Codebox_ IDE is available at http://localhost:8010/ide/ which should open
right to your tools folder and which lets you open a real terminal. This
terminal lets you run ``planemo``, build ``Dockerfile`` s, manage Galaxy,
etc... right from the web browser. For wider monitors -
http://localhost:8010/planemo/ will display the Codebox_ IDE and Galaxy side
by side.

Building the Appliance
======================

These applicances are built using the `planemo-machine
<https://github.com/jmchilton/planemo-machine>`_ project which can be used to
build customized recipes of this nature or even appliance for cloud
environments such as Amazon Web Services and Google Compute Engine.

.. _Galaxy: (http://galaxyproject.org/)
.. _Docker: https://www.docker.com/
.. _linuxbrew: https://github.com/Homebrew/linuxbrew
.. _Vagrant: https://www.vagrantup.com/
.. _Codebox: https://www.codebox.io/
