Running workflows on a server or Slurm cluster
==============================================

Some data cannot be uploaded to a public Galaxy instance because of legal,
policy, or storage constraints. Planemo can launch a temporary Galaxy instance
on a server or cluster login node and run workflow jobs either locally or
through a workload manager such as Slurm. This avoids maintaining a permanent
Galaxy service while keeping data inside the local environment.

This guide uses the ``tutorial.ga`` workflow and ``tutorial-job.yml`` job file
introduced in `The Basics`_ and covers:

#. running the workflow on one server;
#. generating a Slurm job configuration;
#. submitting jobs through Slurm with or without TPV; and
#. choosing a database for concurrent workflows.

Requirements
------------

The server or cluster needs:

* Planemo, installed as described in :doc:`readme`;
* a Slurm client configured for the cluster;
* a DRMAA implementation for Slurm;
* Conda, Docker, or Singularity/Apptainer for tool dependencies; and
* a filesystem visible to both the login node and compute nodes.

PostgreSQL is optional but recommended for workflows that run many jobs
concurrently.

Slurm and DRMAA
~~~~~~~~~~~~~~~

Galaxy's Slurm runner uses the Distributed Resource Management Application API
(DRMAA). Ask the cluster administrators whether a DRMAA library is already
available. If it is not, `slurm-drmaa`_ is the recommended implementation.

For example, build slurm-drmaa 1.1.4 in a user-writable directory with:

.. code-block:: console

    $ wget https://github.com/natefoo/slurm-drmaa/releases/download/1.1.4/slurm-drmaa-1.1.4.tar.gz
    $ tar -xzf slurm-drmaa-1.1.4.tar.gz
    $ cd slurm-drmaa-1.1.4
    $ ./configure --prefix="$PWD/dist"
    $ make
    $ make install

Note the absolute path to the installed ``libdrmaa`` library. You will add it
to the generated job configuration below.

.. _slurm-drmaa: https://github.com/natefoo/slurm-drmaa

Tool dependencies
~~~~~~~~~~~~~~~~~

Galaxy can resolve tool dependencies with Conda or run tools in Docker or
Singularity/Apptainer containers. Containers are generally more reproducible,
but the available runtime and mount configuration depend on the cluster. The
examples below use Singularity because it is commonly available on HPC systems.

Database
~~~~~~~~

Planemo uses SQLite by default. SQLite is sufficient for small runs, but a
workflow with many simultaneous jobs can encounter database-locking errors.
Use an existing PostgreSQL service when one is available. Alternatively,
Planemo can launch a temporary PostgreSQL instance in a Singularity/Apptainer
container with ``--database_type postgres_singularity``.

Run on one server
-----------------

Run the example workflow locally with:

.. code-block:: console

    $ planemo run tutorial.ga tutorial-job.yml \
        --download_outputs \
        --output_directory . \
        --output_json output.json

Planemo launches a temporary Galaxy instance, runs the workflow, downloads its
outputs, and records their paths in ``output.json``.

Run through Slurm
-----------------

Generate the job configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use :doc:`commands/job_config_init` instead of maintaining a complete Galaxy
job configuration by hand:

.. code-block:: console

    $ planemo job_config_init \
        --runner slurm \
        --singularity \
        --galaxy_version 25.0

The command writes ``job_conf.yml`` in the current directory. Review that file
and customize the generated Slurm runner and environment:

* uncomment ``drmaa_library_path`` and set it to the absolute path of the
  cluster's ``libdrmaa`` library;
* set ``native_specification`` to the partition, time, memory, and other Slurm
  options required by the cluster; and
* adjust the generated Singularity command, volumes, or environment variables
  when the cluster requires them.

The generator sets ``tmp_dir: true`` by default. This asks Galaxy to manage a
temporary directory for each job and normally removes the need to set a shared
``TMPDIR`` manually. The Galaxy job working directories and input data still
need to reside on storage visible to the compute nodes.

Run the workflow
~~~~~~~~~~~~~~~~

Pass the generated file to ``planemo run``:

.. code-block:: console

    $ planemo run tutorial.ga tutorial-job.yml \
        --download_outputs \
        --output_directory . \
        --output_json output.json \
        --job_config_file job_conf.yml

For a temporary PostgreSQL database, add:

.. code-block:: console

    $ planemo run tutorial.ga tutorial-job.yml \
        --download_outputs \
        --output_directory . \
        --output_json output.json \
        --job_config_file job_conf.yml \
        --database_type postgres_singularity

This database mode requires a working Singularity/Apptainer installation and
storage suitable for the database container.

Slurm with TPV
--------------

`Total Perspective Vortex (TPV)`_ maps individual tools to job destinations
and resource requirements. Planemo can generate the surrounding Galaxy and TPV
configuration while retaining the Slurm and container settings from the
previous example:

.. code-block:: console

    $ planemo job_config_init \
        --runner slurm \
        --tpv \
        --singularity \
        --galaxy_version 25.0

Edit the generated ``job_conf.yml`` to configure ``drmaa_library_path`` and
the ``tpvdb_slurm`` destination for the local cluster. In particular, set its
``native_specification`` template to translate TPV's ``{cores}``, ``{mem}``,
and other resource values into the Slurm options used by the site. The
generated configuration includes the shared Galaxy TPV database at
``https://gxy.io/tpv/db.yml``; local configuration can override entries from
that database.

Even when the scheduler does not require ``--ntasks``, make sure the native
specification communicates the requested core count. Otherwise a
multithreaded tool may be scheduled with only one core.

For configuration details and examples, see the `TPV documentation`_ and the
Galaxy Training Network's `job destination tutorial`_.

.. _Total Perspective Vortex (TPV): https://total-perspective-vortex.readthedocs.io/
.. _TPV documentation: https://total-perspective-vortex.readthedocs.io/en/latest/topics/tpv_by_example.html
.. _job destination tutorial: https://training.galaxyproject.org/training-material/topics/admin/tutorials/job-destinations/tutorial.html

Troubleshooting
---------------

Database is locked
~~~~~~~~~~~~~~~~~~

Errors containing ``sqlite3.OperationalError: database is locked`` usually
mean the workflow has more concurrent database activity than SQLite can
comfortably handle. Re-run with PostgreSQL, using either an existing database
or ``--database_type postgres_singularity``.

Galaxy exits during tool installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A later ``bioblend.ConnectionError`` can mean the temporary Galaxy process
exited while Planemo was waiting for a Tool Shed installation; it is not
necessarily a network failure in BioBlend. Inspect Planemo's Galaxy log for the
earlier server error. On a slow cluster, also check filesystem performance,
available memory, and the scheduler or container runtime logs.

Port already in use
~~~~~~~~~~~~~~~~~~~

If an earlier Planemo or Galaxy process is still using the configured port,
stop that process or select another port with ``--port``. Use the cluster's
normal process-monitoring tools to identify old ``planemo`` and Galaxy
processes before terminating them.
