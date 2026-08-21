Running workflow on a server or cluster
============================================

In some situations, users may be unable to upload data to public Galaxy instances due to legal or policy restrictions.
Even if they have access to local cluster systems, these environments often do not permit hosting web servers or lack
dedicated personnel to maintain a Galaxy instance.
Planemo addresses this challenge by enabling users to launch a temporary local Galaxy instance, allowing workflow
execution without the need for a permanent server setup.
Jobs can be executed directly on the server or submitted to a workload manager such as Slurm.

This guide explains how to run a workflow:
 1. On a standalone server
 2. On a Slurm cluster using DRMAA


Requirements
---------------------

- planemo: install instructions are available at :doc:`readme`
- DRMAA library
- Tools dependencies
- Database (optional)

DRMAA library
~~~~~~~~~~~~~~~~~~
To enable Galaxy to submit jobs to a cluster via Slurm, you need the DRMAA
(Distributed Resource Management Application API) library, which provides a standard interface 
for job submission. If your system does not already have the DRMAA library installed, you must
compile it manually. For Slurm, the recommended implementation is 
`slurm-drmaa <https://github.com/natefoo/slurm-drmaa>`__.

Follow these steps to download and build the DRMAA library for Slurm:

:: 
   
   # Download and extract source code
   wget https://github.com/natefoo/slurm-drmaa/releases/download/1.1.4/slurm-drmaa-1.1.4.tar.gz 
   tar -xvzf  slurm-drmaa-1.1.4.tar.gz
   cd slurm-drmaa-1.1.4 && mkdir dist
   
   # Configure build
   ./autogen.sh
   ./configure --prefix=${PWD}/dist
   
   # Compile and install
   make
   make install
   
   # print path to library
   LIBDRMAA=`ls ${PWD}/dist/lib/libdrmaa.so`
   echo "Library: ${LIBDRMAA}"


Tool dependencies
~~~~~~~~~~~~~~~~~~
To run a tool/workflow on a server or cluster, the required tools must
be available on the system. Galaxy can manage dependencies
using Conda, Docker or Singularity. While Conda is an option, Docker or 
Singularity is recommended for its greater reliability. One of the following
environments need to be available on the system:

- Conda
- Docker
- Singularity


Database
~~~~~~~~~~~~~~~~~~
Galaxy tracks jobs and their statuses using a database, with SQLite as the default option.
However, when running complex workflows, SQLite can encounter issues if multiple jobs attempt
to update the database simultaneously. This can lead to database locking and premature workflow
termination.

To prevent these issues, it's recommended to use a PostgreSQL database instead of SQLite. If a
PostgreSQL database isn't readily available, Planemo can launch a temporary PostgreSQL instance
using Singularity by adding the ``--database_type postgres_singularity`` to the run command.


Running on workflow
------------------
Planemo can be configured in many ways to suit different user needs. In this guide, we'll demonstrate
how to run a workflow on both a local laptop/server and a SLURM cluster, using the same example workflow
introduced in `the basic <http://127.0.0.1:8000/running.html#the-basics>`_

Laptop/server
~~~~~~~~~~~~~~~~~~
::

  planemo run tutorial.ga tutorial-job.yml \
    --download_outputs \
    --output_directory . \
    --output_json output.json \
    --biocontainers

This command runs the workflow locally and attempts to resolve tool dependencies using containers. If the
workflow completes successfully, the output files will be downloaded to the current directory.

Slurm
~~~~~~~~~~~~~~~~~~

To run the workflow on a SLURM cluster, you’ll need to configure the ``job_config.yaml`` file to set up a job
runner that uses the DRMAA library for job submission.


**job_config.yaml**
::
   
  runners:
    local:
      load: galaxy.jobs.runners.local:LocalJobRunner
    slurm:
      load: galaxy.jobs.runners.slurm:SlurmJobRunner
      drmaa_library_path: PATH/slurm-drmaa-1.1.4/dist/lib/libdrmaa.so

  execution:
    default: slurm
    environments:
      local:
        runner: local
      slurm:
        runner: slurm
        native_specification: "-p PARTITION_NAME --time=02:00:00 --nodes=1"
        singularity_enabled: true
        singularity_volumes: $defaults
  tools:
    - class: local
      environment: local

This example runs the workflow on a SLURM cluster, using Singularity to resolve tool dependencies. However, it
can be easily adapted to use Docker or Conda instead. For more details on configuring the ``job_config.yaml`` file,
refer to the `job configuration documentation <https://docs.galaxyproject.org/en/master/admin/jobs.html>`__

**Note:** before using the example remember to update ``PATH`` and ``PARTITION_NAME`` values to match your system’s configuration.

Run using the default SQLite database
^^^^^^^^^^^^^^
:: 
  
  planemo run tutorial.ga tutorial-job.yml \
    --download_outputs \
    --output_directory . \
    --output_json output.json \
    --job_config_file job_config.yaml \
    --biocontainers

Run using a temporary PostgreSQL database
^^^^^^^^^^^^^^
:: 
  
  planemo run tutorial.ga tutorial-job.yml \
    --download_outputs \
    --output_directory . \
    --output_json output.json \
    --job_config_file job_config.yaml \
    --biocontainers \
    --database_type postgres_singularity

Slurm with TPV
~~~~~~~~~~~~~~~~~~
To make it easier to configure resource used by the different tools we can 
use `TPV (Total Perspective Vortex) <https://total-perspective-vortex.readthedocs.io/en/latest/>`_ to
set the resources used by each individual tool. This well require some changes to the 
``job_config.yaml``.

**job_config.yaml**
::

  runners:
    local:
      load: galaxy.jobs.runners.local:LocalJobRunner

    drmaa:
      load: galaxy.jobs.runners.drmaa:DRMAAJobRunner
      drmaa_library_path: PATH/slurm-drmaa-1.1.4/dist/lib/libdrmaa.so

  execution:
    default: tpv
    environments:
      local:
        runner: local
      tpv:
        runner: dynamic
        function: map_tool_to_destination
        rules_module: tpv.rules
        tpv_config_files:
        - https://gxy.io/tpv/db.yml
        - PATH/destinations.yaml

  tools:
    - class: local
      environment: local

**destinations.yaml**
::

  destinations:
    tpvdb_drmaa:
      runner: drmaa
      params:
        native_specification: "-p PARTION_NAME --time=02:00:00 --nodes=1  --ntasks={cores} --ntasks-per-node={cores}"
        singularity_enabled: true
        singularity_volumes: $defaults

Like the previous example, this setup submits jobs to a SLURM cluster and uses Singularity to resolve dependencies. The
key difference is the use of **TPV**, which allows you to define resource requirements—such as the number of cores used
by bwa mem—based on settings from the shared configuration a ``https://gxy.io/tpv/db.yml``

These defaults can be customized either by providing an entirely new db.yaml file or by overriding specific tool settings
in a separate YAML file. For more details, see the , see `using-the-shared-database <https://total-perspective-vortex.readthedocs.io/en/latest/topics/tpv_by_example.html#using-the-shared-database>`__ section of the TPV documentation.

**Note:** even if your SLURM scheduler does not use ntasks, you should still set it when a tool is intended to use more
than one core. If not specified, Galaxy will default to using a single core for that tool.

Troubleshooting
------------------

**Temp direcory not shared between nodes**

Galaxy uses a temporary directory when creating and running jobs.
This directory must be accessible to both the server that creates
the job and the compute node that executes it. If this folder is
not shared, local jobs will succeed, while those submitted to a
separate compute node will fail. To resolve this issue, configure
a shared folder as the temporary directory.

::
  
  TMPDIR=PATH_TO_SHARED_TEMPORARY_FOLDER planemo run ...

**Locked database**

While running jobs, Galaxy tracks their status using a database. 
The default database is SQLite, which may encounter issues when 
handling a high number of concurrent jobs. In such cases, you may see errors like this:


::
  
  cursor.execute(statement, parameters)
  sqlalchemy.exc.OperationalError: (sqlite3.OperationalError) database is locked
  [SQL: SELECT kombu_queue.id AS kombu_queue_id, kombu_queue.name AS kombu_queue_name 
  FROM kombu_queue 
  WHERE kombu_queue.name = ?
    LIMIT ? OFFSET ?]

The solution is to switch to a more robust database, such as PostgreSQL.

**Multiple galaxy instances**

When you run ``planemo run`` a temporary Galaxy instance is created and 
in some cases that instance my not be properly shut down. This can cause 
the following error:

:: 

  bioblend.ConnectionError: Unexpected HTTP status code: 500: Internal Server Error: make sure that you don't already have a running galaxy instance.

You will be able to see if instances are still running using the following commands:

::
  
    # Find running galaxy instances
    ps aux | grep galaxy
    # Find running planemo commands
    ps aux | grep planemo
