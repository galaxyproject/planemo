Galaxy Workflows on a local cluster
============================================

Some users may be unable to upload data to publicly available Galaxy
instances due to legal restrictions. While they might have access to
locally hosted cluster systems, these systems may not allow hosting
web servers or have the necessary personnel to maintain a Galaxy
instance. In such cases, Planemo provides a solution by spinning up
a temporary local Galaxy instance, allowing users to execute
workflows. Jobs can then be run on the same server or submitted to a
workload manager, such as slurm.

Software requirements
---------------------
Planemo
~~~~~~~~~~~~~~~~~~

If planemo isn't available on the system, it can be easily installed
in a virtual environment.
::
   
   python -m venv planemo_venv
   . planemo_venv/bin/activate
   pip install planemo

Additional instructions to install planemo can be found at
`install planemo <https://planemo.readthedocs.io/en/latest/installation.html>`__.

DRMAA library
~~~~~~~~~~~~~~~~~~
To run workflows on a cluster, Galaxy requires access to DRMAA libraries,
which allow communication with the workload manager.
If your system does not include DRMAA libraries, you will need to compile
them manually. For SLURM, the source code and compilation instructions
are available on GitHub:: `slurm-drmaa <https://github.com/natefoo/slurm-drmaa>`__.

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
be available on the system. If they are not, Galaxy can manage dependencies
using Conda, Docker or Singularity. While Conda is an option, Docker or 
Singularity is recommended for its greater reliability. One of the following
environments need to be available on the system:

- Complete set of used tools installed in the system
- Conda
- Docker
- Singularity


Database
------------------
Galaxy uses a database to track jobs and their statuses, with SQLite as the default 
option. However, when running more complex workflows, errors can occur if multiple 
jobs are updated simultaneously, leading to database locking and causing the workflow 
to terminate prematurely. This issue can be easily avoided by switching to a 
PostgreSQL database instead of SQLite. If the user doesn't have access to a PostgreSQL
database, planemo can spin up a temporary PostgreSQL instance using Singularity by adding
a ``--database_type postgres_singularity`` to the run comman.

Configuration
------------------
To apply the settings needed for running a workflow/tool you need to configure a job_conf.yaml
file, some of the settings can be applied using the command line, for example:

- use temporary postgres database: ``--database_type postgres_singularity``
- use containers for tools: ``--biocontainers``

Job configuration
~~~~~~~~~~~~~~~~~~

The job configuration file must specify the path to the DRMAA library and define basic
SLURM settings using the `native_specification` tag. For convenience, we will also 
configure the Singularity options.

::
   
  runners:
    local:
      load: galaxy.jobs.runners.local:LocalJobRunner
      # modify the number of threads working on local jobs here
      # workers: 4
    slurm:
      load: galaxy.jobs.runners.slurm:SlurmJobRunner
      # modify below to specify a specific drmaa shared library for slurm
      drmaa_library_path: /scratch/10171/smeds/run_planemo_using_drmaa/slurm-drmaa-1.1.4/dist/lib/libdrmaa.so

  execution:
    default: slurm
    environments:
      local:
        runner: local
        singularity_enabled: true
        singularity_cmd: singularity
        singularity_volumes: $defaults
      slurm:
        runner: slurm
        native_specification: "-p icx,skx --time=02:00:00 --nodes=1"
        singularity_enabled: true
        singularity_cmd: singularity
        singularity_volumes: $defaults
        require_container: true
  tools:
    - class: local
      environment: local


If your system has a different setup, Galaxy provides `galaxy-job-config-init <https://pypi.org/project/galaxy-job-config-init/>`__
, a tool for generating job_config files tailored to various environments and workload managers.
More examples of how to setup the job_configuration file can be found at `job_conf.sample.yaml <https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/config/sample/job_conf.sample.yml>`__.

Running 
------------------
Planemo can be run with many different settings, depending on the user's needs. This example
will be running a workflow using slurm and SQLite or PostgreSQL.

Example workflow
~~~~~~~~~~~~~~~~~~
Here is the ``job_conf.yaml`` file for the example workflow. It has been modified
to handle two specific tools differently, as they require Conda to run.

::
  
  runners:
    local:
      load: galaxy.jobs.runners.local:LocalJobRunner
      # modify the number of threads working on local jobs here
      # workers: 4
    slurm:
      load: galaxy.jobs.runners.slurm:SlurmJobRunner
      # modify below to specify a specific drmaa shared library for slurm
      drmaa_library_path: PATH_TO_LIBRARY/slurm-drmaa-1.1.4/dist/lib/libdrmaa.so

  execution:
    default: slurm
    environments:
      local:
        runner: local
        singularity_enabled: true
        singularity_cmd: singularity
        singularity_volumes: $defaults
      slurm:
        runner: slurm
        native_specification: "-p QUEUE_NAME --time=02:00:00 --nodes=1"
        singularity_enabled: true
        singularity_cmd: singularity
        singularity_volumes: $defaults
        require_container: true
      slurm_conda:
        runner: slurm
        native_specification: "-p QUEUE_NAME --time=02:00:00 --nodes=1"
        conda_enabled: true
        require_container: false
  tools:
    - id: cat1
      environment: slurm_conda
    - id: random_lines1
      environment: slurm_conda
    - class: local
      environment: local


Run using the default SQLite database
:: 
  
  git clone https://github.com/usegalaxy-eu/workflow-testing.git
  cd workflow-testing/example3
  python3.9 -m venv planemo_venv && source planemo_venv/bin/activate
  pip install planemo
  # Make sure the temporary folder is shared between servers.
  planemo run tutorial.ga tutorial-job.yml \
    --download_outputs \
    --output_directory . \
    --output_json output.json \
    --job_config_file job_config.yaml

Run using a temporary PostgreSQL database
:: 
  
  git clone https://github.com/usegalaxy-eu/workflow-testing.git
  cd workflow-testing/example3
  python3.9 -m venv venv && source venv/bin/activate
  pip install planemo
  # Make sure that the singulariy command is available
  # Make sure the temporary folder is shared between servers.
  planemo run tutorial.ga tutorial-job.yml \
    --download_outputs \
    --output_directory . \
    --output_json output.json \
    --job_config_file job_config.yaml \
    --database_type postgres_singularity


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
