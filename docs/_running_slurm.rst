Running workflows on a server or Slurm cluster
==============================================

Legal, policy, or storage constraints can keep data off a public Galaxy
instance. Planemo can instead launch a Galaxy instance on a server or cluster
login node and run workflow jobs there, either locally or through a workload
manager such as Slurm. The data stays in the local environment and no permanent
Galaxy service has to be maintained.

This guide uses the ``tutorial.ga`` workflow and ``tutorial-job.yml`` job file
introduced in `The Basics`_ and covers:

#. running the workflow on one server;
#. reusing its database and installed tools through a Planemo profile;
#. generating a Slurm job configuration; and
#. submitting jobs through Slurm with or without TPV.

Requirements
------------

The server or cluster needs:

* Planemo, installed as described in :doc:`readme`;
* a Slurm client configured for the cluster;
* a DRMAA implementation for Slurm;
* Conda, Docker, or Apptainer for tool dependencies; and
* a filesystem visible to both the login node and compute nodes.

Slurm and DRMAA
~~~~~~~~~~~~~~~

Galaxy's Slurm runner uses the Distributed Resource Management Application API
(DRMAA). Ask the cluster administrators whether a DRMAA library is already
available. If it is not, Planemo can download and build the supported
`slurm-drmaa`_ release:

.. code-block:: console

    $ planemo slurm_init

The command copies ``libdrmaa.so`` into the Planemo workspace, which defaults
to ``~/.planemo/libdrmaa.so``. It requires a working compiler toolchain and the
Slurm development files. See :doc:`commands/slurm_init` for command details.

Note the absolute path to this or an administrator-provided ``libdrmaa``
library. You will add it to the generated job configuration below.

.. _slurm-drmaa: https://github.com/natefoo/slurm-drmaa

Tool dependencies
~~~~~~~~~~~~~~~~~

Galaxy can resolve tool dependencies with Conda or run tools in Docker or
Apptainer containers. Containers are generally more reproducible, but the
available runtime and mount configuration depend on the cluster. The examples
below use Apptainer because it is commonly available on HPC systems.

Apptainer was previously called Singularity, and Planemo still uses the older
name for the option and for the generated configuration keys. Select the runtime
with ``--singularity``. The generated configuration runs
``singularity_cmd: singularity``; set that to ``apptainer`` on clusters that do
not also provide a ``singularity`` executable.

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

Reuse Galaxy state with a profile
---------------------------------

Each bare ``planemo run`` starts Galaxy from scratch: an empty database, and a
fresh installation of whatever Tool Shed tools the workflow needs, into a
directory discarded along with the run. A *profile* keeps the database, the
installed tools, and their data tables and data managers in the Planemo
workspace instead, so later runs reuse them. Create one before going further.

.. code-block:: console

    $ planemo profile_create slurm_cluster
    Profile [slurm_cluster] created.

Pass the profile to subsequent runs:

.. code-block:: console

    $ planemo run tutorial.ga tutorial-job.yml \
        --profile slurm_cluster \
        --download_outputs \
        --output_directory . \
        --output_json output.json

The profile lives in ``~/.planemo/profiles/slurm_cluster``. ``planemo
profile_list`` lists existing profiles and ``planemo profile_delete`` removes
one along with its database.

Run through Slurm
-----------------

Generate the job configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than maintain a complete Galaxy job configuration by hand, generate one
into the profile with :doc:`commands/profile_job_config_init`:

.. code-block:: console

    $ planemo profile_job_config_init slurm_cluster \
        --runner slurm \
        --singularity

Runs that pass ``--profile slurm_cluster`` pick the configuration up
automatically. Without a profile, :doc:`commands/job_config_init` takes the same
options and writes ``job_conf.yml`` into the current directory, which
``planemo run --job_config_file job_conf.yml`` then reads.

Review the generated file and adjust the Slurm runner and environment:

* uncomment ``drmaa_library_path`` and set it to the absolute path of the
  cluster's ``libdrmaa`` library;
* set ``native_specification`` to the partition, time, memory, and other Slurm
  options required by the cluster; and
* change the generated ``singularity_cmd``, volumes, or environment variables
  where the cluster requires it.

The generator sets ``tmp_dir: true``, which asks Galaxy to manage a temporary
directory per job and usually removes any need to set a shared ``TMPDIR``.
Galaxy job working directories and input data must still live on storage the
compute nodes can see.

Run the workflow
~~~~~~~~~~~~~~~~

.. code-block:: console

    $ planemo run tutorial.ga tutorial-job.yml \
        --profile slurm_cluster \
        --download_outputs \
        --output_directory . \
        --output_json output.json

Slurm with TPV
--------------

`Total Perspective Vortex (TPV)`_ maps individual tools to job destinations
and resource requirements. Planemo can generate the surrounding Galaxy and TPV
configuration while retaining the Slurm and container settings from the
previous example:

.. code-block:: console

    $ planemo profile_job_config_init slurm_cluster \
        --runner slurm \
        --tpv \
        --singularity

Planemo refuses to generate into a profile that already holds a job
configuration, so delete the profile's ``job_conf.yml`` first or keep the TPV
variant in a second profile.

In the generated ``job_conf.yml``, set ``drmaa_library_path`` and adapt the
``tpvdb_slurm`` destination to the local cluster. In particular its
``native_specification`` template has to translate TPV's ``{cores}``, ``{mem}``,
and other resource values into the Slurm options the site uses. The generated
configuration also pulls in the shared Galaxy TPV database at
``https://gxy.io/tpv/db.yml``, whose entries local configuration can override.

Make sure the native specification passes on the requested core count even when
the scheduler does not require ``--ntasks``. Otherwise a multithreaded tool can
be scheduled with a single core.

For configuration details and examples, see the `TPV documentation`_ and the
Galaxy Training Network's `job destination tutorial`_.

.. _Total Perspective Vortex (TPV): https://total-perspective-vortex.readthedocs.io/
.. _TPV documentation: https://total-perspective-vortex.readthedocs.io/en/latest/topics/tpv_by_example.html
.. _job destination tutorial: https://training.galaxyproject.org/training-material/topics/admin/tutorials/job-destinations/tutorial.html

Troubleshooting
---------------

Port already in use
~~~~~~~~~~~~~~~~~~~

If an earlier Planemo or Galaxy process still holds the configured port, stop
that process or pick another port with ``--port``. Identify leftover ``planemo``
and Galaxy processes with the cluster's usual process-monitoring tools before
terminating them.
