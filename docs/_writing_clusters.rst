Cluster Usage
==============================

-----------------------------------------------------------------------------------------------------
Developing for Clusters - ``GALAXY_SLOTS``, ``GALAXY_MEMORY_MB``, and ``GALAXY_MEMORY_MB_PER_SLOT``
-----------------------------------------------------------------------------------------------------

``GALAXY_SLOTS`` is a special environment variable that is set in a Galaxy
tool's runtime environment. If the tool you are working on allows configuring
the number of processes or threads that should be spawned, this variable
should be used.

For example, for the `StringTie tool
<https://github.com/galaxyproject/tools-iuc/blob/main/tools/stringtie/stringtie.xml>`__
the binary ``stringtie`` can take an argument ``-p`` that allows specification
of the number of threads to be used. The Galaxy tool sets this up as follows

::

    stringtie "$input_bam" -o "$output_gtf" -p "\${GALAXY_SLOTS:-1}" ...

Here we use ``\${GALAXY_SLOTS:-Z}`` instead of a fixed value (Z being an
integer representing a default value in non-Galaxy contexts). The
backslash here is because this value is interpreted at runtime as
environment variable - not during command building time as a templated
value. Now server administrators can configure how many processes the
tool should be allowed to use.

For information on how server administrators can configure this value for
a particular tool, check out `the Galaxy admin documentation
<https://galaxyproject.org/admin/config/galaxy_slots/>`__.

Analogously ``GALAXY_MEMORY_MB`` and ``GALAXY_MEMORY_MB_PER_SLOT`` are special 
environment variables in a Galaxy tool's runtime environment that can be used
to specify the amount of memory that a tool can use overall and per slot,
respectively. 

For an example see the `samtools sort tool <https://github.com/galaxyproject/tools-iuc/blob/main/tool_collections/samtools/samtools_sort/samtools_sort.xml>`__ which allows to specify the
total memory with the -m parameter.

-----------------------------------------------
Test Against Clusters - ``--job_config_file``
-----------------------------------------------

The various commands that start Galaxy servers (``serve``,
``test``, ``shed_serve``, ``shed_test``, etc...) allow specification of
a Galaxy job configuration XML file (e.g. ``job_conf.xml``).

For instance, Slurm_ is a popular distributed reource manager (DRM) in the
Galaxy community. The following ``job_conf.xml`` tells Galaxy to run all jobs
using Slurm_ and allocate ``2`` cores for each job.

.. literalinclude:: writing/job_conf_slurm.xml
   :language: xml

If this file is named ``planemo_job_conf.xml`` and resides in one's home
directory (``~``), Planemo can ``test`` or ``serve`` using this configuration
with the following commands.

::

    $ planemo test --job_config_file ~/planemo_job_conf.xml .
    $ planemo serve --job_config_file ~/planemo_job_conf.xml .

For general information on configuring Galaxy to communicate with clusters
check out `this page
<https://wiki.galaxyproject.org/Admin/Config/Performance/Cluster>`__ on the
Galaxy wiki and for information regarding configuring job configuration XML
files in particular check out `the example
<https://github.com/galaxyproject/galaxy/blob/dev/config/job_conf.xml.sample_advanced>`__
distributed with Galaxy.

--------------------------------------------
Generated Job and TPV Configurations
--------------------------------------------

When Planemo starts a Galaxy instance for testing or serving and no
``--job_config_file`` is given, it generates a ``job_conf.yml`` from the other options -
the job runner, and whether and how containers are enabled. ``planemo job_config_init`` and
``planemo profile_job_config_init`` write such a configuration into the current directory
or into a profile instead, from where it can be reviewed and reused with
``--job_config_file``.

Mapping individual tools to destinations is what `Total Perspective Vortex (TPV)`_ does.
``--tpv`` adds it to the generated configuration: the job configuration then references the
shared Galaxy TPV database at ``https://gxy.io/tpv/db.yml`` plus a single inline destination
(``tpvdb_local``, ``tpvdb_slurm``, ...) that carries the runner and container settings of
the run. That database assigns tools to destinations of its own, which the site that runs
the jobs has to provide - it is what :ref:`writing_slurm` uses to give jobs resource
requirements.

``--tpv_config_file``, which implies ``--tpv``, does not bring in the shared database.
Instead Planemo generates two files next to the job configuration of the Galaxy instance it
manages (with ``planemo job_config_init``: in the current directory):

``job_conf.yml``
    the Galaxy job configuration with a single ``tpv`` environment, which references the
    generated configuration below and the files given with ``--tpv_config_file``, and

``planemo_tpv.yml``
    one destination ``planemo_default`` that takes over the runner, container, and volume
    settings of the job configuration, and a ``default`` tool entity that all tools not
    mentioned elsewhere inherit from.

``--tpv_config_file`` can be passed several times. The files are layered on top of the
generated configuration in the order they are given, and a configuration only changes the
tools it mentions. What Planemo generates is self-contained: it references no shared
database, so mapping jobs needs no network access and every tool is assigned to a
destination that really exists. It also sets no resource requirements of its own - add them
to the ``default`` tool entity of your own configuration, or download the shared Galaxy TPV
database (``curl -o db.yml https://gxy.io/tpv/db.yml``) and pass the local file with
``--tpv_config_file``.

TPV is a Galaxy dependency rather than a Planemo one. Planemo installs
`total-perspective-vortex`_ into the virtualenv of the Galaxy instance it manages when it is
missing, so the first run with ``--tpv`` takes a bit longer. Only Galaxy instances managed
by Planemo are affected - with ``--galaxy_url`` the job configuration of the target instance
applies.

.. _Slurm: http://slurm.schedmd.com/
.. _Total Perspective Vortex (TPV): https://total-perspective-vortex.readthedocs.io/
.. _total-perspective-vortex: https://pypi.org/project/total-perspective-vortex/
