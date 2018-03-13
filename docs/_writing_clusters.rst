Cluster Usage
==============================

-----------------------------------------------------------------------------------------------------
Developing for Clusters - ``GALAXY_SLOTS``, ``GALAXY_MEMORY_MB``, and ``GALAXY_MEMORY_MB_PER_SLOT``
-----------------------------------------------------------------------------------------------------

``GALAXY_SLOTS`` is a special environment variable that is set in a Galaxy
tool's runtime environment. If the tool you are working on allows configuring
the number of processes or threads that should be spawned, this variable
should be used.

For example, the StringTie (tool available `here
<https://github.com/galaxyproject/tools-iuc/blob/master/tools/stringtie/stringtie.xml>`__)
binary ``stringtie`` can take an argument ``-p`` that allows specification
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
a particular tool, check out `the Galaxy wiki
<https://wiki.galaxyproject.org/Admin/Config/GALAXY_SLOTS>`__.

Analogously ``GALAXY_MEMORY_MB`` and ``GALAXY_MEMORY_MB_PER_SLOT`` are special 
environment variables in a Galaxy tool's runtime environment that can be used
to specify the amount of memory that a tool can use overall and per slot, 
respectively. 

For an example see the samtools sort tool (`here https://github.com/galaxyproject/tools-iuc/blob/master/tool_collections/samtools/samtools_sort/samtools_sort.xml`__) which allows to specify the 
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

.. _Slurm: http://slurm.schedmd.com/
