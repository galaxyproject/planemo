Cluster Usage
==============================

------------------------------
``GALAXY_SLOTS``
------------------------------

``GALAXY_SLOTS`` is a special environment variable that is set in a Galaxy
tool's runtime environment. If the tool you are working on allows configuring
the number of processes or threads that should be spawned, this variable
should be used.

For example, the StringTie (tool available `here
<https://github.com/galaxyproject/tools-iuc/blob/master/tools/stringtie/stringtie.xml>`__)
binary ``stringtie`` can take an argument ``-p`` that allows specification
of the number of threads to be used. The Galaxy tool sets this up as follows::

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

.. _stringtie: https://ccb.jhu.edu/software/stringtie/
