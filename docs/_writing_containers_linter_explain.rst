
This last linter indicates that indeed a container has been registered
that is compatible with this tool -- ``quay.io/biocontainers/seqtk:1.2--1``.
We didn't do any extra work to build this container for this tool, all
Bioconda_ recipes are packaged into containers and registered on quay.io_
as part of the BioContainers_ project.

This tool can be tested using ``planemo test`` in its BioContainer
Docker container using the flag ``--biocontainers`` as shown below.

.. _BioContainers: http://biocontainers.pro/
.. _quay.io: https://quay.io
.. _Bioconda: https://github.com/bioconda/bioconda-recipes
