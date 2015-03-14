How do I...
====================================================

This section contains a number of smaller topics with links and examples meant
to provide relatively concrete answers for specific tool development
scenarios.

------------------------------------------
\.\.\. deal with index/reference data?
------------------------------------------

Galaxy's concept of `data tables
<https://wiki.galaxyproject.org/Admin/Tools/Data%20Tables>`__ are meant to
provide tools with access reference datasets or index data not tied to
particular histories or users. A common example would be FASTA files for
various genomes or mapper-specific indices of those files (e.g. a BWA index
for the hg19 genome).

Galaxy `data managers
<https://wiki.galaxyproject.org/Admin/Tools/DataManagers>`__ are specialized
tools designed to populate tool data tables.


------------------------------------------
\.\.\. cite tools without an obvious DOI?
------------------------------------------

In the absence of an obvious DOI_, tools may contain embedded BibTeX_ directly.

Futher reading:

- `bibtex.xml <https://github.com/jmchilton/galaxy/blob/dev/test/functional/tools/bibtex.xml>`__ (test tool with a bunch of random examples)
- `bwa-mem.xml <https://github.com/jmchilton/bwa-mem/commit/0425264039950bfd9ded06997a08cc8b4ee1ad8f>`__ (BWA-MEM tool by Anton Nekrutenko demonstrating citation of an arXiv article)

--------------------------------------------------
\.\.\. declare a Docker container for my tool?
--------------------------------------------------

Galaxy tools can be decorated to with ``container`` tags indicated Docker
container ids that the tools can run inside of.

The longer term plan for the Tool Shed ecosystem is to be able to
automatically build Docker containers for tool dependency descriptions and
thereby obtain this Docker functionality for free and in a way that is
completely backward compatible with non-Docker deployments.

Further reading:

- `Complete tutorial <https://github.com/apetkau/galaxy-hackathon-2014>`__
  on Github by Aaron Petkau. Covers installing Docker, building a Dockerfile_, publishing to `Docker Hub`_, annotating tools and configuring Galaxy.
- `Another tutorial <https://www.e-biogenouest.org/groups/guggo>`__
  from the Galaxy User Group Grand Ouest.
- Landing page on the `Galaxy Wiki <https://wiki.galaxyproject.org/Admin/Tools/Docker>`__
- Impementation details on `Pull Request #401 <https://bitbucket.org/galaxy/galaxy-central/pull-request/401/allow-tools-and-deployers-to-specify>`__

--------------------------------------------------
\.\.\. do extra validation of parameters?
--------------------------------------------------

Tool parameters support a ``validator`` element (`syntax
<https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax#A.3Cvalidator.3E_tag_set>`__)
to perform validation of a single parameter. More complex validation across
parameters can be performed using arbitrary Python functions using the
``code`` file syntax but this feature should be used sparingly.

Further reading:

- `validator <https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax#A.3Cvalidator.3E_tag_set>`__
  XML tag syntax on the Galaxy wiki.
- `fastq_filter.xml <https://github.com/galaxyproject/tools-devteam/blob/master/tool_collections/galaxy_sequence_utils/fastq_filter/fastq_filter.xml>`__
  (a FASTQ filtering tool demonstrating validator constructs)
- `gffread.xml <https://github.com/galaxyproject/tools-devteam/blob/master/tool_collections/cufflinks/gffread/gffread.xml>`__
  (a tool by Jim Johnson demonstrating using regular expressions with ``validator`` tags)
- `code_file.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/code_file.xml>`__,
  `code_file.py <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/code_file.py>`__
  (test files demonstrating defining a simple constraint in Python across
  two parameters)
- `deseq2 tool <https://github.com/bgruening/galaxytools/tree/master/tools/deseq2>`__
  by Björn Grüning demonstrating advanced ``code`` file validation.

.. _DOI: http://www.doi.org/
.. _BibTeX: http://www.bibtex.org/
.. _Dockerfile: https://docs.docker.com/reference/builder/
.. _Docker Hub: https://hub.docker.com/