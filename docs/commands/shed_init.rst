
``shed_init`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_init``. This help message can be generated with ``planemo shed_init
--help``.

**Usage**::

    planemo shed_init [OPTIONS] PROJECT

**Help**

Bootstrap a new Tool Shed configuration (.shed.yml) file.

This file is used by other ``planemo`` commands such as ``shed_lint``,
``shed_create``, ``shed_upload``, and ``shed_diff`` to manage repositories
in a Galaxy Tool Shed.

**Options**::


      --from_workflow PATH            Attempt to generate repository dependencies
                                      from specified workflow.
      --description TEXT              Specify repository description for
                                      .shed.yml.
      --long_description TEXT         Specify repository long_description for
                                      .shed.yml.
      --remote_repository_url TEXT    Specify repository remote_repository_url for
                                      .shed.yml.
      --homepage_url TEXT             Specify repository homepage_url for
                                      .shed.yml.
      --category [Assembly|ChIP-seq|Combinatorial Selections|Computational chemistry|Convert Formats|Data Managers|Data Source|Fasta Manipulation|Fastq Manipulation|Genome-Wide Association Study|Genomic Interval Operations|Graphics|Imaging|Metabolomics|Metagenomics|Micro-array Analysis|Next Gen Mappers|NLP|Ontology Manipulation|Phylogenetics|Proteomics|RNA|SAM|Sequence Analysis|Statistics|Systems Biology|Text Manipulation|Tool Dependency Packages|Tool Generators|Transcriptomics|Variant Analysis|Visualization|Web Services]
                                      Specify repository category for .shed.yml
                                      (may specify multiple).
      --owner TEXT                    Tool Shed repository owner (username).
      --name TEXT                     Tool Shed repository name (defaults to the
                                      inferred tool directory name).
      -f, --force                     Overwrite existing files if present.
      --help                          Show this message and exit.
    
