
``bioc_tool_init`` command
======================================

This section is auto-generated from the help text for the planemo command
``bioc_tool_init``. This help message can be generated with ``planemo bioc_tool_init
--help``.

**Usage**::

    planemo bioc_tool_init [OPTIONS]

**Help**

Generate a bioconductor tool outline from supplied arguments.
**Options**::


      -i, --id TEXT             Short identifier for new tool (no whitespace)
      -f, --force               Overwrite existing tool if present.
      -t, --tool PATH           Output path for new tool (default is <id>.xml)
      -n, --name TEXT           Name for new R/Bioconductor tool (user facing).
      -d, --description TEXT    Short description for new tool (user facing)
      -c, --command TEXT        Command potentially including cheetah variables
                                ()(e.g. 'seqtk seq -a $input > $output')
      --rscript PATH            Name of an R script from which to create a Tool
                                definition file. Requires use of --input and
                                --output arguments. (e.g. --rscript 'file.R')
      --rversion TEXT           R version this tool requries, if not given, the tool
                                defaults to R 3.2.1, (eg: --rversion 'R 3.2.1').
                                This option adds the R dependency in the tool
                                requirements.
      --input TEXT              An input description (e.g. input.fasta)
      --output TEXT             An output location (e.g. output.bam), the Galaxy
                                datatype is inferred from the extension.
      --requirement TEXT        Name of the R/Bioconductor package. Requirements
                                will be set using Bioconda. (e.g. --requirement
                                'affy')
      --help_text TEXT          Help text (reStructuredText)
      --doi TEXT                Supply a DOI (http://www.doi.org/) easing citation
                                of the tool for Galxy users (e.g. 10.1101/014043).
      --cite_url TEXT           Supply a URL for citation.
      --version TEXT            Tool XML version.
      --help_from_command TEXT  Auto populate help from supplied command.
      --test_case               For use with --example_commmand, generate a tool
                                test case from the supplied example.
      --macros                  Generate a macros.xml for reuse across many tools.
      --named_output TEXT       Create a named output for use with command block for
                                example specify --named_output=output1.bam and then
                                use '-o $output1' in your command block.
      --bioconda_path TEXT      Path to bioconda repository. If left empty, path
                                will be made in home directory.
      --help                    Show this message and exit.
    
