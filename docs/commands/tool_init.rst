
``tool_init`` command
===============================

This section is auto-generated from the help text for the planemo command
``tool_init``. This help message can be generated with ``planemo tool_init
--help``.

**Usage**::

    planemo tool_init [OPTIONS]

**Help**

Generate a tool outline from supplied arguments.

**Options**::


      --test_case               For use with --example_commmand, generate a tool
                                test case from the supplied example.
      --container TEXT          Add a Docker image identifier for this tool.
      --requirement TEXT        Add a tool requirement package (e.g. 'seqtk' or
                                'seqtk@1.68').
      --help_from_command TEXT  Auto populate help from supplied command.
      --help_text TEXT          Help text (reStructuredText)
      --named_output TEXT       Create a named output for use with command block
                                for example specify --named_output=output1.bam and
                                then use '-o $ouput1' in your command block.
      --output TEXT             An output location (e.g. output.bam), the Galaxy
                                datatype is inferred from the extension.
      --input TEXT              An input description (e.g. input.fasta)
      --version_command TEXT    Command to print version (e.g. 'seqtk --version')
      --example_output TEXT     For use with --example_command, replace input file
                                (e.g. 2.fastq with a tool output).
      --example_input TEXT      For use with --example_command, replace input file
                                (e.g. 2.fastq with a data input parameter).
      --example_command TEXT    Example to command with paths to build Cheetah
                                template from (e.g. 'seqtk seq -a 2.fastq >
                                2.fasta'). Option cannot be used with --command,
                                should be used --example_input and
                                --example_output.
      -c, --command TEXT        Command potentially including cheetah variables
                                ()(e.g. 'seqtk seq -a $input > $output')
      -d, --description TEXT    Short description for new tool (user facing)
      --version TEXT            Tool XML version.
      -n, --name TEXT           Name for new tool (user facing)
      -t, --tool PATH           Output path for new tool (default is <id>.xml)
      -f, --force               Overwrite existing tool if present.
      -i, --id TEXT             Short identifier for new tool (no whitespace)
      --help                    Show this message and exit.
    
