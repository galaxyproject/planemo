
``output_schema`` command
========================================

This section is auto-generated from the help text for the planemo command
``output_schema``. This help message can be generated with ``planemo output_schema
--help``.

**Usage**::

    planemo output_schema [OPTIONS]

**Help**

Export JSON Schemas for Planemo machine-readable outputs.

Available schemas:

cli-command-metadata: output from `planemo cli_metadata --command NAME`.
cli-metadata: output from `planemo cli_metadata`.
invocation-download-manifest: output from `planemo invocation_download --output_json PATH`.
run-outputs: output from `planemo run --output_json PATH`.
test-report: output from test report JSON writers such as `planemo test --test_output_json PATH`.
**Options**::


      --schema [cli-command-metadata|cli-metadata|invocation-download-manifest|run-outputs|test-report]
                                      Only export one schema.
      --help                          Show this message and exit.
    
