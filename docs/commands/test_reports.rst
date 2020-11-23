
``test_reports`` command
======================================

This section is auto-generated from the help text for the planemo command
``test_reports``. This help message can be generated with ``planemo test_reports
--help``.

**Usage**::

    planemo test_reports [OPTIONS] FILE_PATH

**Help**

Generate human readable tool test reports.

Creates reports in various formats  (HTML, text, markdown)
from the structured test output (tool_test_output.json).

**Options**::


      --test_output PATH           Output test report (HTML - for humans) defaults
                                   to tool_test_output.html.
    
      --test_output_text PATH      Output test report (Basic text - for display in
                                   CI)
    
      --test_output_markdown PATH  Output test report (Markdown style - for humans &
                                   computers)
    
      --test_output_xunit PATH     Output test report (xunit style - for CI systems
      --test_output_junit PATH     Output test report (jUnit style - for CI systems
      --help                       Show this message and exit.
    
