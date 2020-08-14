{% from 'macros.tmpl' import render_steps %}
{% if title %}
# {{ title }}

{% endif %}
## Test Summary
{% set state = namespace(found=false) %}
{% set state.success = raw_data.results.total - raw_data.results.errors - raw_data.results.failures - raw_data.results.skips | default(0) %}
{% set state.error = raw_data.results.errors | default(0) %}
{% set state.failure = raw_data.results.failures | default(0) %}
{% set state.skipped = raw_data.results.skipped | default(0) %}

| Test State | Count |
| ---------- | ----- |
| Total      | {{ raw_data.results.total | default(0)  }} |
| Passed     | {{ state.success }} |
| Error      | {{ state.error }} |
| Failure    | {{ state.failure }} |
| Skipped    | {{ state.skipped }} |


{% set display_job_attributes = {'command_line': 'Command Line', 'exit_code': 'Exit Code', 'stderr': 'Standard Error', 'stdout': 'Standard Output'} %}
{% for status, desc in {'error': 'Errored', 'failure': 'Failed', 'success': 'Passed'}.items() if state[status]%}
<details><summary>{{ desc }} Tests</summary>
{%   for test in raw_data.tests %}
{%     if test.data.status == status %}
{%       if test.data.status == 'success' %}

* <details><summary>&#9989; {{ test.id }}</summary>

{%       else %}

* <details><summary>&#10060; {{ test.id }}</summary>

{%       endif %}
{%       if test.data.output_problems %}
    #### Problems
{%       endif %}
{% if test.data.execution_problem %}
    ##### Execution Problem
    ```
    {{test.data.execution_problem|indent}}
    ```
{% endif %}
{%       for problem in test.data.output_problems %}
    ```
    {{problem|indent}}
    ```
{%       endfor %}
{%     if test.data.job %}
{%       for key, description in display_job_attributes.items() %}
{%         if test.data.job[key] not in ("", None) %}
    #### {{ description }}:

    ```console
    {{ test.data.job[key]|string|indent}}
    ```
{%         endif %}
{%       endfor %}
{%       if test.data.job.params %}
   #### Job Parameters:

     | Job parameter | Parameter value |
     | ------------- | --------------- |
{%         for key, value in test.data.job.params.items() %}
     | {{ key }} | ` {{ value }} ` |
{%         endfor %}
{%       endif %}

{%     endif %}
{%     if test.data.invocation_details %}

    #### Workflow invocation details

{{render_steps(test.data.invocation_details.values(), display_job_attributes)}}

{% endif %}
{% endif %}
{% endfor %}

</details>
{% endfor %}
