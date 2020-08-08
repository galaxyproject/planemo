# {{ title }}

## Test Summary

| Test State | Count |
| ---------- | ----- |
| Total      | {{ raw_data.results.total | default(0)  }} |
| Passed     | {{ raw_data.results.total - raw_data.results.errors - raw_data.results.failures - raw_data.results.skips | default(0)  }} |
| Error      | {{ raw_data.results.errors | default(0) }} |
| Failure    | {{ raw_data.results.failures | default(0) }} |
| Skipped    | {{ raw_data.results.skipped | default(0) }} |


{% for state, desc in {'error': 'Errored', 'failure': 'Failed', 'success': 'Passed'}.items() %}
<details><summary>{{ desc }} Tests</summary>
{%   for test in raw_data.tests %}
{%     if test.data.status == state %}
{%       if test.data.status == 'success' %}

* <details><summary>:white_check_mark: {{ test.id }}</summary>

{%       else %}

* <details><summary>:x: {{ test.id }}</summary>

{%       endif %}

    #### Problems
{%       for problem in test.data.output_problems %}
    ```console
    {{problem|indent}}
    ```
{%       endfor %}
{%     if test.data.job %}

    #### Command Line:

    ```console
    {{ test.data.job.command_line|indent}}
    ```
    exited with code {{ test.data.job.exit_code }}.

 {%    if test.data.job.stdout %}

    #### `stderr`

    ```
    {{ test.data.job.stderr|indent }}
    ```

{%     endif %}
{%     if test.data.job.stdout %}

    #### `stdout`

    ```console
    {{ test.data.job.stdout|indent }}
    ```
{%     endif %}
{%     endif %}
{%     if test.data.invocation_details %}

    #### Workflow invocation details

    <details><summary>Steps</summary>

{%    for step_data in test.data.invocation_details.values() %}
    {{step_data.order_index}}. **{{step_data.workflow_step_label or (step_data.jobs[0].tool_id if step_data.jobs[0] else 'Unlabelled step')}}**:
      step_state: {{step_data.state}}
{%       if step_data.jobs %}
        <details><summary>jobs:</summary>
{%       for job in step_data.jobs %}
        - job {{loop.index}}:

         | Job property | Value |
         | ------------ | ----- |
{%          for key, value in job.items() %}
{%          if value %}
         | {{key}} | `{{value}}` |
{%          endif %}
{%          endfor %}
{%       endfor %}
        </details>
{%       endif %}
  </details>
{%       endfor %}
  </details>
</details>
{% endif %}
{% endif %}
{% endfor %}

</details>
{% endfor %}
