# {{ title }}

## Executive Summary

| Test State | Count |
| ---------- | ----- |
| Total      | {{ raw_data.summary.num_tests | default(0)  }} |
| Passed     | {{ raw_data.summary.num_tests - raw_data.summary.num_errors - raw_data.summary.num_failures - raw_data.summary.num_skips | default(0)  }} |
| Error      | {{ raw_data.summary.num_errror | default(0) }} |
| Failure    | {{ raw_data.summary.num_failure | default(0) }} |
| Skipped    | {{ raw_data.summary.num_skipped | default(0) }} |


<details>
  <summary>
    <h2>Detailed Results</h2>
  </summary>
{% for test in raw_data.tests %}
{% if test.data.status == 'success' %}
### :white_check_mark: {{ test.id }}
{% else %}
### :x: {{ test.id }}
Test Error! (State: {{ test.data.status }})
#### Problems

{% for problem in test.data.output_problems %}
```console
{{problem}}
```
{% endfor %}

{%if test.data.job %}
Command Line:

```console
{{ test.data.job.command_line}}
```

exited with code {{ test.data.job.exit_code }}.

{% if test.data.job.stdout %}
#### `stderr`

```console
{{ test.data.job.stderr}}
```

{% endif %}
{% if test.data.job.stdout %}
#### `stdout`

```console
{{ test.data.job.stdout}}
```

{%- endif -%}
{%- endif -%}
{%- endif -%}

{%if test.data.invocation_details %}

#### Workflow invocation details

<details><summary><h2>Steps</h2></summary>
{%for step_data in test.data.invocation_details.values() %}
{{step_data.order_index}}. **{{step_data.workflow_step_label or (step_data.jobs[0].tool_id if step_data.jobs[0] else 'Unlabelled step')}}**:
  step_state: {{step_data.state}}
  {% if step_data.jobs %}
    <details><summary>jobs:</summary>
  {% for job in step_data.jobs %}
    - job {{loop.index}}:

     | Job property | Value |
     | ------------ | ----- |
     {% for key, value in job.items() %}
     {%- if value %}| {{key}} | `{{value}}` |
     {% endif -%}
     {%- endfor -%}
  {% endfor %}
    </details>
  {% endif %}
{%- endfor -%}
</details>
{%- endif -%}
{%- endfor %}
</details>
