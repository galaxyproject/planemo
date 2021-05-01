{% from 'macros.tmpl' import render_invocation_details, render_job_parameters, render_steps %}
{% if title %}
# {{ execution_type }} {{ title }}

{% endif %}
## {{ execution_type }} Summary
{% set state = namespace(found=false) %}
{% set state.success = raw_data.results.total - raw_data.results.errors - raw_data.results.failures - raw_data.results.skips | default(0) %}
{% set state.error = raw_data.results.errors | default(0) %}
{% set state.failure = raw_data.results.failures | default(0) %}
{% set state.skipped = raw_data.results.skipped | default(0) %}

{% if raw_data.results.total %}
<div class="progress">
  <div class="progress-bar progress-bar-success" style="width: {{ (state.success / raw_data.results.total) * 100 }}%" aria-valuenow="{{ state.success }}" aria-valuemin="0" aria-valuemax="{{ raw_data.results.total }}" data-toggle="tooltip" title="{{state.success}} Passed">
  </div>
  <div class="progress-bar progress-bar-warning" style="width: {{ (state.skipped / raw_data.results.total) * 100 }}%" aria-valuenow="{{ state.skipped }}" aria-valuemin="0" aria-valuemax="{{ raw_data.results.total }}" data-toggle="tooltip" title="{{state.skipped}} Skipped">
  </div>
  <div class="progress-bar progress-bar-danger" style="width: {{ ((state.error + state.failure) / raw_data.results.total) * 100 }}%" aria-valuenow="{{ state.error + state.failure }}" aria-valuemin="0" aria-valuemax="{{ raw_data.results.total }}" title="{{state.error + state.failure}} Failed or Errored">
  </div>
</div>
{% endif %}

| {{ execution_type }} State | Count |
| ---------- | ----- |
| Total      | {{ raw_data.results.total | default(0)  }} |
| Passed     | {{ state.success }} |
| Error      | {{ state.error }} |
| Failure    | {{ state.failure }} |
| Skipped    | {{ state.skipped }} |


{% set display_job_attributes = {'command_line': 'Command Line', 'exit_code': 'Exit Code', 'stderr': 'Standard Error', 'stdout': 'Standard Output'} %}
{% for status, desc in {'error': 'Errored', 'failure': 'Failed', 'success': 'Passed'}.items() if state[status]%}
<details><summary>{{ desc }} {{ execution_type }}s</summary>
{%   for execution in raw_data.executions %}
{%     if execution.data.status == status %}
{%       if execution.data.status == 'success' %}

* <details class="rcorners light-green"><summary class="light-green">&#9989; {{ execution.id }}</summary><div class="padded">

{%       else %}

* <details class="rcorners light-red"><summary class="light-red">&#10060; {{ execution.id }}</summary><div class="padded">

{%       endif %}
{%       if execution.data.output_problems %}
    **Problems**:
{%       endif %}
{%       for problem in execution.data.output_problems %}
    * ```
      {{problem|indent(6)}}
      ```
{%       endfor %}
{%       if execution.data.execution_problem %}
    **Execution Problem:**
    * ```
      {{execution.data.execution_problem|indent(6)}}
      ```
{%       endif %}
{%     if execution.data.job %}
{%       for key, description in display_job_attributes.items() %}
{%         if execution.data.job[key] not in ("", None) %}
    **{{ description }}:**

    * ```console
      {{ execution.data.job[key]|string|indent(6)}}
      ```
{%         endif %}
{%       endfor %}
{{render_job_parameters(execution.data.job)}}

{%     endif %}
{%     if execution.data.invocation_details %}

    #### Workflow invocation details

{{render_steps(execution.data.invocation_details.steps.values(), display_job_attributes)}}

{{render_invocation_details(execution.data.invocation_details.details)}}


  </div></details>

{% endif %}
{% endif %}
{% endfor %}

</details>
{% endfor %}
