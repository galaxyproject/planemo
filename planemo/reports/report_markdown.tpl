# {{ title }}

## Executive Summary

| Test State | Count |
| ---------- | ----- |
| Total      | {{ raw_data.summary.num_tests | default(0)  }} |
| Passed     | {{ raw_data.summary.num_tests - raw_data.summary.num_errors - raw_data.summary.num_failures - raw_data.summary.num_skips | default(0)  }} |
| Error      | {{ raw_data.summary.num_errror | default(0) }} |
| Failure    | {{ raw_data.summary.num_failure | default(0) }} |
| Skipped    | {{ raw_data.summary.num_skipped | default(0) }} |


## Detailed Results
{% for test in raw_data.tests %}
### {{ test.id }}
{% if test.data.status == 'success' %}
Job Passed
{% else %}
Job Error! (State: {{ test.data.status }})

Command Line:

```console
{{ test.data.job.command_line}}
```

exited with code {{ test.data.job.exit_code }}.

#### Problems

{% for problem in test.data.output_problems %}
```console
{{problem}}
```
{% endfor %}

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

{% endif %}
{% endif %}
{% endfor %}
