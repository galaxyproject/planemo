<?xml version="1.0" encoding="UTF-8"?>
<testsuite name="{{ raw_data.suitename }}"
           tests="{{ raw_data.results.total }}"
           errors="{{ raw_data.results.errors }}"
           failures="{{ raw_data.results.failures }}"
           skip="{{ raw_data.results.skips }}">
    {% for testcase in raw_data.tests %}
    <testcase classname="{{ testcase.id }}" name="{{ testcase.data.test_index }}" time="{{ testcase.data.time_seconds }}">
        {% if testcase.data.status != 'success' %}
            <error type="error" message="Tool exit code: {{ testcase.data.job.exit_code }}"><![CDATA[
                {{ testcase.data | tojson(indent=True) }}
            ]]></error>
        {% endif %}
        <system-out><![CDATA[
        {{ testcase.data.job.stdout }}
        ]]></system-out>
        <system-err><![CDATA[
        {{ testcase.data.job.stderr }}
        ]]></system-err>
    </testcase>
    {% endfor %}
</testsuite>
