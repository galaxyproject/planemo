<?xml version="1.0" encoding="UTF-8"?>
<testsuite name="planemo-{{ suitename }}"
           tests="{{ results.total }}"
           errors="{{ results.errors }}"
           failures="{{ results.failures }}"
           skip="{{ results.skips }}">
    {% for testcase in tests %}
    <testcase classname="{{ testcase.classname }}" name="{{ testcase.name }}" time="{{ testcase.time }}">
        {% if 'errorType' in testcase %}
            <error type="planemo.{{ testcase.errorType }}" message="{{ testcase.errorMessage }}">
            {{ testcase.errorContent }}
            </error>
        {% endif %}
    </testcase>
    {% endfor %}
</testsuite>
