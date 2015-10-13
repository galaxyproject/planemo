<?xml version="1.0" encoding="UTF-8"?>
<testsuite name="planemoupdate"
           tests="{{ results.total }}"
           errors="{{ results.errors }}"
           failures="{{ results.failures }}"
           skip="{{ results.skips }}">
    {% for testcase in tests %}
    <testcase classname="{{ testcase.classname }}" name="recursive_update" time="{{ testcase.time }}">
        {% if testcase.result == 1 %}
            <error type="planemo.FailedUpdate" message="Failed to update repository">
            </error>
        {% endif %}
        {% if testcase.result == 2 %}
            <error type="planemo.RepoDoesNotExit" message="Target repository does not exist">
            </error>
        {% endif %}
        {% if testcase.result == 3 %}
            <skipped type="planemo.FailedMetadata" message="Repo updated but metadata was not">
            </skipped>
        {% endif %}
        {% if testcase.result > 3 %}
            <error type="planemo.UpdateError" message="Error executing update">
            </error>
        {% endif %}
    </testcase>
    {% endfor %}
</testsuite>
