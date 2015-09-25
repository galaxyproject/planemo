<?xml version="1.0" encoding="UTF-8"?>
<testsuite name="planemodiff"
           tests="{{ results.total }}"
           errors="{{ results.errors }}"
           failures="{{ results.failures }}"
           skip="{{ results.skips }}">
    {% for testcase in tests %}
    <testcase classname="{{ testcase.classname }}" name="recursive_diff" time="{{ testcase.time }}">
        {% if testcase.result == 1 %}
            <error type="planemo.Different" message="Repository is different">
                <!-- TODO: copy of diff output -->
            </error>
        {% endif %}
        {% if testcase.result == 2 %}
            <error type="planemo.RepoDoesNotExit" message="Target repository does not exist">
                <!-- TODO: copy of diff output -->
            </error>
        {% endif %}
        {% if testcase.result > 2 %}
            <error type="planemo.DiffError" message="Error executing diff">
                <!-- TODO: copy of diff output -->
            </error>
        {% endif %}
    </testcase>
    {% endfor %}
</testsuite>
