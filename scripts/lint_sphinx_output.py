from __future__ import print_function
import sys
import re

IGNORE_PATTERNS = [
    re.compile(r'nonlocal image URI found'),
    re.compile(r'included in any toctree'),
]


def warning_line(line):
    if 'WARNING' not in line:
        return False
    for pat in IGNORE_PATTERNS:
        if pat.search(line):
            return False
    return True


def main(argv=None):
    if argv is None:
        argv = sys.argv
    sphinx_output = sys.stdin.read()
    warning_lines = filter(warning_line, sphinx_output.splitlines())
    for line in warning_lines:
        print(line)
    sys.exit(1 if warning_lines else 0)


if __name__ == "__main__":
    main()
