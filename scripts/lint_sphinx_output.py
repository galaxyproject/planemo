from __future__ import print_function

import re
import sys

IGNORE_PATTERNS = [
    re.compile(r'nonlocal image URI found'),
    re.compile(r'included in any toctree'),
]


def warning_line(line):
    if 'WARNING' not in line:
        return False
    if 'docs/tests' in line:  # Doesn't actually show up in docs so don't lint.
        return False
    for pat in IGNORE_PATTERNS:
        if pat.search(line):
            return False
    return True


def main(argv=None):
    if argv is None:
        argv = sys.argv
    sphinx_output = sys.stdin.read()
    warning_lines = [_ for _ in sphinx_output.splitlines() if warning_line(_)]
    for line in warning_lines:
        print(line)
    sys.exit(1 if warning_lines else 0)


if __name__ == "__main__":
    main()
