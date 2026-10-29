"""Fail if HISTORY.rst's current .devN section has no bullet entries."""

import re
import sys
from pathlib import Path

HISTORY = Path(__file__).resolve().parent.parent / "HISTORY.rst"


def main() -> int:
    text = HISTORY.read_text()
    m = re.search(r"-{3,}\n(\d+\.\d+\.\d+\.dev\d+)\n-{3,}\n(.*?)(?=\n-{3,}\n|\Z)", text, re.DOTALL)
    if not m:
        print(f"ERROR: no .devN section found in {HISTORY}", file=sys.stderr)
        return 1
    version, body = m.group(1), m.group(2)
    bullets = [line for line in body.splitlines() if line.lstrip().startswith("*")]
    if not bullets:
        print(f"ERROR: {version} section in HISTORY.rst has no entries; run 'make add-history'", file=sys.stderr)
        return 1
    print(f"OK: {version} has {len(bullets)} entr{'y' if len(bullets) == 1 else 'ies'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
