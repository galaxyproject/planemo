import ast
import os.path
import re
import sys

from packaging.version import Version

source_dir = sys.argv[1]

_version_re = re.compile(r"__version__\s+=\s+(.*)")

with open(os.path.join(source_dir, "__init__.py")) as f:
    version_match = _version_re.search(f.read())
assert version_match
version = ast.literal_eval(version_match.group(1))

version_obj = Version(version)
# Strip .devN
print(version_obj.base_version)
