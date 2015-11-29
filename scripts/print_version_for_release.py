from __future__ import print_function
from distutils.version import LooseVersion
import ast
import re

_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('planemo/__init__.py', 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read().decode('utf-8')).group(1)))

# Strip .devN
version_tuple = LooseVersion(version).version[0:3]
print(".".join(map(str, version_tuple)))
