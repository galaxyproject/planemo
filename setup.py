#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ast
import os
import re
import sys
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.version_info < (3, 5):
    sys.stderr.write("ERROR: planemo requires at least Python Version 3.5\n")
    sys.exit(1)

# Allow installer to turn off dependency on lxml by setting the environment variable
# PLANEMO_REQUIRE_LXML to "0". lxml should be considered optional if xmllint is
# available on the PATH - but python doesn't really provide me a fantastic way to
# express that.
DEFAULT_PLANEMO_REQUIRE_LXML = 1
PLANEMO_REQUIRE_LXML = os.environ.get("PLANEMO_REQUIRE_LXML", "%d" % DEFAULT_PLANEMO_REQUIRE_LXML) != "0"

SOURCE_DIR = "planemo"

_version_re = re.compile(r'__version__\s+=\s+(.*)')


with open('%s/__init__.py' % SOURCE_DIR, 'rb') as f:
    init_contents = f.read().decode('utf-8')

    def get_var(var_name):
        pattern = re.compile(r'%s\s+=\s+(.*)' % var_name)
        match = pattern.search(init_contents).group(1)
        return str(ast.literal_eval(match))

    version = get_var("__version__")
    PROJECT_NAME = get_var("PROJECT_NAME")
    PROJECT_URL = get_var("PROJECT_URL")
    PROJECT_AUTHOR = get_var("PROJECT_AUTHOR")
    PROJECT_EMAIL = get_var("PROJECT_EMAIL")

TEST_DIR = 'tests'
PROJECT_DESCRIPTION = 'Command-line utilities to assist in building tools for the Galaxy project (http://galaxyproject.org/).'
PACKAGES = [
    'planemo',
    'planemo.cwl',
    'planemo.commands',
    'planemo.conda_verify',
    'planemo.database',
    'planemo.engine',
    'planemo.galaxy',
    'planemo.galaxy.test',
    'planemo.linters',
    'planemo.reports',
    'planemo.shed',
    'planemo.shed2tap',
    'planemo.test',
    'planemo.training',
    'planemo.xml',
]
ENTRY_POINTS = '''
    [console_scripts]
    planemo=planemo.cli:planemo
'''
PACKAGE_DATA = {
    'planemo': [
        'xml/xsd/repository_dependencies.xsd',
        'xml/xsd/tool_dependencies.xsd',
        'reports/*',
        'scripts/*',
    ],
}
PACKAGE_DIR = {
    SOURCE_DIR: SOURCE_DIR,
}

with open('README.rst') as fh:
    readme = fh.read()
with open('HISTORY.rst') as fh:
    history = fh.read().replace('.. :changelog:', '')

if os.path.exists("requirements.txt"):
    with open("requirements.txt") as fh:
        requirements = [r for r in fh.read().split("\n") if ";" not in r]
    if not PLANEMO_REQUIRE_LXML:
        requirements.remove("lxml")
else:
    # In tox, it will cover them anyway.
    requirements = []

test_requirements = [
    # TODO: put package test requirements here
]


setup(
    name=PROJECT_NAME,
    version=version,
    description=PROJECT_DESCRIPTION,
    long_description=readme + '\n\n' + history,
    long_description_content_type="text/x-rst",
    author=PROJECT_AUTHOR,
    author_email=PROJECT_EMAIL,
    url=PROJECT_URL,
    packages=PACKAGES,
    entry_points=ENTRY_POINTS,
    package_data=PACKAGE_DATA,
    package_dir=PACKAGE_DIR,
    include_package_data=True,
    install_requires=requirements,
    license="AFL",
    zip_safe=False,
    keywords='planemo',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Environment :: Console',
        'License :: OSI Approved :: Academic Free License (AFL)',
        'Operating System :: POSIX',
        'Topic :: Software Development',
        'Topic :: Software Development :: Code Generators',
        'Topic :: Software Development :: Testing',
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    test_suite=TEST_DIR,
    tests_require=test_requirements
)
