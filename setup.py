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

if sys.version_info < (2, 7):
    sys.stderr.write("ERROR: planemo requires at least Python Version 2.7\n")
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
    'planemo.bioconda_scripts',
    'planemo.cwl',
    'planemo.cwl.cwl2script',
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
    'planemo.xml',
]
ENTRY_POINTS = '''
    [console_scripts]
    planemo=planemo.cli:planemo
'''
PACKAGE_DATA = {
    'planemo_ext': [
        'tool_factory_2/rgToolFactory2.xml',
        'tool_factory_2/rgToolFactory2.py',
        'tool_factory_2/getlocalrpackages.py',
    ],
    'planemo': [
        'xml/xsd/repository_dependencies.xsd',
        'xml/xsd/tool_dependencies.xsd',
        'xml/xsd/tool/galaxy.xsd',
        'reports/*',
        'scripts/*',
    ],
}
PACKAGE_DIR = {
    SOURCE_DIR: SOURCE_DIR,
    'planemo_ext': 'planemo_ext'
}

readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

if os.path.exists("requirements.txt"):
    requirements = [ r for r in open("requirements.txt").read().split("\n") if ";" not in r ]
    py27_requirements = [ r.split(";", 1)[0].strip() for r in open("requirements.txt").read().split("\n") if ";" in r ]
    if not PLANEMO_REQUIRE_LXML:
        requirements.remove("lxml")
else:
    # In tox, it will cover them anyway.
    requirements = []
    py27_requirements = []

test_requirements = [
    # TODO: put package test requirements here
]


setup(
    name=PROJECT_NAME,
    version=version,
    description=PROJECT_DESCRIPTION,
    long_description=readme + '\n\n' + history,
    author=PROJECT_AUTHOR,
    author_email=PROJECT_EMAIL,
    url=PROJECT_URL,
    packages=PACKAGES,
    entry_points=ENTRY_POINTS,
    package_data=PACKAGE_DATA,
    package_dir=PACKAGE_DIR,
    include_package_data=True,
    install_requires=requirements,
    extras_require={
        ':python_version=="2.7"': py27_requirements,
    },
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
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
    ],
    test_suite=TEST_DIR,
    tests_require=test_requirements
)
