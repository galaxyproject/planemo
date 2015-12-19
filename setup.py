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

PROJECT_NAME = "planemo"
SOURCE_DIR = "planemo"
TEST_DIR = 'tests'
PROJECT_DESCRIPTION = 'Command-line utilities to assist in building tools for the Galaxy project (http://galaxyproject.org/).'
PROJECT_AUTHOR = 'Galaxy Project and Community'
PROJECT_EMAIL = 'jmchilton@gmail.com'
PROJECT_URL = 'https://github.com/galaxyproject/planemo'
PACKAGES = [
    'planemo',
    'planemo.cwl',
    'planemo.commands',
    'planemo.galaxy_test',
    'planemo.linters',
    'planemo.reports',
    'planemo.shed',
    'planemo.shed2tap',
    'planemo.xml',
    'planemo_ext',
    'planemo_ext.cwl2script',
    'planemo_ext.galaxy',
    'planemo_ext.galaxy.eggs',
    'planemo_ext.galaxy.tools',
    'planemo_ext.galaxy.tools.linters',
    'planemo_ext.galaxy.tools.deps',
    'planemo_ext.galaxy.tools.deps.resolvers',
    'planemo_ext.galaxy.util',
]
ENTRY_POINTS = '''
    [console_scripts]
    planemo=planemo.cli:planemo
'''
PACKAGE_DATA = {
    'planemo_ext': [
        'galaxy/util/docutils_template.txt',
        'tool_factory_2/rgToolFactory2.xml',
        'tool_factory_2/rgToolFactory2.py',
        'tool_factory_2/getlocalrpackages.py',
    ],
    'planemo': [
        'xml/xsd/repository_dependencies.xsd',
        'xml/xsd/tool_dependencies.xsd',
        'xml/xsd/tool/galaxy.xsd',
        'xml/xsd/tool/citation.xsd',
        'xml/xsd/tool/citations.xsd',
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
    requirements = open("requirements.txt").read().split("\n")
else:
    # In tox, it will cover them anyway.
    requirements = []

# Only import cwltool for Python 2.7.
if sys.version_info[0] == 2 and sys.version_info[1] >= 7:
    requirements.append("cwltool")


test_requirements = [
    # TODO: put package test requirements here
]

_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('%s/__init__.py' % SOURCE_DIR, 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read().decode('utf-8')).group(1)))

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
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite=TEST_DIR,
    tests_require=test_requirements
)
