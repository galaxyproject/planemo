#!/usr/bin/env python
# -*- coding: utf-8 -*-


import re
import ast
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = [
    'Click',
    'six',
    'pyyaml',
    'jinja2',
    'docutils',
    'PyGithub',
    'bioblend',
]

test_requirements = [
    # TODO: put package test requirements here
]

_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('planemo/__init__.py', 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read().decode('utf-8')).group(1)))

setup(
    name='planemo',
    version=version,
    description='Command-line utilities to assist in building tools for the Galaxy project (http://galaxyproject.org/).',
    long_description=readme + '\n\n' + history,
    author='Galaxy Project and Community',
    author_email='jmchilton@gmail.com',
    url='https://github.com/galaxyproject/planemo',
    packages=[
        'planemo',
        'planemo.commands',
        'planemo.reports',
        'planemo_ext',
        'planemo_ext.galaxy',
        'planemo_ext.galaxy.eggs',
        'planemo_ext.galaxy.tools',
        'planemo_ext.galaxy.tools.linters',
        'planemo_ext.galaxy.tools.deps',
        'planemo_ext.galaxy.tools.deps.resolvers',
        'planemo_ext.galaxy.util',
    ],
    data_files=[('tool_factory_2', ['tool_factory_2/rgToolFactory2.xml',
                                    'tool_factory_2/rgToolFactory2.py',
                                    'tool_factory_2/getlocalrpackages.py']),
                ('planemo_ext/galaxy/util', ['planemo_ext/galaxy/util/docutils_template.txt']),
               ],
    entry_points='''
        [console_scripts]
        planemo=planemo.cli:planemo
    ''',
    package_dir={'planemo': 'planemo',
                 'planemo_ext': 'planemo_ext'},
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
    test_suite='tests',
    tests_require=test_requirements
)
