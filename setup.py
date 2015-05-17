#!/usr/bin/env python
# -*- coding: utf-8 -*-


import ast
import re
import sys
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
    'glob2',
]

# Latest stable bioblend does not support Python 3, setup dev dependency.
dependency_links = []
if sys.version_info[0] != 2:
    dependency_links = ['git://github.com/galaxyproject/bioblend#egg=bioblend']


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
        'planemo.galaxy_test',
        'planemo.linters',
        'planemo.reports',
        'planemo.shed',
        'planemo.shed2tap',
        'planemo.xml',
        'planemo_ext',
        'planemo_ext.galaxy',
        'planemo_ext.galaxy.eggs',
        'planemo_ext.galaxy.tools',
        'planemo_ext.galaxy.tools.linters',
        'planemo_ext.galaxy.tools.deps',
        'planemo_ext.galaxy.tools.deps.resolvers',
        'planemo_ext.galaxy.util',
    ],
    entry_points='''
        [console_scripts]
        planemo=planemo.cli:planemo
    ''',
    package_data={'planemo_ext': ['galaxy/util/docutils_template.txt',
                                  'tool_factory_2/rgToolFactory2.xml',
                                  'tool_factory_2/rgToolFactory2.py',
                                  'tool_factory_2/getlocalrpackages.py',
                                 ],
                  'planemo': ['xml/xsd/repository_dependencies.xsd',
                              'xml/xsd/tool_dependencies.xsd',
                              'xml/xsd/tool/galaxy.xsd',
                              'xml/xsd/tool/citation.xsd',
                              'xml/xsd/tool/citations.xsd',
                              'reports/*',
                              'scripts/*',
                             ],
                 },
    package_dir={'planemo': 'planemo',
                 'planemo_ext': 'planemo_ext'},
    include_package_data=True,
    install_requires=requirements,
    dependency_links=dependency_links,
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
