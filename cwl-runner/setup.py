#!/usr/bin/env python

import os

from setuptools import setup

SETUP_DIR = os.path.dirname(__file__)
README = os.path.join(SETUP_DIR, "README")

scripts = ["cwl-runner"]

with open(README) as fh:
    readme_content = fh.read()

setup(
    name="cwl_runner",
    version="1.0",
    description="Common workflow language interpreter implementation (for Galaxy + Planemo)",
    long_description=readme_content,
    author="John Chilton",
    author_email="jmchilton@gmail.com",
    url="https://github.com/galaxyproject/planemo",
    download_url="https://github.com/galaxyproject/planemo",
    license="MIT",
    install_requires=["planemo"],
    scripts=scripts,
    zip_safe=True,
)
