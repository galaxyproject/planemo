from setuptools import setup


def read_requirements():
    with open("requirements.txt") as req:
        return req.read().splitlines()


setup(
    install_requires=read_requirements(),
)
