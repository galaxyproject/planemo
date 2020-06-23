from setuptools import setup

setup(
    name='cmd_autoupdate',
    version='0.1',
    py_modules=['cmd_autoupdate'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        cmd_autoupdate=cmd_autoupdate:cli
    ''',
)