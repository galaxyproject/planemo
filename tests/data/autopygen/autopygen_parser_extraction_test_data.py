def import_discovery():
    import sys  # noqa: F401, F841

    import os  # noqa: F401, F841

    import shoud_not_be_imported as snbi  # noqa: F401, F841

    from os import path as asd  # noqa: F401, F841

    import argparse as parser  # noqa: F401, F841


def parser_discovery_and_replacement():
    import argparse
    parser = argparse.ArgumentParser()  # noqa: F401, F841


def group_discovery():
    import argparse

    parser = argparse.ArgumentParser()
    group = parser.add_argument_group("test")
    group.add_argument()


def argument_creation():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--test")
