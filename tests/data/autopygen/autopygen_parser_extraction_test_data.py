def import_discovery():
    import sys
    import os
    import shoud_not_be_imported as snbi
    from os import path as asd
    import argparse as parser


def parser_discovery_and_replacement():
    import argparse

    parser = argparse.ArgumentParser()


def group_discovery():
    import argparse

    parser = argparse.ArgumentParser()
    group = parser.add_argument_group("test")
    group.add_argument()


def argument_creation():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--test")
