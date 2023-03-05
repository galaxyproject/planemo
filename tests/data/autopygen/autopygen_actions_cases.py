def none_positional():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("test")


def none():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test")


def store():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store")


def store_const():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_const", const=42)


def store_const_text():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_const", const="asd")


def store_true():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_true")


def append():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="append")


def append_const():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="append_const", const=444)


def count():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="count", default=0)


def version():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action="version", version="2.0")


def extend():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="extend", nargs="+", type=str)
