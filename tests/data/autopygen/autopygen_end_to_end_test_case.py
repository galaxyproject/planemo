# borrowed with permission from https://github.com/hexylena/argparse2tool/blob/main/examples/example.py
# /usr/bin/python3
import argparse

parser = argparse.ArgumentParser(
    description="Process some integers.",
    prefix_chars="-+",
    epilog="here's some epilog text",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("keyword", metavar="Q", type=str, nargs=1, help="action keyword")

parser.add_argument("integers", metavar="N", type=int, nargs="+", help="an integer for the accumulator")

parser.add_argument(
    "--sum",
    "-s",
    dest="accumulate",
    action="store_const",
    const=sum,
    default=max,
    help="sum the integers (default: find the max)",
)

parser.add_argument("--foo", nargs="?", help="foo help")
parser.add_argument("--bar", nargs="*", default=[1, 2, 3], help="BAR!")
parser.add_argument("--true", action="store_true", help="Store a true")
parser.add_argument("--false", action="store_false", help="Store a false")
parser.add_argument("--append", action="append", help="Append a value")

parser.add_argument("--nargs2", nargs=2, help="nargs2")

parser.add_argument("--mode", choices=["rock", "paper", "scissors"], default="scissors")


parser.add_argument("--version", action="version", version="2.0")
args = parser.parse_args()
