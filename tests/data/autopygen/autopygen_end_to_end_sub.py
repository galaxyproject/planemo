# /usr/bin/python3
import argparse

parent = argparse.ArgumentParser(description='Process some integers.', prefix_chars='-+', epilog="here's some epilog text", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


subparsers = parent.add_subparsers()

parser_foo = subparsers.add_parser('foo')
parser_bar = subparsers.add_parser('bar')

parser_foo.add_argument('keyword', metavar='Q', type=str, nargs=1, help='action keyword')

parser_foo.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')

parser_foo.add_argument('--sum', '-s', dest='accumulate', action='store_const',
                    const=sum, default=max, help='sum the integers (default: find the max)')

parser_foo.add_argument('--foo', nargs='?', help='foo help')
parser_foo.add_argument('--bar', nargs='*', default=[1, 2, 3], help='BAR!')
parser_foo.add_argument('--true', action='store_true', help='Store a true')
parser_bar.add_argument('--false', action='store_false', help='Store a false')
parser_bar.add_argument('--append', action='append', help='Append a value')

parser_bar.add_argument('--nargs2', nargs=2, help='nargs2')

parser_bar.add_argument('--mode', choices=['rock', 'paper', 'scissors'], default='scissors')


parser_bar.add_argument('--version', action='version', version='2.0')
args = parent.parse_args()
