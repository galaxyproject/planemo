import collections
import sys

from planemo.conda_verify.const import (
    DLL_TYPES,
    MAGIC_HEADERS,
)


def get_object_type(data):
    head = data[:4]
    if head not in MAGIC_HEADERS:
        return None
    lookup = MAGIC_HEADERS.get(head)
    if lookup == 'DLL':
        pos = data.find('PE\0\0')
        if pos < 0:
            return "<no PE header found>"
        i = ord(data[pos + 4]) + 256 * ord(data[pos + 5])
        return "DLL " + DLL_TYPES.get(i)
    elif lookup.startswith('MachO'):
        return lookup
    elif lookup == 'ELF':
        return "ELF" + {'\x01': '32', '\x02': '64'}.get(data[4])


def get_bad_seq(s):
    for seq in ('--', '-.', '-_',
                '.-', '..', '._',
                '_-', '_.'):  # but '__' is fine
        if seq in s:
            return seq
    return None


def all_ascii(data, allow_CR=False):
    newline = [10]  # LF
    if allow_CR:
        newline.append(13)  # CF
    for c in data:
        n = ord(c) if sys.version_info[0] == 2 else c
        if not (n in newline or 32 <= n < 127):
            return False
    return True


class memoized(object):
    """Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value


if __name__ == '__main__':
    print(sys.version)
    print(all_ascii(b'Hello\x00'), all_ascii(b"Hello World!"))
