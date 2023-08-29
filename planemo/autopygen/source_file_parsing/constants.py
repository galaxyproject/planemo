# stdlib_module_names is only available from python 3.10
# the inner try-except is essentially needed for mypy
# since we assume that stdlib-list is available
try:
    from sys import stdlib_module_names
except ImportError:
    try:
        from stdlib_list import stdlib_list
        stdlib_module_names = stdlib_list()
    except ImportError:
        # Handle the case where stdlib_list is not available
        stdlib_module_names = []

WARNING_STRING = "##!_FIXME_!##"

STD_LIB_MODULE_NAMES = frozenset(stdlib_module_names)
