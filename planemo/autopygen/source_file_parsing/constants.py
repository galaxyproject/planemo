# stdlib_module_names is only available from python 3.10
try:
    from sys import stdlib_module_names
except ImportError:
    from stdlib_list import stdlib_list
    stdlib_module_names = stdlib_list()
WARNING_STRING = "##!_FIXME_!##"

STD_LIB_MODULE_NAMES = frozenset(stdlib_module_names)
