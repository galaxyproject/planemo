import sys
import traceback
from planemo.io import error

from galaxy.tools import loader_directory


def load_tool_elements_from_path(path, recursive):
    return loader_directory.load_tool_elements_from_path(
        path,
        load_exception_handler,
        recursive,
    )


def load_exception_handler(path, exc_info):
    error("Error loading tool with path %s" % path)
    traceback.print_exception(*exc_info, limit=1, file=sys.stderr)
