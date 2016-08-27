"""Check an output file from a generalize artifact test."""

import os

from galaxy.tools.verify import verify


def check_output(runnable, output_properties, test_properties, **kwds):
    """Use galaxy-lib to check a test output.

    Return a list of strings describing the problems encountered,
    and empty list indicates no problems were detected.

    Currently this will only ever return at most one detected problem because
    of the way galaxy-lib throws exceptions instead of returning individual
    descriptions - but this may be enhanced in the future.
    """
    get_filename = _test_filename_getter(runnable)
    path = output_properties["path"]
    output_content = open(path, "rb").read()
    expected_file = test_properties.get("file", None)
    job_output_files = kwds.get("job_output_files", None)
    item_label = "Output with path %s" % path
    problems = []
    try:
        verify(
            item_label,
            output_content,
            attributes=test_properties,
            filename=expected_file,
            get_filename=get_filename,
            keep_outputs_dir=job_output_files,
            verify_extra_files=None,
        )
    except AssertionError as e:
        problems.append(str(e))

    return problems


def _test_filename_getter(runnable):

    def get_filename(name):
        artifact_directory = os.path.dirname(runnable.path)
        return os.path.join(artifact_directory, name)

    return get_filename


__all__ = [
    "check_output",
]
