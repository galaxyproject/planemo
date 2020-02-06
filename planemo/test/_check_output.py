"""Check an output file from a generalize artifact test."""

import os

from galaxy.tool_util.verify import verify
from galaxy.util import unicodify


def check_output(runnable, output_properties, test_properties, **kwds):
    """Use galaxy-tool-util to check a test output.

    Return a list of strings describing the problems encountered,
    and empty list indicates no problems were detected.

    Currently this will only ever return at most one detected problem because
    of the way galaxy-tool-util throws exceptions instead of returning individual
    descriptions - but this may be enhanced in the future.
    """
    if "elements" in test_properties:
        checker = _check_output_collection
    else:
        checker = _check_output_file
    return checker(runnable, output_properties, test_properties, **kwds)


def _check_output_collection(runnable, output_properties, test_properties, **kwds):
    data_collection = self._get("dataset_collections/%s" % output_collection_id, data={"instance_type": "history"}).json()

    def verify_dataset(element, element_attrib, element_outfile):
        pass

    verify_collection(output_collection_def, data_collection, verify_dataset)

def _check_output_file(runnable, output_properties, test_properties, **kwds):
    get_filename = _test_filename_getter(runnable)
    path = output_properties["path"]
    with open(path, "rb") as fh:
        output_content = fh.read()
    # Support Galaxy-like file location (using "file") or CWL-like ("path" or "location").
    expected_file = test_properties.get("file", None)
    if expected_file is None:
        expected_file = test_properties.get("path", None)
    if expected_file is None:
        expected_file = test_properties.get("location", None)

    job_output_files = kwds.get("job_output_files", None)
    item_label = "Output with path %s" % path
    problems = []
    if "asserts" in test_properties:
        # TODO: break fewer abstractions here...
        from galaxy.tool_util.parser.yaml import __to_test_assert_list
        test_properties["assert_list"] = __to_test_assert_list(test_properties["asserts"])
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
        problems.append(unicodify(e))

    return problems


def _test_filename_getter(runnable):

    def get_filename(name):
        artifact_directory = os.path.dirname(runnable.path)
        return os.path.join(artifact_directory, name)

    return get_filename


__all__ = (
    "check_output",
)
