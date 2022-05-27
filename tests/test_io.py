"""Test utilities from :module:`planemo.io`."""
import tempfile

from planemo import io
from .test_utils import assert_equal


def test_io_capture():
    """Test :func:`planemo.io.conditionally_captured_io`."""
    with io.conditionally_captured_io(True, tee=False) as capture:
        io.warn("Problem...")
    assert_equal(capture[0]["data"], "Problem...")

    with io.conditionally_captured_io(True, tee=False) as capture:
        io.shell("echo 'Problem...'")
    assert_equal(capture[0]["data"], "echo 'Problem...'")
    assert_equal(capture[1]["data"], "Problem...")

    with io.conditionally_captured_io(True, tee=False) as capture:
        io.communicate("echo 'Problem...'")
    assert_equal(capture[0]["data"], "echo 'Problem...'")
    assert_equal(capture[1]["data"], "Problem...")

    with io.conditionally_captured_io(False, tee=False) as capture:
        io.communicate("echo 'Test...'")

    assert capture is None


def test_filter_paths():
    """Test :func:`planemo.io.filter_paths`."""
    test_cwd = "/a/b"

    def assert_filtered_is(paths, expected, **kwds):
        result = io.filter_paths(paths, cwd=test_cwd, **kwds)
        assert result == expected, f"paths [{result}] arent't expected [{expected}]"

    assert_filtered_is([], [], exclude=["/a"])
    assert_filtered_is(["/a/c"], [], exclude=["/a"])
    assert_filtered_is(["/b"], ["/b"], exclude=["/a"])
    assert_filtered_is(["/a/b/c"], [], exclude=["c"])
    with tempfile.NamedTemporaryFile(mode="w+") as tmp:
        tmp.write("#exclude c\n\nc\n")
        tmp.flush()
        assert_filtered_is(["/a/b/c", "/a/b/d"], ["/a/b/d"], exclude_from=[tmp.name])
