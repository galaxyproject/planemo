"""Test utilities from :module:`planemo.io`."""
from .test_utils import (
    io,
    assert_equal,
)


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
