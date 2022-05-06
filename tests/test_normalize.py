""" Integration tests for normalize command.
"""

import os

from .test_utils import (
    CliTestCase,
    TEST_TOOLS_DIR,
)


class NormalizeTestCase(CliTestCase):
    def test_shed_diff(self):
        with self._isolate() as f:
            fail_order = os.path.join(TEST_TOOLS_DIR, "fail_order.xml")
            # Doesn't pass linting before normalization...
            self._check_exit_code(["lint", fail_order], exit_code=1)
            result = self._check_exit_code(["normalize", fail_order])
            good_path = os.path.join(f, "good.xml")
            with open(good_path, "w") as good_f:
                good_f.write(result.output)
            # Does after!
            self._check_exit_code(["lint", good_path])
