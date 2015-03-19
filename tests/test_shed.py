"""
"""

import os
import string
import random

from bioblend import toolshed

from .test_utils import (
    skip_unless_environ,
    CliTestCase,
)

from planemo.shed import username

SHED_TEMPLATE = string.Template("""owner: ${owner}
name: ${name}
description: Planemo test repository.
homepage_url: https://planemo.readthedocs.org/
remote_repository_url: https://github.com/galaxyproject/planemo/
type: tool_dependency_definition
categories:
  - cooltools
""")


class ShedTestCase(CliTestCase):

    @skip_unless_environ("TEST_TOOL_SHED_API_KEY")
    def test_shed(self):
        shed_url = os.environ.get("TEST_TOOL_SHED_URL",
                                  "http://localhost:9009")
        shed_api_key = os.environ.get("TEST_TOOL_SHED_API_KEY")
        tsi = toolshed.ToolShedInstance(shed_url, key=shed_api_key)
        owner = username(tsi)
        name = "planemotestrepo%d" % random.randint(0, 1000000)
        with self._isolate():
            shed_yml_contents = SHED_TEMPLATE.safe_substitute(
                owner=owner,
                name=name,
            )
            open(".shed.yml", "w").write(shed_yml_contents)
            init_cmd = [
                "shed_create",
                "--shed_key", shed_api_key,
                "--shed_target", shed_url
            ]
            self._check_exit_code(init_cmd)
