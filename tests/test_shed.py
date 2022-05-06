"""
"""

import os
import random
import string

import yaml
from bioblend import toolshed

from planemo import io
from planemo.shed import username
from .test_utils import (
    CliTestCase,
    skip_unless_environ,
    TEST_DIR,
)

SHED_TEMPLATE = string.Template(
    """owner: ${owner}
name: ${name}
description: Planemo test repository.
homepage_url: https://planemo.readthedocs.org/
remote_repository_url: https://github.com/galaxyproject/planemo/
type: tool_dependency_definition
categories:
  - cooltools
"""
)


class ShedTestCase(CliTestCase):
    @skip_unless_environ("TEST_TOOL_SHED_API_KEY")
    def test_shed(self):
        shed_url = os.environ.get("TEST_TOOL_SHED_URL", "http://localhost:9009")
        shed_api_key = os.environ.get("TEST_TOOL_SHED_API_KEY")
        tsi = toolshed.ToolShedInstance(shed_url, key=shed_api_key)
        owner = username(tsi)
        name = "planemotestrepo%d" % random.randint(0, 1000000)
        with self._isolate():
            shed_yml_contents = SHED_TEMPLATE.safe_substitute(
                owner=owner,
                name=name,
            )
            io.write_file(".shed.yml", shed_yml_contents)
            test_path = os.path.join(TEST_DIR, "tool_dependencies_good_1.xml")
            with open(test_path) as fh:
                contents = fh.read()
            io.write_file("tool_dependencies.xml", contents)
            init_cmd = ["shed_create", "--shed_key", shed_api_key, "--shed_target", shed_url]
            self._check_exit_code(init_cmd)
            with open(".shed.yml") as fh:
                contents_dict = yaml.safe_load(fh)
            contents_dict["description"] = "Update test repository."
            io.write_file(".shed.yml", yaml.dump(contents_dict))
            update_cmd = ["shed_update", "--shed_key", shed_api_key, "--shed_target", shed_url]
            self._check_exit_code(update_cmd)
            upload_cmd = ["shed_upload", "--shed_key", shed_api_key, "--shed_target", shed_url]
            self._check_exit_code(upload_cmd)
            download_cmd = ["shed_download", "--shed_target", shed_url, "--destination", "shed_download.tar.gz"]
            self._check_exit_code(download_cmd)
