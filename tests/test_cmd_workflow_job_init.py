"""Tests for the ``workflow_job_init`` command."""

import os

import yaml

from .test_utils import CliTestCase


class CmdWorkflowJobInitTestCase(CliTestCase):
    def test_plain_init(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["workflow_job_init", "0_basic_native.yml"]
            self._check_exit_code(init_cmd)
            job_path = os.path.join(f, "0_basic_native_job.yml")
            assert os.path.exists(job_path)
            with open(job_path) as stream:
                job = yaml.safe_load(stream)
            assert isinstance(job, dict)
            assert "the_input" in job
            assert job.get("the_input").get("path") == "todo_test_data_path.ext"
            # Check that comments with type and doc are present
            with open(job_path) as stream:
                content = stream.read()
            assert "# type: data, doc: input doc" in content

    def test_cannot_overwrite(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["workflow_job_init", "0_basic_native.yml"]
            job_path = os.path.join(f, "0_basic_native_job.yml")
            with open(job_path, "w") as f:
                f.write("already exists")
            self._check_exit_code(init_cmd, exit_code=1)

    def test_force(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["workflow_job_init", "--force", "0_basic_native.yml"]
            job_path = os.path.join(f, "0_basic_native_job.yml")
            with open(job_path, "w") as f:
                f.write("already exists")
            self._check_exit_code(init_cmd, exit_code=0)
            with open(job_path) as stream:
                job = yaml.safe_load(stream)
            assert isinstance(job, dict)
            assert "the_input" in job

    def test_integer_parameter(self):
        """Test that integer parameters are handled correctly."""
        with self._isolate_with_test_data("") as f:
            init_cmd = ["workflow_job_init", "wf9-int-input.gxwf.yml"]
            self._check_exit_code(init_cmd)
            job_path = os.path.join(f, "wf9-int-input.gxwf-job.yml")
            assert os.path.exists(job_path)

            with open(job_path) as stream:
                content = stream.read()
            # Check integer parameter comment
            assert "# type: int" in content
            # Check file parameter comment
            assert "# type: data" in content

            with open(job_path) as stream:
                job = yaml.safe_load(stream)
            # Integer params should have placeholder value
            assert job.get("lines") == "todo_param_value"
            # File params should have class and path
            assert job.get("input1").get("class") == "File"

    def test_various_parameter_types(self):
        """Test workflow with various parameter types including optional, defaults, and format."""
        with self._isolate_with_test_data("") as f:
            init_cmd = ["workflow_job_init", "wf_param_types.gxwf.yml"]
            self._check_exit_code(init_cmd)
            job_path = os.path.join(f, "wf_param_types.gxwf_job.yml")
            assert os.path.exists(job_path)

            with open(job_path) as stream:
                content = stream.read()

            # Check required file input
            assert "# type: data, doc: A required file input" in content

            # Check optional file input - should have optional: true
            assert "optional: true" in content
            assert "doc: An optional file input" in content

            # Check format restriction
            assert "format: fastq" in content

            # Check integer with default - should show default value
            assert "default: 10" in content

            # Check boolean with default
            assert "type: boolean" in content
            assert "default: True" in content or "default: true" in content

            with open(job_path) as stream:
                job = yaml.safe_load(stream)

            # Integer with default should use the default value
            assert job.get("int_with_default") == 10

            # Boolean with default should use the default value
            assert job.get("bool_param") is True

            # Integer without default should have placeholder
            assert job.get("int_param") == "todo_param_value"

            # String without default should have placeholder
            assert job.get("string_param") == "todo_param_value"

            # All file inputs should be present
            assert "required_file" in job
            assert "optional_file" in job
            assert "fastq_input" in job

    def test_boolean_parameter(self):
        """Test that boolean parameters are handled correctly."""
        with self._isolate_with_test_data("") as f:
            init_cmd = ["workflow_job_init", "wf18_simple_conditional.yml"]
            self._check_exit_code(init_cmd)
            job_path = os.path.join(f, "wf18_simple_conditional_job.yml")
            assert os.path.exists(job_path)

            with open(job_path) as stream:
                content = stream.read()
            # Check boolean parameter comment
            assert "# type: boolean" in content

            with open(job_path) as stream:
                job = yaml.safe_load(stream)
            # Boolean without default should use false
            assert job.get("should_run") is False

    def test_collection_types(self):
        """Test that collection types generate appropriate sample entries."""
        with self._isolate_with_test_data("") as f:
            init_cmd = ["workflow_job_init", "wf_collection_types.gxwf.yml"]
            self._check_exit_code(init_cmd)
            job_path = os.path.join(f, "wf_collection_types.gxwf_job.yml")
            assert os.path.exists(job_path)

            with open(job_path) as stream:
                content = stream.read()

            # Check collection_type is shown in comments
            assert "collection_type: list" in content
            assert "collection_type: paired" in content
            assert "collection_type: list:paired" in content

            with open(job_path) as stream:
                job = yaml.safe_load(stream)

            # Check list collection
            list_input = job.get("list_input")
            assert list_input.get("class") == "Collection"
            assert list_input.get("collection_type") == "list"
            assert len(list_input.get("elements")) == 1
            assert list_input["elements"][0]["identifier"] == "todo_element_name"

            # Check paired collection - should have forward and reverse
            paired_input = job.get("paired_input")
            assert paired_input.get("class") == "Collection"
            assert paired_input.get("collection_type") == "paired"
            assert len(paired_input.get("elements")) == 2
            identifiers = [e["identifier"] for e in paired_input["elements"]]
            assert "forward" in identifiers
            assert "reverse" in identifiers

            # Check list:paired collection - should have nested structure
            list_paired_input = job.get("list_paired_input")
            assert list_paired_input.get("class") == "Collection"
            assert list_paired_input.get("collection_type") == "list:paired"
            assert len(list_paired_input.get("elements")) == 1
            nested = list_paired_input["elements"][0]
            assert nested.get("class") == "Collection"
            assert nested.get("type") == "paired"
            assert len(nested.get("elements")) == 2
            nested_identifiers = [e["identifier"] for e in nested["elements"]]
            assert "forward" in nested_identifiers
            assert "reverse" in nested_identifiers
