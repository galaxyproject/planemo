import pathlib
import subprocess
import sys
from typing import (
    List,
    Optional,
)

SCHEMA = pathlib.Path(__file__).parent / "test_file_schema.json"

def validate_schema(test_files: List[str], verbose: bool = False) -> Optional[str]:
    """
    Runs check_jsonschema on `test_files`.

    Returns validation failure message if validation failed.
    """
    check_args = [sys.executable, "-m", "check_jsonschema", "--schemafile", str(SCHEMA), *test_files]
    if verbose:
        check_args.append("--verbose")
    result = subprocess.run(check_args, capture_output=True, text=True)
    if result.returncode:
        return result.stdout
    