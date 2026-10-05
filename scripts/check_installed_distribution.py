#!/usr/bin/env python
"""Smoke-test an installed wheel from outside the source checkout."""

import argparse
import importlib.metadata
import subprocess
import tempfile
from pathlib import Path

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("distribution", choices=["planemo", "planemo-cli"])
    args = parser.parse_args()
    distribution = importlib.metadata.distribution(args.distribution)
    if args.distribution == "planemo-cli":
        expected = {"planemo-cli"}
        for value in distribution.requires:
            requirement = Requirement(value)
            if requirement.marker is None or requirement.marker.evaluate():
                expected.add(canonicalize_name(requirement.name))
                assert importlib.metadata.version(requirement.name) in requirement.specifier, value
        installed = {canonicalize_name(package.metadata["Name"]) for package in importlib.metadata.distributions()}
        assert installed == expected, f"Installed dependency closure differs: {installed ^ expected}"
        try:
            importlib.metadata.distribution("planemo")
        except importlib.metadata.PackageNotFoundError:
            pass
        else:
            raise AssertionError("Both alternative distributions are installed")
    with tempfile.TemporaryDirectory(prefix="planemo-installed-test-") as directory:
        subprocess.run(["planemo", "--version"], cwd=directory, check=True)
        subprocess.run(["planemo", "--help"], cwd=directory, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(
            ["planemo", "lint", str(ROOT / "tests/data/tools/ok_conditional.xml")], cwd=directory, check=True
        )
        subprocess.run(
            [
                "planemo",
                "test_reports",
                str(ROOT / "tests/data/issue381.json"),
                "--test_output",
                str(Path(directory) / "report.html"),
            ],
            cwd=directory,
            check=True,
        )
        assert (Path(directory) / "report.html").stat().st_size > 0


if __name__ == "__main__":
    main()
