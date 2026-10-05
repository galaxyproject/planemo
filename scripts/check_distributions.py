#!/usr/bin/env python
"""Verify the built distributions have identical payloads and complete runtime pins."""

import argparse
from email.parser import BytesParser
from pathlib import Path
from zipfile import ZipFile

from packaging.markers import default_environment
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name


def wheel_contents(path):
    with ZipFile(path) as archive:
        metadata_name = next(name for name in archive.namelist() if name.endswith(".dist-info/METADATA"))
        metadata = BytesParser().parsebytes(archive.read(metadata_name))
        payload = {name: archive.read(name) for name in archive.namelist() if ".dist-info/" not in name}
        entry_points = archive.read(metadata_name.replace("METADATA", "entry_points.txt"))
    return metadata, payload, entry_points


def check_distributions(output):
    normal, payload, entry_points = wheel_contents(next(output.glob("planemo-*.whl")))
    cli, cli_payload, cli_entry_points = wheel_contents(next(output.glob("planemo_cli-*.whl")))
    assert normal["Name"] == "planemo"
    assert cli["Name"] == "planemo-cli"
    assert normal["Version"] == cli["Version"]
    assert normal["Requires-Python"] == cli["Requires-Python"]
    assert payload == cli_payload, "Distribution code and assets differ"
    assert entry_points == cli_entry_points
    for resource in (
        "planemo/xml/xsd/repository_dependencies.xsd",
        "planemo/xml/xsd/tool_dependencies.xsd",
        "planemo/reports/report_html.tpl",
    ):
        assert resource in payload, f"Missing installed resource: {resource}"
    requirements = [Requirement(value) for value in cli.get_all("Requires-Dist", [])]
    assert requirements, "CLI distribution has no dependencies"
    for requirement in requirements:
        assert len(requirement.specifier) == 1 and next(iter(requirement.specifier)).operator == "==", requirement
        assert canonicalize_name(requirement.name) != "planemo", "CLI must contain its own code"
    expected_pins = {
        str(Requirement(line))
        for line in (Path(__file__).resolve().parents[1] / "requirements-cli.txt").read_text().splitlines()
        if line and not line.startswith("#")
    }
    assert {str(requirement) for requirement in requirements} == expected_pins
    for minor in range(10, 15):
        for system, platform in (("Linux", "linux"), ("Darwin", "darwin"), ("Windows", "win32")):
            environment = dict(
                default_environment(),
                python_version=f"3.{minor}",
                python_full_version=f"3.{minor}.0",
                platform_system=system,
                sys_platform=platform,
            )
            active = [r for r in requirements if r.marker is None or r.marker.evaluate(environment)]
            pins = {canonicalize_name(r.name): next(iter(r.specifier)).version for r in active}
            assert len(pins) == len(active), f"Overlapping pins for {environment}"
            for value in normal.get_all("Requires-Dist", []):
                requirement = Requirement(value)
                if requirement.marker is None or requirement.marker.evaluate(environment):
                    version = pins.get(canonicalize_name(requirement.name))
                    assert version is not None and version in requirement.specifier, requirement
    print("Both distributions contain identical code/assets and valid runtime pins")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parents[1] / "dist")
    args = parser.parse_args()
    check_distributions(args.outdir.resolve())


if __name__ == "__main__":
    main()
