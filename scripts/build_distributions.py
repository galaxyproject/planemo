#!/usr/bin/env python
"""Build Planemo and its pinned standalone CLI from the same source distribution."""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import tomlkit

ROOT = Path(__file__).resolve().parents[1]


def cli_project(source):
    path = source / "pyproject.toml"
    config = tomlkit.loads(path.read_text())
    config["project"]["name"] = "planemo-cli"
    del config["project"]["dependencies"]
    config["project"]["dynamic"].append("dependencies")
    config["tool"]["setuptools"]["dynamic"]["dependencies"] = {"file": "requirements-cli.txt"}
    path.write_text(tomlkit.dumps(config))
    # Source lists and metadata must be regenerated for the new distribution name.
    for metadata in source.glob("*.egg-info"):
        shutil.rmtree(metadata)
    (source / "PKG-INFO").unlink(missing_ok=True)


def build_distributions(output):
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, "-m", "build", "--outdir", str(output), str(ROOT)], check=True)
    sdists = list(output.glob("planemo-*.tar.gz"))
    if len(sdists) != 1:
        raise RuntimeError("Expected exactly one Planemo source distribution; use a clean output directory")
    with tempfile.TemporaryDirectory(prefix="planemo-cli-build-") as temporary:
        # This archive was produced immediately above by our own build backend.
        shutil.unpack_archive(str(sdists[0]), temporary)
        source = next(Path(temporary).iterdir())
        cli_project(source)
        subprocess.run([sys.executable, "-m", "build", "--outdir", str(output), str(source)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=ROOT / "dist")
    args = parser.parse_args()
    build_distributions(args.outdir.resolve())


if __name__ == "__main__":
    main()
