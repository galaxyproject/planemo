#!/usr/bin/env python
"""
Script to create a new version by incrementing version numbers and updating files.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

try:
    from packaging.version import Version
except ImportError:
    print("Error: packaging library is required. Install with: pip install packaging")
    sys.exit(1)

PROJECT_DIRECTORY = Path(__file__).parent.parent


class VersionBumper:
    """Handles version bumping operations."""

    def __init__(self, project_dir: Path):
        self.project_dir = project_dir
        self.history_path = project_dir / "HISTORY.rst"

    def get_current_version(self, source_dir: str) -> str:
        """Extract current version from __init__.py file."""
        init_path = self.project_dir / source_dir / "__init__.py"
        if not init_path.exists():
            raise FileNotFoundError(f"Cannot find {init_path}")

        content = init_path.read_text(encoding="utf-8")
        version_match = re.search(r'__version__ = ["\']([^"\']+)["\']', content)
        if not version_match:
            raise ValueError(f"Cannot find version in {init_path}")

        return version_match.group(1)

    def increment_version(self, version_str: str, bump_type: str = "patch") -> str:
        """
        Increment version based on bump type.

        Args:
            version_str: Current version string (e.g., "1.2.3" or "1.2.3.dev0")
            bump_type: Type of version bump ("major", "minor", "patch")

        Returns:
            New version string
        """
        # Remove .dev0 suffix if present
        clean_version = version_str.split(".dev")[0]

        try:
            Version(clean_version)
        except Exception as e:
            raise ValueError(f"Invalid version format '{clean_version}': {e}")

        # Parse version components
        parts = clean_version.split(".")
        if len(parts) < 3:
            # Pad to at least 3 components
            parts.extend(["0"] * (3 - len(parts)))

        major, minor, patch = int(parts[0]), int(parts[1]), int(parts[2])

        if bump_type == "major":
            major += 1
            minor = 0
            patch = 0
        elif bump_type == "minor":
            minor += 1
            patch = 0
        elif bump_type == "patch":
            patch += 1
        else:
            raise ValueError(f"Invalid bump type: {bump_type}. Use 'major', 'minor', or 'patch'")

        return f"{major}.{minor}.{patch}"

    def update_history_file(self, new_version: str) -> None:
        """Add new development version section to HISTORY.rst."""
        if not self.history_path.exists():
            raise FileNotFoundError(f"Cannot find {self.history_path}")

        history_content = self.history_path.read_text(encoding="utf-8")

        # Find the insertion point after ".. to_doc"
        to_doc_marker = ".. to_doc\n"
        if to_doc_marker not in history_content:
            raise ValueError("Cannot find '.. to_doc' marker in HISTORY.rst")

        new_section = f"""
---------------------
{new_version}.dev0
---------------------

    """

        updated_content = history_content.replace(to_doc_marker, to_doc_marker + new_section)

        self.history_path.write_text(updated_content, encoding="utf-8")
        print(f"Updated {self.history_path} with new version section")

    def update_init_file(self, source_dir: str, new_version: str) -> None:
        """Update version in __init__.py file."""
        init_path = self.project_dir / source_dir / "__init__.py"
        if not init_path.exists():
            raise FileNotFoundError(f"Cannot find {init_path}")

        content = init_path.read_text(encoding="utf-8")

        # Update version string
        updated_content = re.sub(r'(__version__ = ["\'])[^"\']+(["\'])', rf"\g<1>{new_version}.dev0\g<2>", content)

        if updated_content == content:
            raise ValueError("Failed to update version in __init__.py")

        init_path.write_text(updated_content, encoding="utf-8")
        print(f"Updated {init_path} with version {new_version}.dev0")

    def commit_changes(self, source_dir: str, new_version: str) -> None:
        """Commit the version changes to git."""
        files_to_commit = ["HISTORY.rst", f"{source_dir}/__init__.py"]

        try:
            cmd = ["git", "commit", "-m", f"Starting work on {new_version}"] + files_to_commit
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"Committed changes: Starting work on {new_version}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to commit changes: {e}")
            print(f"stdout: {e.stdout}")
            print(f"stderr: {e.stderr}")
            raise

    def bump_version(self, source_dir: str, bump_type: str = "patch") -> str:
        """
        Complete version bump process.

        Args:
            source_dir: Source directory containing __init__.py
            bump_type: Type of version bump ("major", "minor", "patch")

        Returns:
            New version string
        """
        print(f"Bumping {bump_type} version in {source_dir}/")

        # Get current version
        current_version = self.get_current_version(source_dir)
        print(f"Current version: {current_version}")

        # Calculate new version
        new_version = self.increment_version(current_version, bump_type)
        print(f"New version: {new_version}")

        # Update files
        self.update_history_file(new_version)
        self.update_init_file(source_dir, new_version)

        # Commit changes
        self.commit_changes(source_dir, new_version)

        return new_version


def create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Create a new version by incrementing version numbers and updating files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s planemo                    # Bump patch version (1.2.3 -> 1.2.4)
  %(prog)s planemo --major            # Bump major version (1.2.3 -> 2.0.0)
  %(prog)s planemo --minor            # Bump minor version (1.2.3 -> 1.3.0)
  %(prog)s planemo --patch            # Bump patch version (1.2.3 -> 1.2.4)
        """,
    )

    parser.add_argument("source_dir", help="Source directory containing __init__.py (e.g., 'planemo')")

    # Version bump type (mutually exclusive)
    bump_group = parser.add_mutually_exclusive_group()
    bump_group.add_argument(
        "--major", action="store_const", const="major", dest="bump_type", help="Bump major version (X.y.z -> X+1.0.0)"
    )
    bump_group.add_argument(
        "--minor", action="store_const", const="minor", dest="bump_type", help="Bump minor version (x.Y.z -> x.Y+1.0)"
    )
    bump_group.add_argument(
        "--patch",
        action="store_const",
        const="patch",
        dest="bump_type",
        help="Bump patch version (x.y.Z -> x.y.Z+1) [default]",
    )

    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without making changes")

    # Set default bump type to patch
    parser.set_defaults(bump_type="patch")

    return parser


def main() -> None:
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()

    try:
        bumper = VersionBumper(PROJECT_DIRECTORY)

        if args.dry_run:
            current_version = bumper.get_current_version(args.source_dir)
            new_version = bumper.increment_version(current_version, args.bump_type)
            print(f"DRY RUN: Would bump {args.bump_type} version")
            print(f"Current version: {current_version}")
            print(f"New version would be: {new_version}")
            return

        new_version = bumper.bump_version(args.source_dir, args.bump_type)
        print(f"\n✅ Successfully bumped version to {new_version}.dev0")

    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
