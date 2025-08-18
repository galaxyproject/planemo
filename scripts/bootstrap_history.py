#!/usr/bin/env python
# Little script to make HISTORY.rst more easy to format properly, lots TODO
# pull message down and embed, use arg parse, handle multiple, etc...
import os
import re
import subprocess
import sys
import textwrap
from urllib.parse import urljoin

import requests

PROJECT_DIRECTORY = os.path.join(os.path.dirname(__file__), "..")
new_path = [PROJECT_DIRECTORY]
new_path.extend(sys.path[1:])  # remove scripts/ from the path
sys.path = new_path

import planemo as project  # noqa: E402

PROJECT_AUTHOR = project.PROJECT_AUTHOR
PROJECT_NAME = project.PROJECT_NAME
PROJECT_URL = f"https://github.com/{PROJECT_AUTHOR}/{PROJECT_NAME}"
PROJECT_API = f"https://api.github.com/repos/{PROJECT_AUTHOR}/{PROJECT_NAME}/"


def get_last_release_tag():
    """Get the last release tag based on the current version in __init__.py"""
    version = project.__version__
    # Remove .dev0 suffix if present
    if ".dev" in version:
        version = version.split(".dev")[0]

    # Parse version components
    parts = version.split(".")
    if len(parts) >= 3:
        major, minor, patch = parts[:3]
        # Decrement patch version to get last release
        last_patch = max(0, int(patch) - 1)
        return f"{major}.{minor}.{last_patch}"
    return version


def get_merge_commits_since_tag(tag):
    """Get merge commits since the specified tag"""
    try:
        result = subprocess.run(
            ["git", "log", "--merges", "--oneline", f"{tag}..HEAD"], capture_output=True, text=True, check=True
        )
        return result.stdout.strip().split("\n") if result.stdout.strip() else []
    except subprocess.CalledProcessError:
        return []


def extract_pr_info_from_merge(merge_line):
    """Extract PR number and author from merge commit message"""
    # Match pattern: "Merge pull request #1234 from author/branch"
    match = re.match(r"[a-f0-9]+\s+Merge pull request #(\d+) from ([^/]+)/", merge_line)
    if match:
        pr_number = match.group(1)
        author = match.group(2)
        return pr_number, author
    return None, None


def generate_acknowledgements():
    """Generate acknowledgement lines for merge commits since last release"""
    tag = get_last_release_tag()
    merge_commits = get_merge_commits_since_tag(tag)

    acknowledgements = []
    for merge in merge_commits:
        if merge.strip():
            pr_number, author = extract_pr_info_from_merge(merge)
            if pr_number and author:
                try:
                    # Get PR details from GitHub API
                    api_url = urljoin(PROJECT_API, f"pulls/{pr_number}")
                    req = requests.get(api_url).json()
                    title = req.get("title", "")
                    login = req["user"]["login"]

                    # Format acknowledgement line
                    title_clean = title.rstrip(".")
                    ack_line = f"* {title_clean} (thanks to `@{login}`_). `Pull Request {pr_number}`_"
                    acknowledgements.append(ack_line)

                    # Add GitHub link
                    github_link = f".. _Pull Request {pr_number}: {PROJECT_URL}/pull/{pr_number}"
                    acknowledgements.append(github_link)
                except Exception as e:
                    print(f"Error processing PR {pr_number}: {e}", file=sys.stderr)

    return acknowledgements


def main(argv):
    history_path = os.path.join(PROJECT_DIRECTORY, "HISTORY.rst")
    with open(history_path, encoding="utf-8") as fh:
        history = fh.read()

    def extend(from_str, line):
        from_str += "\n"
        return history.replace(from_str, from_str + line + "\n")

    # Check if we should generate acknowledgements for merge commits
    if len(argv) > 1 and argv[1] == "--acknowledgements":
        acknowledgements = generate_acknowledgements()
        if acknowledgements:
            print("Generated acknowledgement lines:")
            for ack in acknowledgements:
                if ack.startswith("*"):
                    print(ack)
                    history = extend(".. to_doc", ack)
                elif ack.startswith(".."):
                    print(ack)
                    history = extend(".. github_links", ack)

            with open(history_path, "w", encoding="utf-8") as fh:
                fh.write(history)
            print(f"\nAcknowledgements added to {history_path}")
        else:
            print("No merge commits found since last release.")
        return

    if len(argv) < 2:
        print("Usage: python bootstrap_history.py <identifier> [message]")
        print("   or: python bootstrap_history.py --acknowledgements")
        return

    ident = argv[1]

    message = ""
    if len(argv) > 2:
        message = argv[2]
    elif not (ident.startswith("pr") or ident.startswith("issue")):
        api_url = urljoin(PROJECT_API, "commits/%s" % ident)
        req = requests.get(api_url).json()
        commit = req["commit"]
        message = commit["message"]
        message = get_first_sentence(message)
    elif requests is not None and ident.startswith("pr"):
        pull_request = ident[len("pr") :]
        api_url = urljoin(PROJECT_API, "pulls/%s" % pull_request)
        req = requests.get(api_url).json()
        message = req["title"]
        login = req["user"]["login"]
        message = message.rstrip(".")
        message += f" (thanks to `@{login}`_)."
    elif requests is not None and ident.startswith("issue"):
        issue = ident[len("issue") :]
        api_url = urljoin(PROJECT_API, "issues/%s" % issue)
        req = requests.get(api_url).json()
        message = req["title"]
    else:
        message = ""

    to_doc = message + " "

    if ident.startswith("pr"):
        pull_request = ident[len("pr") :]
        text = ".. _Pull Request {0}: {1}/pull/{0}".format(pull_request, PROJECT_URL)
        history = extend(".. github_links", text)
        to_doc += f"`Pull Request {pull_request}`_"
    elif ident.startswith("issue"):
        issue = ident[len("issue") :]
        text = ".. _Issue {0}: {1}/issues/{0}".format(issue, PROJECT_URL)
        history = extend(".. github_links", text)
        to_doc += f"`Issue {issue}`_"
    else:
        short_rev = ident[:7]
        text = ".. _{0}: {1}/commit/{0}".format(short_rev, PROJECT_URL)
        history = extend(".. github_links", text)
        to_doc += f"{short_rev}_"

    to_doc = wrap(to_doc)
    history = extend(".. to_doc", to_doc)
    with open(history_path, "w", encoding="utf-8") as fh:
        fh.write(history)


def get_first_sentence(message):
    first_line = message.split("\n")[0]
    return first_line


def wrap(message):
    wrapper = textwrap.TextWrapper(initial_indent="* ")
    wrapper.subsequent_indent = "  "
    wrapper.width = 78
    return "\n".join(wrapper.wrap(message))


if __name__ == "__main__":
    main(sys.argv)
