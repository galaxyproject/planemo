#!/usr/bin/env python
# Little script to make HISTORY.rst more easy to format properly, lots TODO
# pull message down and embed, use arg parse, handle multiple, etc...
import os
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


def main(argv):
    history_path = os.path.join(PROJECT_DIRECTORY, "HISTORY.rst")
    with open(history_path, encoding="utf-8") as fh:
        history = fh.read()

    def extend(from_str, line):
        from_str += "\n"
        return history.replace(from_str, from_str + line + "\n")

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
