#!/usr/bin/env python
# Little script to make HISTORY.rst more easy to format properly, lots TODO
# pull message down and embed, use arg parse, handle multiple, etc...
import os
import sys
try:
    import requests
except ImportError:
    requests = None
import textwrap

from six.moves.urllib.parse import urljoin

PROJECT_DIRECTORY = os.path.join(os.path.dirname(__file__), "..")
new_path = [PROJECT_DIRECTORY]
new_path.extend(sys.path[1:])  # remove scripts/ from the path
sys.path = new_path

import planemo as project

PROJECT_OWNER = project.PROJECT_OWNER
PROJECT_NAME = project.PROJECT_NAME
PROJECT_URL = "https://github.com/%s/%s" % (PROJECT_OWNER, PROJECT_NAME)
PROJECT_API = "https://api.github.com/repos/%s/%s/" % (PROJECT_OWNER, PROJECT_NAME)
AUTHORS_SKIP_CREDIT = ["jmchilton", "mvdbeek", "nsoranzo", "bgruening", "natefoo"]


def main(argv):
    history_path = os.path.join(PROJECT_DIRECTORY, "HISTORY.rst")
    history = open(history_path, "r").read()

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
        pull_request = ident[len("pr"):]
        api_url = urljoin(PROJECT_API, "pulls/%s" % pull_request)
        req = requests.get(api_url).json()
        message = req["title"]
        login = req["user"]["login"]
        if login not in AUTHORS_SKIP_CREDIT:
            message = message.rstrip(".")
            message += " (thanks to `@%s`_)." % req["user"]["login"]
    elif requests is not None and ident.startswith("issue"):
        issue = ident[len("issue"):]
        api_url = urljoin(PROJECT_API, "issues/%s" % issue)
        req = requests.get(api_url).json()
        message = req["title"]
    else:
        message = ""

    to_doc = message + " "

    if ident.startswith("pr"):
        pull_request = ident[len("pr"):]
        text = ".. _Pull Request {0}: {1}/pull/{0}".format(pull_request, PROJECT_URL)
        history = extend(".. github_links", text)
        to_doc += "`Pull Request {0}`_".format(pull_request)
    elif ident.startswith("issue"):
        issue = ident[len("issue"):]
        text = ".. _Issue {0}: {1}/issues/{0}".format(issue, PROJECT_URL)
        history = extend(".. github_links", text)
        to_doc += "`Issue {0}`_".format(issue)
    else:
        short_rev = ident[:7]
        text = ".. _{0}: {1}/commit/{0}".format(short_rev, PROJECT_URL)
        history = extend(".. github_links", text)
        to_doc += "{0}_".format(short_rev)

    to_doc = wrap(to_doc)
    history = extend(".. to_doc", to_doc)
    open(history_path, "w").write(history)


def get_first_sentence(message):
    first_line = message.split("\n")[0]
    return first_line


def wrap(message):
    wrapper = textwrap.TextWrapper(initial_indent="* ")
    wrapper.subsequent_indent = '  '
    wrapper.width = 78
    return "\n".join(wrapper.wrap(message))


if __name__ == "__main__":
    main(sys.argv)
