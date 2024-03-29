"""Module contains code for the Requirement, Reference and some general functions for training."""

import collections

import oyaml as yaml


class Requirement:
    """Class to describe a training requirement."""

    def __init__(self, req_type="internal", topic_name="introduction", title=None, tutorials=None, link=None):
        """Init a Requirement instance."""
        self.type = req_type
        self.topic_name = topic_name
        self.tutorials = tutorials
        self.title = title
        self.link = link

    def init_from_dict(self, metadata):
        """Init from a dictionary generated by export_to_ordered_dict."""
        self.type = metadata["type"]
        if self.type == "internal":
            self.topic_name = metadata["topic_name"]
            if "tutorials" in metadata:
                self.tutorials = metadata["tutorials"]
        else:
            self.title = metadata["title"]
            if self.type == "external":
                self.link = metadata["link"]

    def export_to_ordered_dict(self):
        """Export the requirement into an ordered dictionary."""
        req = collections.OrderedDict()
        req["type"] = self.type
        if self.type == "internal":
            req["topic_name"] = self.topic_name
            if self.tutorials:
                req["tutorials"] = self.tutorials
        else:
            req["title"] = self.title
            if self.type == "external":
                req["link"] = self.link
        return req


def load_yaml(filepath):
    """Load the content of a YAML file to a dictionary."""
    with open(filepath) as m_file:
        content = yaml.safe_load(m_file)
    return content


def save_to_yaml(content, filepath):
    """Save a dictionary to a YAML file."""
    with open(filepath, "w") as stream:
        yaml.safe_dump(
            content,
            stream,
            indent=2,
            default_flow_style=False,
            default_style="",
            explicit_start=True,
            encoding="utf-8",
            allow_unicode=True,
        )
