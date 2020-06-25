"""Autoupdate older conda dependencies in the requirements section."""
from __future__ import absolute_import

import planemo.conda


def begin_update(tool_path, tool_xml):
    version_dict, main_req_dict = find_packages_in_requirements(tool_xml)
    updated_version_dict = get_latest_versions(version_dict)
    return update_requirements(tool_path, tool_xml, updated_version_dict)


def find_packages_in_requirements(tool_xml):
    packages = {}
    main_req_dict = {}
    for element in tool_xml.xml_tree.findall("requirements"):
        for requirement in list(element):
            if requirement.tag == 'requirement' and requirement.attrib.get('type', '') == 'package':
                packages[requirement.text] = requirement.attrib.get('version', '')
    return packages, main_req_dict


def get_latest_versions(version_dict):

    for package in version_dict.keys():
        target = planemo.conda.conda_util.CondaTarget(package)
        search_results = planemo.conda.best_practice_search(target)
        version_dict[package] = search_results[0]['version']
    return version_dict


def update_requirements(tool_path, tool_xml, updated_version_dict):
    for element in tool_xml.xml_tree.findall("requirements"):
        for requirement in list(element):
            if requirement.tag == 'requirement' and requirement.attrib.get('type', '') == 'package':
                requirement.set('version', updated_version_dict[requirement.text])
    tool_xml.xml_tree.write(tool_path)
    return tool_xml


__all__ = (
    "begin_update"
)
