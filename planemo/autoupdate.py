"""Autoupdate older conda dependencies in the requirements section."""
from __future__ import absolute_import

import planemo.conda


def begin_update(tool_path, tool_xml):
    version_dict = find_packages_in_requirements(tool_xml)
    updated_version_dict = get_latest_versions(version_dict)
    macros = get_macros(tool_xml)
    return update_requirements(tool_path, tool_xml, updated_version_dict)


def find_packages_in_requirements(tool_xml):
    packages = {}
    for element in tool_xml.xml_tree.findall("requirements"):
        for requirement in list(element):
            packages[requirement.text] = requirement.attrib.get('version')
            #print(requirement)
    return packages


def get_latest_versions(version_dict):

    for package in version_dict.keys():
        target = planemo.conda.conda_util.CondaTarget(package)
        search_results = planemo.conda.best_practice_search(target)
        version_dict[package] = search_results[0]['version']
    #print(version_dict)
    return version_dict


def get_macros(tool_xml):
    for element in tool_xml.xml_tree.findall("macros"):
        for macro in list(element):
            print(macro.attrib)


def update_requirements(tool_path, tool_xml, updated_version_dict):
    for element in tool_xml.xml_tree.findall("requirements"):
        for requirement in list(element):
            requirement.set('version', updated_version_dict[requirement.text])
    tool_xml.xml_tree.write(tool_path)
    return tool_xml


__all__ = (
    "begin_update"
)
