"""Autoupdate older conda dependencies in the requirements section."""
from __future__ import absolute_import

import re
import xml.etree.ElementTree as ET
import planemo.conda

def begin_update(tool_path, tool_xml):
    tool_xml.xml_tree = ET.parse(tool_path)
    macros = get_tokens(tool_xml)

    version_dict = find_packages_in_requirements(tool_xml)
    updated_main_req = get_latest_versions({version_dict['main_req']: tokens['@TOOL_VERSION@']})
    if updated_main_req[version_dict['main_req']] == tokens['@TOOL_VERSION@']:
        return

    updated_version_dict = get_latest_versions(version_dict['other_reqs'])
    return update_requirements(tool_path, tool_xml, updated_version_dict, updated_main_req)


def find_packages_in_requirements(tool_xml):
    """
    Get all requirements with versions as a dictionary of the form
    {'main_req': main requirement, 'other_reqs': {req1: version, ... }}
    """
    packages = {'other_reqs': {}}
    for element in tool_xml.xml_tree.findall("requirements"):
        for requirement in list(element):
            if requirement.attrib.get('version') == '@TOOL_VERSION@':
                packages['main_req'] = requirement.text
            else:
                packages['other_reqs'][requirement.text] = requirement.attrib.get('version')
    return packages


def get_latest_versions(version_dict):
    """
    Update a dict with current conda versions
    """
    for package in version_dict.keys():
        target = planemo.conda.conda_util.CondaTarget(package)
        search_results = planemo.conda.best_practice_search(target)
        version_dict[package] = search_results[0]['version']
    return version_dict


def get_tokens(tool_xml):
    """
    Get tokens used to define versions
    """
    macros = {}
    for element in tool_xml.xml_tree.findall("macros"):
        for macro in list(element):
            if 'VERSION@' in macro.attrib['name']:
                macros[macro.attrib['name']] = macro.text
    return macros

def update_requirements(tool_path, tool_xml, updated_version_dict, updated_main_req):
    """
    Update requirements to latest conda versions
    and update version tokens
    """
    for element in tool_xml.xml_tree.findall("macros"):
        for token in list(element):
            if token.attrib.get('name') == '@TOOL_VERSION@':
                token.text = updated_main_req.values()[0]
            elif token.attrib.get('name') == '@GALAXY_VERSION@':
                token.text = '0'
            else:
                continue

    for element in tool_xml.xml_tree.findall("requirements"):
        for requirement in list(element):
            if requirement.text not in updated_main_req:
                requirement.set('version', updated_version_dict[requirement.text])
    print(tool_path)
    print(tool_xml)

    write_to_xml(tool_path, tool_xml)
    return tool_xml


def write_to_xml(tool_path, tool_xml):
    """
    Write modified XML to tool_path
    """
    # macros
    m = tool_xml.xml_tree.find('macros')
    m_str = ET.tostring(m).strip()

    # requirements
    r = tool_xml.xml_tree.find('requirements')
    r_str = ET.tostring(r).strip()

    # write to file
    with open(tool_path, 'r+') as f:
        xml_text = f.read()
        xml_text = re.sub('<macros>(.|\n)*</macros>', m_str, xml_text)
        xml_text = re.sub('<requirements>(.|\n)*</requirements>', r_str, xml_text)
        f.seek(0)
        f.truncate()
        f.write(xml_text)


__all__ = (
    "begin_update"
)
