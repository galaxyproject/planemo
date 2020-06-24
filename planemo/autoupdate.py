"""Autoupdate older conda dependencies in the requirements section."""
from __future__ import absolute_import

import planemo.conda


def autoupdate(tool_xml):
    version_dict, main_req_dict = find_packages_in_requirements(tool_xml)
    main_req_dict = get_main_requirement_version(tool_xml, main_req_dict)
    updated_version_dict, main_req_updated_version = get_latest_versions(version_dict, main_req_dict)
    update_requirements(tool_xml, updated_version_dict, main_req_updated_version)


def find_packages_in_requirements(tool_xml):
    packages = {}
    main_req_dict = {}
    main_req = None
    print("fas")
    print(tool_xml)
    for element in tool_xml.getroot().findall("requirements"):
        for requirement in list(element):
            if requirement.tag == 'requirement' and requirement.attrib.get('type', '') == 'package':
                if requirement.attrib.get('version', '') == '@TOOL_VERSION@':
                    main_req = requirement.text
                packages[requirement.text] = requirement.attrib.get('version', '')
    main_req_dict['main_req'] = main_req
    return packages, main_req_dict


def get_main_requirement_version(tool_xml, main_req_dict):
    for element in tool_xml.getroot().findall("macros"):
        for macro in list(element):
            if macro.tag == 'token' and macro.attrib.get('name', '') == '@TOOL_VERSION@':
                main_req_dict['@TOOL_VERSION@'] = macro.text
                return main_req_dict
    return None


def get_latest_versions(version_dict, main_req_dict):

    for package in version_dict.keys():
        target = planemo.conda.conda_util.CondaTarget(package)
        search_results = planemo.conda.best_practice_search(target)
        version_dict[package] = search_results[0]['version']

        main_req_dict['@TOOL_VERSION@'] = version_dict[main_req_dict['main_req']]
    return version_dict, main_req_dict


def update_requirements(tool_xml, updated_version_dict, main_req_updated_version):
    for element in tool_xml.getroot().findall("requirements"):
        for requirement in list(element):
            if requirement.tag == 'requirement' and requirement.attrib.get('type', '') == 'package':
                if requirement.text == main_req_updated_version['main_req']:
                    requirement.set('version', '@TOOL_VERSION@')
                else:
                    requirement.set('version', updated_version_dict[requirement.text])

    for element in tool_xml.getroot().findall("macros"):
        for macro in list(element):
            if macro.tag == 'token' and macro.attrib.get('name', '') == '@TOOL_VERSION@':
                macro.text = main_req_updated_version['@TOOL_VERSION@']
