"""Autoupdate older conda dependencies in the requirements section."""
from __future__ import absolute_import

import re
import xml.etree.ElementTree as ET
import planemo.conda

from galaxy.tool_util.deps import conda_util
from planemo.io import error, info

def autoupdate(ctx, tool_path, **kwds):
    """
    Autoupdate an XML file
    """
    xml_tree = ET.parse(tool_path)
    requirements = find_requirements(xml_tree)

    for macro_import in requirements['imports']:
        # recursively check macros
        macro_requirements = autoupdate('/'.join(tool_path.split('/')[:-1] + [macro_import]), dry_run)
        for requirement in macro_requirements:
            if requirement not in requirements:
                requirements[requirement] = macro_requirements[requirement]

    if not requirements.get('@TOOL_VERSION@'):
        # if tool_version is not specified, finish without changes
        error("The @TOOL_VERSION@ token is not specified in {}. This is required for autoupdating.".format(tool_path))
        return requirements
    updated_main_req = get_latest_versions(ctx, {requirements.get('main_req'): requirements.get('@TOOL_VERSION@')}, **kwds)
    if updated_main_req[requirements.get('main_req')] == requirements.get('@TOOL_VERSION@'):
        # check main_req is up-to-date; if so, finish without changes
        info("No updates required or made to {}.".format(tool_path))
        return requirements

    if kwds.get('dry_run'):
        error("Update required to {}! Tool main requirement has version {}, newest conda version is {}".format(tool_path, requirements.get('@TOOL_VERSION@'), updated_main_req[requirements.get('main_req')]))
        return requirements

    # if main_req is not up-to-date, update everything
    updated_version_dict = get_latest_versions(ctx, requirements.get('other_reqs'), **kwds)
    update_requirements(tool_path, xml_tree, updated_version_dict, updated_main_req)
    info("Tool {} updated.".format(tool_path))
    return requirements


def find_requirements(xml_tree):
    """
    Get all requirements with versions as a dictionary of the form
    {'main_req': main requirement, 'other_reqs': {req1: version, ... },
        'imports: ['macros.xml'], '*VERSION@': '...'}
    """
    requirements = {'other_reqs': {}, 'imports': []}

    # get tokens
    for token in xml_tree.iter("token"):
        if 'VERSION@' in token.attrib.get('name', ''):
            requirements[token.attrib['name']] = token.text
    for macro_import in xml_tree.iter("import"):
        requirements['imports'].append(macro_import.text)

    # get requirements
    for requirement in xml_tree.iter("requirement"):
        if requirement.attrib.get('version') == '@TOOL_VERSION@':
            requirements['main_req'] = requirement.text
        else:
            requirements['other_reqs'][requirement.text] = requirement.attrib.get('version')
    return requirements


def get_latest_versions(ctx, version_dict, **kwds):
    """
    Update a dict with current conda versions for tool requirements
    """
    for package in version_dict.keys():
        conda_context = planemo.conda.build_conda_context(ctx, **kwds)
        target = planemo.conda.conda_util.CondaTarget(package)
        search_results = conda_util.best_search_result(target, conda_context=conda_context)
        if search_results[1]:
            version_dict[package] = search_results[0]['version']
        else:
            error('Requirements could not be found!')
            exit()
    return version_dict


def update_requirements(tool_path, xml_tree, updated_version_dict, updated_main_req):
    """
    Update requirements to latest conda versions
    and update version tokens
    """

    tags_to_update = {'tokens': [], 'requirements': []}

    for token in xml_tree.iter("token"):
        if token.attrib.get('name') == '@TOOL_VERSION@':
            #  check this
            token.text = list(updated_main_req.values())[0]
        elif token.attrib.get('name') == '@GALAXY_VERSION@':
            token.text = '0'
        else:
            continue
        tags_to_update['tokens'].append(ET.tostring(token, encoding="unicode").strip())

    if '@GALAXY_VERSION@' not in [n.attrib.get('name') for n in xml_tree.iter('token')]:
        tags_to_update['update_tool'] = True

    for requirement in xml_tree.iter("requirement"):
        if requirement.text not in updated_main_req:
            requirement.set('version', updated_version_dict[requirement.text])
        tags_to_update['requirements'].append(ET.tostring(requirement, encoding="unicode").strip())
    write_to_xml(tool_path, xml_tree, tags_to_update)
    return xml_tree


def write_to_xml(tool_path, xml_tree, tags_to_update):
    """
    Write modified XML to tool_path
    """
    with open(tool_path, 'r+') as f:
        xml_text = f.read()
        for token in tags_to_update['tokens']:
            xml_text = re.sub('{}>.*<{}'.format(*re.split('>.*<', token)), token, xml_text)

        for requirement in tags_to_update['requirements']:
            xml_text = re.sub('{}version=".*"{}'.format(*re.split('version=".*"', requirement)), requirement, xml_text)

        # if '@GALAXY_VERSION@' not in tags_to_update['tokens']:
        if tags_to_update.get('update_tool'):
            # update the version directly in the tool tag
            xml_text = re.sub('version="@TOOL_VERSION@\+galaxy.*"', 'version="@TOOL_VERSION@+galaxy0"', xml_text)

        f.seek(0)
        f.truncate()
        f.write(xml_text)

__all__ = (
    "autoupdate"
)
