"""Autoupdate older conda dependencies in the requirements section."""
from __future__ import absolute_import

import collections
import re
import xml.etree.ElementTree as ET

from galaxy.tool_util.deps import conda_util

import planemo.conda
from planemo.io import error, info


def find_macros(xml_tree):
    """
    Get macros from the XML tree
    """
    macros = []
    for macro_import in xml_tree.iter("import"):
        macros.append(macro_import.text)
    return macros


def get_requirements(xml_tree):
    """
    Get requirements from the XML tree
    """
    requirements = {}
    main_req = None
    for requirement in xml_tree.iter("requirement"):
        if requirement.attrib.get('version') == '@TOOL_VERSION@':
            main_req = requirement.text
        else:
            requirements[requirement.text] = {'tag': ET.tostring(requirement, encoding="unicode").strip(), 'text': requirement.attrib.get('version')}
    return requirements, main_req


def get_tokens(xml_tree):
    """
    Get tokens from the XML tree
    """
    tokens = {}
    for token in xml_tree.iter("token"):
        tokens[token.attrib['name']] = {'tag': ET.tostring(token, encoding="unicode").strip(), 'text': token.text}
    return tokens


def check_conda(tool_name, tool_version, ctx, **kwds):
    """
    Update a dict with current conda versions for tool requirements
    """
    conda_context = planemo.conda.build_conda_context(ctx, **kwds)
    target = planemo.conda.conda_util.CondaTarget(tool_name)
    search_results = conda_util.best_search_result(target, conda_context=conda_context)
    if search_results[0]['version'] == tool_version:
        return None
    else:
        return search_results[0]['version']


def update_xml(tool_path, xml_tree, tags_to_update, wrapper_version_token, is_macro=False):
    """
    Write modified XML to tool_path
    """
    def update_token(xml_text, tag, token_value):
        new_tag = '>{}<'.format(token_value).join(re.split('>.*<', tag))
        return re.sub(tag, new_tag, xml_text)

    def update_requirement(xml_text, tag, requirement_value):
        new_tag = 'version="{}"'.format(requirement_value).join(re.split('version=".*"', tag))
        return re.sub(tag, new_tag, xml_text)

    with open(tool_path, 'r+') as f:
        xml_text = f.read()
        for tag_to_update in tags_to_update:
            if tag_to_update['type'] == 'token':
                xml_text = update_token(xml_text, tag_to_update['tag'], tag_to_update['value'])
            if tag_to_update['type'] == 'requirement':
                xml_text = update_requirement(xml_text, tag_to_update['tag'], tag_to_update['value'])
        if not wrapper_version_token and not is_macro:
            # i.e. @GALAXY_VERSION@ not specified so update the version directly in the tool tag
            tool_tag = re.sub('version="@TOOL_VERSION@.*?"', 'version="@TOOL_VERSION@+galaxy0"',
                              re.findall('<tool .*version="@TOOL_VERSION@.*">', xml_text)[0])
            xml_text = re.sub('<tool .*version="@TOOL_VERSION@.*">', tool_tag, xml_text)
        f.seek(0)
        f.truncate()
        f.write(xml_text)


def autoupdate(ctx, tool_path, modified_files=set(), **kwds):
    """
    Autoupdate an XML file
    """
    # create a dict of all files that need editing - wrapper plus macros
    xml_files = {tool_path: ET.parse(tool_path)}

    # get name of token which defines the wrapper version; if just an integer, None
    versions = xml_files[tool_path].getroot().attrib.get('version')
    if versions:
        versions = versions.split('+galaxy')
        if versions[0] != '@TOOL_VERSION@':
            error('Tool version does not contain @TOOL_VERSION@ as required by autoupdate.')
            return
        elif len(versions) == 1:
            wrapper_version_token = None
        else:
            if versions[1][0] == versions[1][-1] == '@':
                wrapper_version_token = versions[1]
            else:
                wrapper_version_token = None
    else:
        wrapper_version_token = None

    # add macros to xml_files
    for macro in find_macros(xml_files[tool_path]):
        macro_path = '/'.join(tool_path.split('/')[:-1] + [macro])
        xml_files[macro_path] = ET.parse(macro_path)

    # create a dict of requirements
    requirements = {}
    main_req = None
    for k, v in xml_files.items():
        file_reqs, file_main_req = get_requirements(v)
        requirements[k] = file_reqs
        if file_main_req:
            if main_req:
                error('Multiple requirements use the token @TOOL_VERSION@!')
            main_req = (file_main_req, k)
    if not main_req:
        error('No requirement uses the token @TOOL_VERSION@!')

    # create a dict of tokens
    tokens = {}
    need_to_update = False
    xml_to_update = collections.defaultdict(list)
    for k, v in xml_files.items():
        tokens[k] = get_tokens(v)
        # check if it is @TOOL_VERSION@ and if so do check_conda
        if '@TOOL_VERSION@' in tokens[k]:
            check_main_req = check_conda(main_req[0], tokens[k]['@TOOL_VERSION@']['text'], ctx, **kwds)
            if check_main_req:
                xml_to_update[k].append({'type': 'token', 'tag': tokens[k]['@TOOL_VERSION@']['tag'], 'value': check_main_req})
                need_to_update = True
    if not need_to_update and not (modified_files & set(xml_files)):
        info("No updates required or made to {}.".format(tool_path))
        return  # end here if no update needed

    if kwds.get('dry_run'):
        error("Update required to {}! Tool main requirement has version {}, newest conda version is {}".format(
              tool_path, tokens[k]['@TOOL_VERSION@']['text'], check_main_req))
        return

    # check all requirements
    for k, v in requirements.items():
        for req in v:
            req_check = check_conda(req, v[req]['text'], ctx, **kwds)
            if req_check:
                xml_to_update[k].append({'type': 'requirement', 'tag': v[req]['tag'], 'value': req_check})

    # check all tokens, if wrapper_version_token exists
    if wrapper_version_token:
        for k, v in tokens.items():
            if wrapper_version_token in v:
                xml_to_update[k].append({'type': 'token', 'tag': v[wrapper_version_token]['tag'], 'value': 0})

    # finally, update each file separately
    for k, v in xml_files.items():
        update_xml(k, v, xml_to_update[k], wrapper_version_token, is_macro=(k != tool_path))

    info("Tool {} updated.".format(tool_path))
    return list(xml_files)
