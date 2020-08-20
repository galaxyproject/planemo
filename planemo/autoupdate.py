"""Autoupdate older conda dependencies in the requirements section."""
from __future__ import absolute_import

import collections
import re
import xml.etree.ElementTree as ET

from galaxy.tool_util.deps import conda_util

import planemo.conda
from planemo.io import error, info


def find_macros(xml_tree):
    macros = []
    for macro_import in xml_tree.iter("import"):
        macros.append(macro_import.text)
    if macros:
        return macros
    else:
        return None


def get_requirements(xml_tree):
    requirements = {}
    main_req = None
    for requirement in xml_tree.iter("requirement"):
        if requirement.attrib.get('version') == '@TOOL_VERSION@':
            main_req = requirement.text
        else:
            requirements[requirement.text] = requirement.attrib.get('version')
    return requirements, main_req


def get_tokens(xml_tree):
    tokens = {}
    for token in xml_tree.iter("token"):
        tokens[token.attrib['name']] = token.text
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


def update_xml(tool_path, xml_tree, tags_to_update):
    """
    Write modified XML to tool_path
    """

    print(tool_path)
    print(tags_to_update)

    def update_token(xml_text, token_name, token_value):
        return re.sub('{}>.*<{}'.format(*re.split('>.*<', token)), token, xml_text)

    def update_requirement(xml_text, requirement_name, requirement_value):
        return

    with open(tool_path, 'r+') as f:
        xml_text = f.read()

        for tag_to_update in tags_to_update:
            if tag_to_update['type'] = 'token':
                xml_text = update_token(xml_text, tag_to_update['key'], tag_to_update['value'])
            if tag_to_update['type'] = 'requirement':
                update_requirement(xml_text, tag_to_update['key'], tag_to_update['value'])




    # print(tags_to_update)
    # with open(tool_path, 'r+') as f:
    #     xml_text = f.read()
    #     for token in tags_to_update['tokens']:
    #         xml_text = re.sub('{}>.*<{}'.format(*re.split('>.*<', token)), token, xml_text)

    #     for requirement in tags_to_update['requirements']:
    #         xml_text = re.sub('{}version=".*"{}'.format(*re.split('version=".*"', requirement)), requirement, xml_text)

    #     # if '@GALAXY_VERSION@' not in tags_to_update['tokens']:
    #     if tags_to_update.get('update_tool'):
    #         # update the version directly in the tool tag
    #         xml_text = re.sub(r'version="@TOOL_VERSION@\+galaxy.*?"', 'version="@TOOL_VERSION@+galaxy0"', xml_text)

    #     f.seek(0)
    #     f.truncate()
    #     f.write(xml_text)

def autoupdate(ctx, tool_path, **kwds):
    """
    Autoupdate an XML file
    """
    # create a dict of all files that need editing - wrapper plus macros
    xml_files = {tool_path: ET.parse(tool_path)}
    for macro in find_macros(xml_files[tool_path]):
        macro_path = '/'.join(tool_path.split('/')[:-1] + [macro])
        xml_files[macro_path] = ET.parse(macro_path)

    # create a dict of requirements
    requirements = {} #collections.defaultdict(dict)
    main_req = None
    for k, v in xml_files.items():
        r, mr = get_requirements(v)
        requirements[k] = r
        if mr:
            if main_req:
                error('Multiple requirements use the token @TOOL_VERSION@!')
            main_req = (mr, k)  #{mr: k} 
    if not main_req:
        error('No requirement uses the token @TOOL_VERSION@!')

    # create a dict of tokens
    tokens = {}
    need_to_update = False
    xml_to_update = collections.defaultdict(list)
    for k, v in xml_files.items():
        tokens[k] = get_tokens(v)
        if '@TOOL_VERSION@' in tokens[k]:
            check_main_req = check_conda(main_req[0], tokens[k]['@TOOL_VERSION@'], ctx, **kwds)
            if check_main_req:
                xml_to_update[k].append({'type': 'token', 'key': '@TOOL_VERSION@', 'value': check_main_req})
                need_to_update = True
    if not need_to_update:
        info("No updates required or made to {}.".format(tool_path))
        return  # end here if no update needed

    # check all requirements
    for k, v in requirements.items():
        for req in v:
            req_check = check_conda(req, v[req], ctx, **kwds)
            if req_check:
                xml_to_update[k].append({'type': 'requirement', 'tag': req, 'value': req_check})

    # check all tokens
    for k, v in tokens.items():
        if '@GALAXY_VERSION' in v:
            xml_to_update[k].append({'type': 'token', 'key': '@GALAXY_VERSION@', 'value': 0})

    # finally, update each file separately
    for k, v in xml_files.items():
        update_xml(k, v, xml_to_update[k])


    print(tokens)
    print(requirements)
    print(xml_to_update)
        # requirements[k] = r
        # if mr:
        #     if main_req:
        #         error('Multiple tools use the token @TOOL_VERSION@!')
        #     main_req = {mr: k}  # (mr, k)

    # tokens = get_tokens(xml_files[tool_path]))
    

    # print(requirements)
    # print(main_req)









# def get_tokens():
#     return ...

# def conda_check():
#     return ...

# def update():
#     return