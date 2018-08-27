"""gtdk: Galaxy training development kit."""

import collections
import json
import os
import re
import shutil

import oyaml as yaml
import requests

from planemo import templates
from planemo.bioblend import galaxy
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.io import info
from planemo.runnable import for_path


INDEX_FILE_TEMPLATE = """---
layout: topic
topic_name: {{ topic }}
---
"""


README_FILE_TEMPLATE = """
{{ topic }}
==========

Please refer to the [CONTRIBUTING.md](../../CONTRIBUTING.md) before adding or updating any material
"""


DOCKER_FILE_TEMPLATE = """
# Galaxy - {{ topic_title }}
#
# to build the docker image, go to root of training repo and
#    docker build -t {{ topic_name }} -f topics/{{ topic_name }}/docker/Dockerfile .
#
# to run image:
#    docker run -p "8080:80" -t {{ topic_name }}

FROM bgruening/galaxy-stable

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: {{ topic_title }}"

# prerequisites
RUN pip install ephemeris -U
ADD bin/galaxy-sleep.py /galaxy-sleep.py

# copy the tutorials directory for your topic
ADD topics/{{ topic_name }}/tutorials/ /tutorials/

# install everything for tutorials
ADD bin/docker-install-tutorials.sh /setup-tutorials.sh
ADD bin/mergeyaml.py /mergeyaml.py
RUN /setup-tutorials.sh
"""


INTRO_SLIDES_FILE_TEMPLATE = """---
layout: introduction_slides
logo: "GTN"

title: {{ title }}
type: {{ type }}
contributors:
- contributor
---

### How to fill the slide decks?

Please follow our
[tutorial to learn how to fill the slides]({{ '{{' }} site.baseurl {{ '}}' }}/topics/contributing/tutorials/create-new-tutorial-slides/slides.html)
"""

TUTO_SLIDES_TEMPLATE = """---
layout: tutorial_slides
logo: "GTN"

{{ metadata }}
---

### How to fill the slide decks?

Please follow our
[tutorial to learn how to fill the slides]({{ '{{' }} site.baseurl {{ '}}' }}/topics/contributing/tutorials/create-new-tutorial-slides/slides.html)
"""


TUTO_HAND_ON_TEMPLATE = """---
layout: tutorial_hands_on

{{ metadata }}
---

{{ body }}
"""


INPUT_FILE_TEMPLATE = """
>{{space}}- {{ '{%' }} icon {{icon}} {{ '%}' }} *"{{input_name}}"*: {{input_value}}
"""

INPUT_SECTION = """
>{{space}}- In *"{{section_label}}"*:
"""

INPUT_ADD_REPEAT = """
>{{space}}- Click on *"Insert {{repeat_label}}"*:
"""

INPUT_PARAM = """
>{{space}}- *"{{param_label}}"*: `{{param_value}}`
"""

HANDS_ON_TOOL_BOX_TEMPLATE = """
## Sub-step with **{{tool_name}}**

> ### {{ '{%' }} icon hands_on {{ '%}' }} Hands-on: Task description
>
> 1. **{{tool_name}}** {{ '{%' }} icon tool {{ '%}' }} with the following parameters:{{inputlist}}{{paramlist}}
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {{ '{%' }} icon comment {{ '%}' }} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {{ '{%' }} icon question {{ '%}' }} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {{ '{%' }} icon solution {{ '%}' }} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

"""

TUTO_HAND_ON_BODY_TEMPLATE = """
# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

**Please follow our
[tutorial to learn how to fill the Markdown]({{ '{{' }} site.baseurl {{ '}}' }}/topics/contributing/tutorials/\
create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.

Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {{ '{%' }} icon hands_on {{ '%}' }} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ zenodo_link }}) or from the shared data library
>
>    ```
>    {{ z_file_links }}
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {{ '{%' }} include snippets/import_via_link.md {{ '%}' }}
>    {{ '{%' }} include snippets/import_from_data_library.md {{ '%}' }}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {{ '{%' }} include snippets/change_datatype.md datatype="datatypes" {{ '%}' }}
>
> 5. Add to each database a tag corresponding to ...
>
>    {{ '{%' }} include snippets/add_tag.md {{ '%}' }}
>    
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {{ '{%' }} icon details {{ '%}' }} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:

{{ body }}

## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
"""

SPACE = '    '


def load_yaml(filepath):
    """Load the content of a YAML file to a dictionary."""
    with open(filepath, "r") as m_file:
        content = yaml.load(m_file)
    return content


def save_to_yaml(content, filepath):
    """Save a dictionary to a YAML file."""
    with open(filepath, 'w') as stream:
        yaml.safe_dump(content,
                       stream,
                       indent=2,
                       default_flow_style=False,
                       default_style='',
                       explicit_start=True,
                       encoding='utf-8',
                       allow_unicode=True)


def create_topic(kwds, topic_dir):
    """
    Create the skeleton of a new topic.

    1. create the folder and its structure
    2. update the index.md to match your topic's name
    3. fill the metadata
    4. add a symbolic link to the metadata.yaml from the metadata folder
    """
    # create the folder and its structure
    os.makedirs(topic_dir)
    img_folder = os.path.join(topic_dir, "images")
    os.makedirs(img_folder)
    tuto_folder = os.path.join(topic_dir, "tutorials")
    os.makedirs(tuto_folder)

    # create the index.md and add the topic name
    index_fp = os.path.join(topic_dir, "index.md")
    with open(index_fp, 'w') as index_f:
        index_f.write(
            templates.render(INDEX_FILE_TEMPLATE, **{'topic': kwds["topic_name"]}))

    # create the README file
    readme_fp = os.path.join(topic_dir, "README.md")
    with open(readme_fp, 'w') as readme_f:
        readme_f.write(
            templates.render(README_FILE_TEMPLATE, **{'topic': kwds["topic_title"]}))

    # create the metadata file
    metadata_fp = os.path.join(topic_dir, "metadata.yaml")
    metadata = collections.OrderedDict()
    metadata['name'] = kwds["topic_name"]
    metadata['type'] = kwds["topic_target"]
    metadata['title'] = kwds["topic_title"]
    metadata['summary'] = kwds["topic_summary"]
    metadata['requirements'] = []
    if metadata['type'] == 'use':
        req = collections.OrderedDict()
        req['title'] = "Galaxy introduction"
        req['type'] = "internal"
        req['link'] = "/introduction/"
        metadata['requirements'].append(req)
    metadata['docker_image'] = ""
    metadata['maintainers'] = ["maintainer"]
    if metadata['type'] == 'use':
        metadata['references'] = []
        ref = collections.OrderedDict()
        ref['authors'] = "authors et al"
        ref['title'] = "the title"
        ref['link'] = "link"
        ref['summary'] = "A short explanation of why this reference is useful"
        metadata['references'].append(ref)
    save_to_yaml(metadata, metadata_fp)

    # add a symbolic link to the metadata.yaml
    metadata_dir = "metadata"
    if not os.path.isdir(metadata_dir):
        os.makedirs(metadata_dir)
    os.chdir(metadata_dir)
    os.symlink(os.path.join("..", metadata_fp), "%s.yaml" % kwds["topic_name"])
    os.chdir("..")

    # create Dockerfile
    docker_folder = os.path.join(topic_dir, "docker")
    os.makedirs(docker_folder)
    dockerfile_fp = os.path.join(docker_folder, "Dockerfile")
    with open(dockerfile_fp, 'w') as dockerfile:
        dockerfile.write(
            templates.render(
                DOCKER_FILE_TEMPLATE,
                **{'topic_name': kwds["topic_name"], 'topic_title': kwds["topic_title"]}))

    # create empty introduction slides
    slides_folder = os.path.join(topic_dir, "slides")
    os.makedirs(slides_folder)
    intro_slide_fp = os.path.join(slides_folder, "introduction.html")
    with open(intro_slide_fp, 'w') as intro_slide_f:
        intro_slide_f.write(
            templates.render(
                INTRO_SLIDES_FILE_TEMPLATE,
                **{'title': "Introduction to %s" % kwds["topic_title"], 'type': "introduction"}))


def get_zenodo_record(zenodo_link):
    """Get the content of a Zenodo record."""
    # get the record in the Zenodo link
    if 'doi' in zenodo_link:
        z_record = zenodo_link.split('.')[-1]
    else:
        z_record = zenodo_link.split('/')[-1]
    # get JSON corresponding to the record from Zenodo API
    req = "https://zenodo.org/api/records/%s" % (z_record)
    r = requests.get(req)
    if r:
        req_res = r.json()
    else:
        info("The Zenodo link (%s) seems invalid" % (zenodo_link))
        req_res = {'files': []}
        z_record = None
    return(z_record, req_res)


def get_galaxy_datatype(z_ext, datatype_fp):
    """Get the Galaxy datatype corresponding to a Zenodo file type."""
    g_datatype = ''
    datatypes = load_yaml(datatype_fp)
    if z_ext in datatypes:
        g_datatype = datatypes[z_ext]
    if g_datatype == '':
        g_datatype = '# Please add a Galaxy datatype or update the shared/datatypes.yaml file'
    info("Get Galaxy datatypes: %s --> %s" % (z_ext, g_datatype))
    return g_datatype


def get_files_from_zenodo(z_link, datatype_fp):
    """
    Extract a list of URLs and dictionary describing the files from the JSON output of the Zenodo API.
    """
    z_record, req_res = get_zenodo_record(z_link)

    links = []
    if 'files' not in req_res:
        raise ValueError("No files in the Zenodo record")

    files = []
    for f in req_res['files']:
        file_dict = {'url': '', 'src': 'url', 'ext': '', 'info': z_link}
        if 'type' in f:
            file_dict['ext'] = get_galaxy_datatype(f['type'], datatype_fp)
        if 'links' not in f and 'self' not in f['links']:
            raise ValueError("No link for file %s" % f)
        file_dict['url'] = f['links']['self']
        links.append(f['links']['self'])
        files.append(file_dict)

    return (files, links, z_record)


def init_data_lib(data_lib_filepath):
    """Init the data library dictionary."""
    if os.path.exists(data_lib_filepath):
        data_lib = load_yaml(data_lib_filepath)
    else:
        data_lib = collections.OrderedDict()
    # set default information
    data_lib.setdefault('destination', collections.OrderedDict())
    data_lib['destination']['type'] = 'library'
    data_lib['destination']['name'] = 'GTN - Material'
    data_lib['destination']['description'] = 'Galaxy Training Network Material'
    data_lib['destination']['synopsis'] = 'Galaxy Training Network Material. See https://training.galaxyproject.org'
    data_lib.setdefault('items', [])
    data_lib.pop('libraries', None)
    return data_lib


def prepare_data_library(files, kwds, z_record, tuto_dir):
    """Fill or update the data library file."""
    data_lib_filepath = os.path.join(tuto_dir, "data-library.yaml")
    data_lib = init_data_lib(data_lib_filepath)
    # get topic or create new one
    topic = collections.OrderedDict()
    for item in data_lib['items']:
        if item['name'] == kwds['topic_title']:
            topic = item
    if not topic:
        data_lib['items'].append(topic)
        topic['name'] = kwds['topic_title']
        topic['description'] = kwds['topic_summary']
        topic['items'] = []
    # get tutorial or create new one
    tuto = collections.OrderedDict()
    for item in topic['items']:
        if item['name'] == kwds['tutorial_title']:
            tuto = item
    if not tuto:
        topic['items'].append(tuto)
        tuto['name'] = kwds['tutorial_title']
        tuto['items'] = []
    # get current data library and/or previous data library for the tutorial
    # remove the latest tag of any existing library
    # remove the any other existing library
    if z_record:
        current_data_lib = collections.OrderedDict()
        previous_data_lib = collections.OrderedDict()
        for item in tuto['items']:
            if item['name'] == "DOI: 10.5281/zenodo.%s" % z_record:
                current_data_lib = item
            elif item['description'] == 'latest':
                previous_data_lib = item
                previous_data_lib['description'] = ''
        if not current_data_lib:
            current_data_lib['name'] = "DOI: 10.5281/zenodo.%s" % z_record
            current_data_lib['description'] = 'latest'
            current_data_lib['items'] = []
        current_data_lib['items'] = files

        tuto['items'] = [current_data_lib]
        if previous_data_lib:
            tuto['items'].append(previous_data_lib)

    save_to_yaml(data_lib, data_lib_filepath)


def prepare_data_library_from_zenodo(kwds, tuto_dir):
    """Get the list of URLs of the files on Zenodo and fill the data library file."""
    links = []
    if not kwds['zenodo_link']:
        return links
    files, links, z_record = get_files_from_zenodo(kwds['zenodo_link'], kwds['datatypes'])
    prepare_data_library(files, kwds, z_record, tuto_dir)
    return links


def get_tool_input(tool_desc):
    """
    Get a dictionary with the tool descriptions.

    The labels are the tool parameter name and the value the description
    of the parameter extracted from the show_tool function of bioblend
    """
    tool_inp = collections.OrderedDict()
    for inp in tool_desc["inputs"]:
        tool_inp.setdefault(inp['name'], inp)
    return tool_inp


def get_wf_tool_description(wf, gi):
    """Get a dictionary with description of inputs of all tools in a workflow."""
    tools = {}
    for s in wf['steps']:
        step = wf['steps'][s]
        if not step['input_connections']:
            continue
        try:
            tool_desc = gi.tools.show_tool(step['tool_id'], io_details=True)
        except Exception:
            tool_desc = {'inputs': []}
        tools.setdefault(step['name'], get_tool_input(tool_desc))
    return tools


def get_wf_tool_from_local_galaxy(kwds, wf_filepath, ctx):
    """Server local Galaxy and get the workflow dictionary."""
    assert is_galaxy_engine(**kwds)
    runnable = for_path(wf_filepath)
    with engine_context(ctx, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([runnable]) as config:
            workflow_id = config.workflow_id(wf_filepath)
            wf = config.gi.workflows.export_workflow_dict(workflow_id)
            tools = get_wf_tool_description(wf, config.gi)
    return wf, tools


def get_wf_tools_from_running_galaxy(kwds):
    """Get the workflow dictionary from a running Galaxy instance with the workflow installed on it."""
    gi = galaxy.GalaxyInstance(kwds['galaxy_url'], key=kwds['galaxy_api_key'])
    wf = gi.workflows.export_workflow_dict(kwds['workflow_id'])
    tools = get_wf_tool_description(wf, gi)
    return wf, tools


def get_input_tool_name(step_id, steps):
    """Get the string with the name of the tool that generated an input."""
    inp_provenance = ''
    inp_prov_id = str(step_id)
    if inp_prov_id in steps:
        name = steps[inp_prov_id]['name']
        if 'Input dataset' in name:
            inp_provenance = "(%s)" % name
        else:
            inp_provenance = "(output of **%s** {%% icon tool %%})" % name
    return inp_provenance


def format_inputs(step_inputs, tp_desc, wf_steps, level):
    """Format the inputs of a step."""
    inputlist = ''
    for inp_n, inp in step_inputs.items():
        if inp_n != tp_desc['name']:
            continue
        inps = []
        if isinstance(inp, list):
            # multiple input (not collection)
            icon = 'param-files'
            for i in inp:
                inps.append('`%s` %s' % (
                    i['output_name'],
                    get_input_tool_name(i['id'], wf_steps)))
        else:
            # sinle input or collection
            inp_type = wf_steps[str(inp['id'])]['type']
            if 'collection' in inp_type:
                icon = 'param-collection'
            else:
                icon = 'param-file'
            inps = ['`%s` %s' % (
                inp['output_name'],
                get_input_tool_name(inp['id'], wf_steps))]
        inputlist += templates.render(INPUT_FILE_TEMPLATE, **{
            "icon": icon,
            "input_name": tp_desc['label'],
            "input_value": ', '.join(inps),
            "space": SPACE * level
        })
    return inputlist


def get_wf_step_inputs(step_inp):
    """Get the inputs from a workflow step and format them."""
    step_inputs = {}
    for inp_n, inp in step_inp.items():
        if '|' in inp_n:
            repeat_regex = '(?P<prefix>[^\|]*)_(?P<nb>\d+)\|(?P<suffix>.+).+'
            repeat_search = re.search(repeat_regex, inp_n)
            hier_regex = '(?P<prefix>[^\|]*)\|(?P<suffix>.+)'
            hier_regex = re.search(hier_regex, inp_n)
            if repeat_search and repeat_search.start(0) <= hier_regex.start(0):
                step_inputs.setdefault(repeat_search.group('prefix'), {})
                step_inputs[repeat_search.group('prefix')].setdefault(
                    repeat_search.group('nb'),
                    get_wf_step_inputs({hier_regex.group('suffix'): inp}))
            else:
                step_inputs.setdefault(hier_regex.group('prefix'), {})
                step_inputs[hier_regex.group('prefix')].update(
                    get_wf_step_inputs({hier_regex.group('suffix'): inp}))
        else:
            step_inputs.setdefault(inp_n, inp)
    return step_inputs


def json_load(string):
    """Transform a string into a dictionary."""
    if string is not None and ":" in string and '{' in string:
        return json.loads(string)
    else:
        return string


def get_lower_params(step_params, name):
    """Get the parameters from workflow that are below name in the hierarchy."""
    params = json_load(step_params)
    if isinstance(params, dict) and name in params:
        params = json_load(params[name])
    return params


def get_lower_inputs(step_inputs, name):
    """Get the inputs from workflow that are below name in the hierarchy."""
    inputs = {}
    if isinstance(step_inputs, dict) and name in step_inputs:
        inputs = step_inputs[name]
    else:
        inputs = step_inputs
    return inputs


def format_section_param_desc(step_params, step_inputs, tp_desc, level, wf_steps):
    """Format the description (label and value) for parameters in a section."""
    section_paramlist = ''
    # get sub params and inputs
    params = get_lower_params(step_params, tp_desc['name'])
    inputs = get_lower_inputs(step_inputs, tp_desc['name'])
    # get description of parameters in lower hierarchy
    sub_param_desc = get_param_desc(params, inputs, get_tool_input(tp_desc), level+1, wf_steps)
    if sub_param_desc != '':
        section_paramlist += templates.render(INPUT_SECTION, **{
            'space': SPACE * level, 
            'section_label': tp_desc['title']})
        section_paramlist += sub_param_desc
    return section_paramlist


def format_conditional_param_desc(step_params, step_inputs, tp_desc, level, wf_steps):
    """Format the description (label and value) for parameters in a conditional."""
    conditional_paramlist = ''
    # Get conditional parameter
    test_param = tp_desc['test_param']
    params = get_lower_params(step_params, tp_desc['name'])
    inputs = get_lower_inputs(step_inputs, tp_desc['name'])
    cond_param = get_lower_params(params, test_param['name'])
    print("-")
    print(cond_param)
    print("-")
    print(test_param)
    print("-")
    conditional_paramlist += format_param_desc(
        cond_param,
        step_inputs,
        test_param,
        level,
        wf_steps,
        force_default=True)
    # Get parameters in the when
    for case in tp_desc['cases']:
        if case['value'] == cond_param:
            if len(case['inputs']) > 0:
                conditional_paramlist += get_param_desc(
                    params,
                    inputs,
                    get_tool_input(case),
                    level+1,
                    wf_steps)
    return conditional_paramlist


def format_repeat_param_desc(step_params, step_inputs, tp_desc, level, wf_steps):
    """Format the description (label and value) for parameters in a repeat."""
    repeat_inp_desc = get_tool_input(tp_desc)
    params = get_lower_params(step_params, tp_desc['name'])
    inputs = get_lower_inputs(step_inputs, tp_desc['name'])
    repeat_paramlist = ''
    for r in range(len(params)):
        r_inputs = get_lower_inputs(inputs, str(r))
        r_params = params[r]
        paramlist_in_repeat = get_param_desc(r_params, r_inputs, repeat_inp_desc, level+2, wf_steps)
        if paramlist_in_repeat != '':
            # add first click
            repeat_paramlist += templates.render(INPUT_ADD_REPEAT, **{
                'space': SPACE * (level+1),
                'repeat_label': tp_desc['title']})
            # add description of parameters in the repeat
            repeat_paramlist += templates.render(INPUT_SECTION, **{
                'space': SPACE * (level+1),
                'section_label': "%s: %s" % (r+1, tp_desc['title'])})
            repeat_paramlist += paramlist_in_repeat
    if repeat_paramlist != '':
        repeat_paramlist = templates.render(INPUT_SECTION, **{
            'space': SPACE * level,
            'section_label': tp_desc['title']}) + repeat_paramlist
    return repeat_paramlist


def get_param_value(step_params, tp_desc, force_default=False):
    """Get value of a 'simple' parameter if different from the default value, None otherwise."""
    param_value = ''
    if isinstance(step_params, str) and '"' in step_params:
        step_params = step_params.replace('"', '')
    if tp_desc['value'] == step_params and not force_default:
        param_value = None
    elif tp_desc['type'] == 'boolean':
        if bool(tp_desc['value']) == step_params:
            param_value = None
        else:
            param_value = 'Yes' if step_params else 'No'
    elif tp_desc['type'] == 'select':
        param_value = ''
        for opt in tp_desc['options']:
            if opt[1] == step_params:
                param_value = opt[0]
    elif tp_desc['type'] == 'data_column':
        param_value = "c%s" % step_params
    else:
        param_value = step_params
    return param_value


def format_param_desc(step_params, step_inputs, tp_desc, level, wf_steps, force_default=False):
    """Format the parameter description (label and value) given the type of parameter."""
    paramlist = ''
    if 'type' not in tp_desc:
        raise ValueError("No type for the paramater %s" % tp_desc['name'])
    if tp_desc['type'] == 'data' or tp_desc['type'] == 'data_collection':
        paramlist += format_inputs(step_inputs, tp_desc, wf_steps, level)
    elif tp_desc['type'] == 'section':
        paramlist += format_section_param_desc(step_params, step_inputs, tp_desc, level, wf_steps)
    elif tp_desc['type'] == 'conditional':
        paramlist += format_conditional_param_desc(step_params, step_inputs, tp_desc, level, wf_steps)
    elif tp_desc['type'] == 'repeat':
        paramlist += format_repeat_param_desc(step_params, step_inputs, tp_desc, level, wf_steps)
    else:
        param_value = get_param_value(step_params, tp_desc, force_default)
        if param_value is not None:
            paramlist += templates.render(INPUT_PARAM, **{
                'space': SPACE * level,
                'param_label': tp_desc['label'],
                'param_value': param_value})
    return paramlist


def get_param_desc(step_params, step_inputs, tp_desc, level, wf_steps, should_be_there=False):
    """Parse the parameters of the tool and return a formatted list with the values set in the workflow."""
    paramlist = ''
    for n, tp_d in tp_desc.items():
        if n not in step_params:
            if not should_be_there:
                info("%s not in workflow" % n)
            else:
                raise ValueError("%s not in workflow" % n)
        else:
            step_param = get_lower_params(step_params, n)
            if step_param is None:
                continue
            paramlist += format_param_desc(step_param, step_inputs, tp_d, level, wf_steps)
    return paramlist


def get_handson_box(step, steps, tools):
    """Get the string for an hands-on box based on a step in a workflow."""
    # get input (if none: input step)
    step_inputs = get_wf_step_inputs(step['input_connections'])
    if not step_inputs:
        return ''
    # get params
    step_params = json.loads(step['tool_state'])
    # get tool
    tool_name = step['name']
    tp_desc = tools[tool_name]
    # get formatted param description
    paramlist = get_param_desc(step_params, step_inputs, tp_desc, 1, steps, should_be_there=True)
    context = {"tool_name": tool_name, "paramlist": paramlist}
    return templates.render(HANDS_ON_TOOL_BOX_TEMPLATE, **context)


def init_tuto_metadata(kwds):
    """Init tutorial metadata"""
    metadata = collections.OrderedDict()
    metadata['title'] = kwds["tutorial_title"]
    metadata['zenodo_link'] = kwds["zenodo_link"] if kwds["zenodo_link"] else ''
    metadata['questions'] = [
            "Which biological questions are addressed by the tutorial?",
            "Which bioinformatics techniques is important to know for this type of data?"]
    metadata['objectives'] = [
            "The learning objectives are the goals of the tutorial",
            "They will be informed by your audience and will communicate to them and to yourself what you should focus on during the course",
            "They are single sentence describing what a learner will be able to do once they have done the tutorial",
            "You can use the Bloom's Taxonomy to write effective learning objectives"]
    metadata['time'] = "3H"
    metadata['key_points'] = [
            "The take-home messages",
            "They will appear at the end of the tutorial"]
    metadata['contributors'] = ["contributor1", "contributor2"]
    return metadata


def format_tuto_metadata(metadata):
    """Return the string corresponding to the tutorial metadata"""
    return yaml.safe_dump(metadata,
                           indent=2,
                           default_flow_style=False,
                           default_style='',
                           explicit_start=False)
        

def write_hands_on_tutorial(metadata, body, tuto_dir):
    """Write the tutorial hands-on"""
    m_str = format_tuto_metadata(metadata)
    template = templates.render(TUTO_HAND_ON_TEMPLATE, **{
        "metadata": m_str,
        "body": body
    })

    md_path = os.path.join(tuto_dir, "tutorial.md")
    with open(md_path, 'w') as md:
        md.write(template)


def get_tuto_body(z_file_links, body = None):
    """Get the body for a tutorial"""
    if body is None:
        body = templates.render(HANDS_ON_TOOL_BOX_TEMPLATE, **{
            'tool_name': "My Tool",
            'inputlist': templates.render(INPUT_FILE_TEMPLATE, **{
                'space': 1*SPACE,
                'icon': 'param-file',
                'input_name': 'Input file',
                'input_value': 'File'
                }),
            'paramlist': templates.render(INPUT_PARAM, **{
                'space': 1*SPACE,
                'param_label': 'Parameter',
                'param_value': 'a value'
                })
        })
    return templates.render(TUTO_HAND_ON_BODY_TEMPLATE, **{
        "z_file_links": "\n>    ".join(z_file_links),
        "body": body})


def create_hands_on_tutorial_from_workflow(kwds, z_file_links, tuto_dir, ctx, metadata=None):
    """Create tutorial structure from the workflow file."""
    # load workflow
    if kwds['workflow_id']:
        if not kwds['galaxy_url']:
            raise ValueError("No Galaxy URL given")
        if not kwds['galaxy_api_key']:
            raise ValueError("No API key to access Galaxy given")
        wf, tools = get_wf_tools_from_running_galaxy(kwds)
    else:
        wf, tools = get_wf_tool_from_local_galaxy(kwds, kwds["workflow"], ctx)

    # get hands-on body from the workflow
    body = ''
    for step_id in range(len(wf['steps'].keys())):
        step = wf['steps'][str(step_id)]
        if not step['tool_state']:
            continue
        body += get_handson_box(step, wf['steps'], tools)
    body = get_tuto_body(z_file_links, body)

    # write in the tutorial file with the metadata on the top
    if not metadata:
        metadata = init_tuto_metadata(kwds)
    write_hands_on_tutorial(metadata, body, tuto_dir)


def add_workflow_file(kwds, tuto_dir):
    """Copy or extract workflow file and add it to the tutorial directory"""
    wf_dir = os.path.join(tuto_dir, "workflows")
    # copy / extract workflow
    wf_filepath = os.path.join(wf_dir, "init_workflow.ga")
    if kwds["workflow"]:
        shutil.copy(kwds["workflow"], wf_filepath)
    else:
        gi = galaxy.GalaxyInstance(kwds['galaxy_url'], key=kwds['galaxy_api_key'])
        gi.workflows.export_workflow_to_local_path(kwds['workflow_id'],
                                                   wf_filepath,
                                                   use_default_filename=False)
    # remove empty workflow file if there
    empty_wf_filepath = os.path.join(wf_dir, "empty_workflow.ga")
    if os.path.exists(empty_wf_filepath):
        os.remove(empty_wf_filepath)


def create_tutorial(kwds, tuto_dir, ctx):
    """Create the skeleton of a new tutorial."""
    # create tuto folder and empty files
    os.makedirs(tuto_dir)
    tour_folder = os.path.join(tuto_dir, "tours")
    os.makedirs(tour_folder)
    workflow_folder = os.path.join(tuto_dir, "workflows")
    os.makedirs(workflow_folder)

    metadata = init_tuto_metadata(kwds)

    # extract the data library from Zenodo and the links for the tutorial
    z_file_links = ''
    if kwds["zenodo_link"]:
        info("Create the data library from Zenodo")
        z_file_links = prepare_data_library_from_zenodo(kwds, tuto_dir)

    # create tutorial skeleton from workflow and copy workflow file
    if kwds["hands_on"]:
        if kwds["workflow"] or kwds['workflow_id']:
            info("Create tutorial skeleton from workflow")
            create_hands_on_tutorial_from_workflow(kwds, z_file_links, tuto_dir, ctx)
            add_workflow_file(kwds, tuto_dir)
        else:
            body = get_tuto_body(z_file_links)
            print(body)
            write_hands_on_tutorial(metadata, body, tuto_dir)

    # create slide skeleton
    if kwds["slides"]:
        slide_path = os.path.join(tuto_dir, 'slides.html')
        m_str = format_tuto_metadata(metadata)
        with open(slide_path, 'w') as slide_f:
            slide_f.write(
                templates.render(TUTO_SLIDES_TEMPLATE, **{"metadata": m_str}))


def init(ctx, kwds):
    """Create/update a topic/tutorial"""
    topic_dir = os.path.join("topics", kwds['topic_name'])
    if not os.path.isdir(topic_dir):
        info("The topic %s does not exist. It will be created" % kwds['topic_name'])
        create_topic(kwds, topic_dir)

    if not kwds['tutorial_name']:
        if kwds["slides"]:
            raise Exception("A tutorial name is needed to create the skeleton of a tutorial slide deck")
        if kwds['workflow'] or kwds['workflow_id']:
            raise Exception("A tutorial name is needed to create the skeleton of the tutorial from a workflow")
        if kwds['zenodo_link']:
            raise Exception("A tutorial name is needed to add Zenodo information")
    else:
        tuto_dir = os.path.join(topic_dir, "tutorials", kwds['tutorial_name'])
        if not os.path.isdir(tuto_dir):
            info("The tutorial %s in topic %s does not exist. It will be created." % (kwds['tutorial_name'], kwds['topic_name']))
            create_tutorial(kwds, tuto_dir, ctx)


def get_tuto_info(tuto_dir):
    """Extract the metadata front matter on the top of the tutorial file and its body"""
    tuto_fp = os.path.join(tuto_dir, "tutorial.md")
    with open(tuto_fp, "r") as tuto_f:
        tuto_content = tuto_f.read()

    regex = '^---\n(?P<metadata>[\s\S]*)\n---(?P<body>[\s\S]*)'
    tuto_split_regex = re.search(regex, tuto_content)
    if not tuto_split_regex:
        raise Exception("No metadata found at the top of the tutorial")

    metadata = yaml.load(tuto_split_regex.group("metadata"))
    body = tuto_split_regex.group("body")

    return metadata, body


def check_topic_tuto_exist(kwds):
    """Check that the topic and tutorial are already there."""
    topic_dir = os.path.join("topics", kwds['topic_name'])
    if not os.path.isdir(topic_dir):
        raise Exception("The topic %s does not exists. It should be created" % kwds['topic_name'])

    tuto_dir = os.path.join(topic_dir, "tutorials", kwds['tutorial_name'])
    if not os.path.isdir(tuto_dir):
        raise Exception("The tutorial %s does not exists. It should be created" % kwds['tutorial_name'])

    return topic_dir, tuto_dir


def fill_data_library(ctx, kwds):
    """Fill a data library for a tutorial."""
    topic_dir, tuto_dir = check_topic_tuto_exist(kwds)
    metadata, body = get_tuto_info(tuto_dir)

    # get the zenodo link
    z_link = ''
    if 'zenodo_link' in metadata and metadata['zenodo_link'] != '':
        if kwds['zenodo_link']:
            info("The data library and the metadata will be updated with the new Zenodo link")
            z_link = kwds['zenodo_link']
            metadata['zenodo_link'] = z_link
        else:
            info("The data library will be extracted using the Zenodo link in the metadata")
            z_link = metadata['zenodo_link']
    elif kwds['zenodo_link']:
        info("The data library will be created and the metadata will be filled with the new Zenodo link")
        z_link = kwds['zenodo_link']
        metadata['zenodo_link'] = z_link

    if z_link == '' or z_link is None:
        raise Exception("A Zenodo link should be provided either in the metadata file or as argument of the command")

    # get the topic metadata
    topic_metadata_fp = os.path.join(topic_dir, "metadata.yaml")
    topic_metadata = load_yaml(topic_metadata_fp)

    # extract the data library from Zenodo
    topic_kwds = {
        'topic_title': topic_metadata['title'],
        'topic_summary': topic_metadata['summary'],
        'tutorial_title': metadata['title'],
        'zenodo_link': z_link,
        'datatypes': kwds['datatypes']
    }
    prepare_data_library_from_zenodo(topic_kwds, tuto_dir)

    # update the metadata
    write_hands_on_tutorial(metadata, body, tuto_dir)


def generate_tuto_from_wf(ctx, kwds):
    """Generate the skeleton of a tutorial from a workflow."""
    if kwds["workflow"] or kwds['workflow_id']:
        topic_dir, tuto_dir = check_topic_tuto_exist(kwds)
        metadata, body = get_tuto_info(tuto_dir)
        info("Create tutorial skeleton from workflow")
        create_hands_on_tutorial_from_workflow(kwds, [], tuto_dir, ctx, metadata)
        add_workflow_file(kwds, tuto_dir)
    else:
        exc = "A path to a local workflow or the id of a workflow on a running Galaxy instance should be provided"
        raise Exception(exc)
