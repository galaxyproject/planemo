""" gtdk: Galaxy training development kit """

import collections
import json
import os
import requests
import shutil
import time
import oyaml as yaml


from planemo import templates
from planemo.bioblend import galaxy
from planemo.io import info
from planemo.runnable import for_path
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)


INPUT_FILE_TEMPLATE = """
>   - {{ '{%' }} icon {{icon}} {{ '%}' }} *"{{input_name}}"*: {{input_value}}
"""

INPUT_SECTION = """
>{{space}}- In *"{{section_label}}"*:
"""

INPUT_ADD_REPEAT = """
>{{space}}- Click on *"Insert {{section_label}}"*:
"""

INPUT_PARAM = """
>{{space}}- *"{{param_label}}"*: `{{param_value}}`
"""

HANDS_ON_TOOL_BOX_TEMPLATE = """
# Title of the section usually corresponding to a big step

Description of the step: some background and some theory. Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part. 

<-- Consider adding a detail box to expand the theory -->

> ### {{ '{%' }} icon details {{ '%}' }} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
> 
{: .details}

> ### {{ '{%' }} icon hands_on {{ '%}' }} Hands-on: TODO: task description
>
> 1. **{{tool_name}}** {{ '{%' }} icon tool {{ '%}' }} with the following parameters:{{inputlist}}{{paramlist}}
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
>   <-- Consider adding a comment or tip box -->
>
>   > ### {{ '{%' }} icon comment {{ '%}' }}} Comment
>   >
>   > A comment about the tool or something else. This box can also be in the main text
>   {: .comment}
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->

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


TUTORIAL_TEMPLATE = """---
layout: tutorial_hands_on
topic_name: {{ topic_name }}
tutorial_name: {{ tutorial_name }}
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the 
tutorial (the questions and the objectives). It is nice also to have a 
scheme to sum up the pipeline used during the tutorial. The idea is to 
give to trainees insight into the content of the tutorial and the (theoretical 
and technical) key concepts they will learn.

**Please follow our 
[tutorial to learn how to fill the Markdown]({{ '{{' }} site.baseurl {{ '}}' }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

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

Anywhere you find the word `TODO`, there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {{ '{%' }} icon hands_on {{ '%}' }} Hands-on: Data upload
>
> 1. Import the following files from [Zenodo]({{ zenodo_link }}) or from a data
>    library named `TODO` if available (ask your instructor)
>
>    ```
>    TODO: add the files by the ones on Zenodo here (if not added)
>    TODO: remove the useless files (if added)
>    TODO: so that they can easily be copy-pasted into Galaxy's upload dialog
>    {{ z_file_links }}
>    ```
>
>    > ### {{ '{%' }} icon tip {{ '%}' }} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    >
>    > By default, Galaxy uses the url as the name, so please rename them to something more pleasing.
>    {: .tip}
>
>    > ### {{ '{%' }} icon tip {{ '%}' }} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Analyses of metagenomics data"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
{: .hands_on}

{{ body }}

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
"""


def load_yaml(filepath):
    """Load the content of a YAML file to a dictionary"""
    with open(filepath, "r") as m_file:
        content = yaml.load(m_file)
    return content


def save_to_yaml(content, filepath):
    """Save a dictionary to a YAML file"""
    with open(filepath, 'w') as stream:
        yaml.dump(content,
                  stream,
                  indent=2,
                  default_flow_style=False,
                  default_style='',
                  explicit_start=True)


def get_template_dir():
    """Check and return the templates directory"""
    template_dir = "templates"
    if not os.path.isdir(template_dir):
        raise Exception("This script needs to be run in the training material repository")
    return template_dir


def change_topic_name(topic_name, filepath):
    """Change the topic name in the top metadata of a file"""
    with open(filepath, "r") as in_f:
        content = in_f.read()

    content = content.replace("your_topic", topic_name)
    content = content.replace("your_tutorial_name", "tutorial1")

    with open(filepath, 'w') as out_f:
        out_f.write(content)


def create_topic(kwds, topic_dir, template_dir):
    """Create the skeleton of a new topic:

    1. copy templates
    2. update the index.md to match your topic's name
    3. fill the metadata
    4. add a symbolic link to the metadata.yaml from the metadata folder
    """
    # copy templates
    shutil.copytree(template_dir, topic_dir)

    # update the index.md to match your topic's name
    index_path = os.path.join(topic_dir, "index.md")
    change_topic_name(kwds["topic_name"], index_path)

    # update the metadata file
    metadata_path = os.path.join(topic_dir, "metadata.yaml")

    metadata = load_yaml(metadata_path)
    metadata['name'] = kwds["topic_name"]
    metadata['title'] = kwds["topic_title"]
    metadata['type'] = kwds["topic_target"]
    metadata['summary'] = kwds["topic_summary"]

    save_to_yaml(metadata, metadata_path)

    # update the metadata in top of tutorial.md and slides.html
    tuto_path = os.path.join(topic_dir, "tutorials", "tutorial1")
    hand_on_path = os.path.join(tuto_path, "tutorial.md")
    change_topic_name(kwds["topic_name"], hand_on_path)
    slides_path = os.path.join(tuto_path, "slides.html")
    change_topic_name(kwds["topic_name"], slides_path)

    # add a symbolic link to the metadata.yaml
    os.chdir("metadata")
    os.symlink(os.path.join("..", metadata_path), "%s.yaml" % kwds["topic_name"])
    os.chdir("..")


def update_tuto_file(filepath, keep, topic_name, tutorial_name):
    """Update or delete a tutorial (hands-on or slide) file"""
    if keep:
        with open(filepath, "r") as in_f:
            content = in_f.read()

        content = content.replace("your_topic", topic_name)
        content = content.replace("your_tutorial_name", tutorial_name)

        with open(filepath, 'w') as out_f:
            out_f.write(content)

    elif os.path.isfile(filepath):
        os.remove(filepath)


def update_tutorial(kwds, tuto_dir, topic_dir):
    """Update the metadata information of a tutorial"""
    # update the metadata file to add the new tutorial
    metadata_path = os.path.join(topic_dir, "metadata.yaml")

    metadata = load_yaml(metadata_path)
    found = False
    for mat in metadata["material"]:
        if mat["name"] == kwds["tutorial_name"]:
            mat["name"] = kwds["tutorial_name"]
            mat["title"] = kwds["tutorial_title"]
            mat["hands_on"] = kwds["hands_on"]
            mat["slides"] = kwds["slides"]
            mat["workflows"] = True if kwds["workflow"] else False
            mat["zenodo_link"] = kwds["zenodo"] if kwds["zenodo"] else ''
            found = True
        elif mat["name"] == "tutorial1":
            metadata["material"].remove(mat)

    if not found:
        new_mat = collections.OrderedDict()
        new_mat["title"] = kwds["tutorial_title"]
        new_mat["name"] = kwds["tutorial_name"]
        new_mat["type"] = 'tutorial'
        new_mat["zenodo_link"] = kwds["zenodo"] if kwds["zenodo"] else ''
        new_mat["hands_on"] = kwds["hands_on"]
        new_mat["slides"] = kwds["slides"]
        new_mat["workflows"] = True if kwds["workflow"] else False
        new_mat["galaxy_tour"] = False
        new_mat["questions"] = ['', '']
        new_mat["objectives"] = ['', '']
        new_mat["time_estimation"] = '1d/3h/6h'
        new_mat["key_points"] = ['', '']
        new_mat["contributors"] = ['contributor1', 'contributor2']
        metadata["material"].append(new_mat)

    save_to_yaml(metadata, metadata_path)

    # update the metadata in top of tutorial.md or remove it if not needed
    hand_on_path = os.path.join(tuto_dir, "tutorial.md")
    update_tuto_file(hand_on_path, kwds["hands_on"], kwds["topic_name"], kwds["tutorial_name"])

    # update the metadata in top of slides.md or remove it if not needed
    slides_path = os.path.join(tuto_dir, "slides.html")
    update_tuto_file(slides_path, kwds["slides"], kwds["topic_name"], kwds["tutorial_name"])


def get_zenodo_record(zenodo_link):
    """Get the content of a Zenodo record"""
    # get the record in the Zenodo link
    if 'doi' in zenodo_link:
        z_record = zenodo_link.split('.')[-1]
    else:
        z_record = zenodo_link.split('/')[-1]

    # get JSON corresponding to the record from Zenodo API
    req = "https://zenodo.org/api/records/%s" % (z_record)
    r = requests.get(req)
    r.raise_for_status()
    req_res = r.json()

    return(z_record, req_res)


def prepare_data_library(files, kwds, z_record, tuto_dir):
    """Fill or update the data library file"""
    data_lib_filepath = os.path.join(tuto_dir, "data-library.yaml")

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


def get_galaxy_datatype(z_ext, kwds):
    """Get the Galaxy datatype corresponding to a Zenodo file type"""
    g_datatype = '' 
    datatypes = load_yaml(kwds['datatypes'])
    if z_ext in datatypes:
        g_datatype = datatypes[z_ext]
    if g_datatype == '':
        g_datatype = '# Please add a Galaxy datatype or update the shared/datatypes.yaml file' 
    info("Get Galaxy datatypes: %s --> %s" %(z_ext, g_datatype))
    return g_datatype


def extract_from_zenodo(kwds, tuto_dir):
    """Get the list of URLs of the files on Zenodo and fill the data library file"""
    links = []
    if not kwds['zenodo']:
        return links

    z_record, req_res = get_zenodo_record(kwds['zenodo'])

    # extract the URLs from the JSON
    if 'files' not in req_res:
        raise ValueError("No files in the Zenodo record")

    files = []
    for f in req_res['files']:
        file_dict = {'url':'', 'src': 'url', 'ext': '', 'info': kwds['zenodo']}
        if 'type' in f:
            file_dict['ext'] = get_galaxy_datatype(f['type'], kwds)
        if 'links' not in f and 'self' not in f['links']:
            raise ValueError("No link for file %s" % f)
        file_dict['url'] = f['links']['self']
        links.append(f['links']['self'])
        files.append(file_dict)

    # prepare the data library dictionary
    prepare_data_library(files, kwds, z_record, tuto_dir)
    
    return links


def get_input_tool_name(step_id, steps):
    """Get the string with the name of the tool that generated an input"""
    inp_provenance = ''
    inp_prov_id = str(step_id)
    if inp_prov_id in steps:
        inp_provenance = "(output of **%s** {%% icon tool %%})" % steps[inp_prov_id]['name']
    return inp_provenance


def get_tool_input(tool_desc):
    """Get a dictionary with label being the tool parameter name and the value the description
    of the parameter extracted from the show_tool function of bioblend"""
    tool_inp = collections.OrderedDict()
    for inp in tool_desc["inputs"]:
        tool_inp.setdefault(inp['name'], inp)
    return tool_inp


def format_inputs(wf_inputs, tp_desc, wf_steps):
    inputlist = ''
    for inp_n, inp in wf_inputs.items():
        if inp_n != tp_desc['name']:
            continue
        inps = []
        if isinstance(inp, list): # multiple input (not collection)
            icon = 'param-files'
            for i in inp:
                inps.append('`%s` %s' % (
                    i['output_name'],
                    get_input_tool_name(i['id'], wf_steps)))
        else: # sinle input
            icon = 'param-file'
            inps = ['`%s` %s' % (
                inp['output_name'],
                get_input_tool_name(inp['id'], wf_steps))]
        context = {
            "icon": icon,
            "input_name": tp_desc['label'],
            "input_value": ', '.join(inps)
        }
        inputlist += templates.render(INPUT_FILE_TEMPLATE, **context)
    return inputlist


def format_param_desc(wf_params, wf_inputs, tp_desc, level, wf_steps):
    """Format the parameter description (label and value) given the type of parameter"""
    paramlist = ''
    if 'type' not in tp_desc:
        raise ValueError("No type for the paramater %s" % tp_desc['name'])
    if tp_desc['type'] == 'data':
        paramlist += format_inputs(wf_inputs, tp_desc, wf_steps)
    elif tp_desc['type'] == 'data_collection':
        info("data_collection parameters are currently not supported")
    elif tp_desc['type'] == 'section':
        context = {'space': '   ' * level, 'section_label': tp_desc['title']}
        sub_param_desc = get_param_desc(wf_params, wf_inputs, get_tool_input(tp_desc), level+1, wf_steps)
        if sub_param_desc != '':
            paramlist += templates.render(INPUT_SECTION, **context)
            paramlist += sub_param_desc
    elif tp_desc['type'] == 'conditional':
        test_param = tp_desc['test_param']
        paramlist += format_param_desc(wf_params[test_param['name']], wf_inputs, test_param, level, wf_steps)
        for case in tp_desc['cases']:
            if case['value'] == wf_params[test_param['name']]:
                if len(case['inputs']) > 1:
                    paramlist += get_param_desc(wf_params, wf_inputs, get_tool_input(case), level+1, wf_steps)
    elif tp_desc['type'] == 'repeat':
        repeat_inp_desc = get_tool_input(tp_desc)
        context = {'space': '   ' * level, 'section_label': tp_desc['title']}
        paramlist += templates.render(INPUT_SECTION, **context)
        for r in range(len(wf_params)):
            context = {
                'space': '   ' * (level+1),
                'section_label': "%s: %s" % (r+1, tp_desc['title'])}
            paramlist += templates.render(INPUT_SECTION, **context)
            paramlist += get_param_desc(wf_params[r], wf_inputs, repeat_inp_desc, level+2, wf_steps)
            if r < len(wf_params) - 1:
                context = {'space': '   ' * (level+1), 'section_label': tp_desc['title']}
                paramlist += templates.render(INPUT_ADD_REPEAT, **context)
    else:
        if tp_desc['value'] == wf_params:
            return paramlist
        elif tp_desc['type'] == 'boolean':
            if bool(tp_desc['value']) == wf_params:
                return paramlist
            param_value = 'Yes' if wf_params else 'No'
        elif tp_desc['type'] == 'select':
            param_value = ''
            for opt in tp_desc['options']:
                if opt[1] == wf_params:
                    param_value = opt[0]
        elif tp_desc['type'] == 'data_column':
            param_value = "c%s" % wf_params
        else:
            param_value = wf_params
        context = {
            'space': '   ' * level,
            'param_label': tp_desc['label'],
            'param_value': param_value}
        paramlist += templates.render(INPUT_PARAM, **context)
    return paramlist


def get_param_desc(wf_params, wf_inputs, tp_desc, level, wf_steps):
    """Parse the parameters of the tool and return the formatted list of the
    parameters and values set in the workflow"""
    paramlist = ''
    for n, tp_d in tp_desc.items():
        if n not in wf_params:
            raise ValueError("%s not in workflow" % n)
        wf_param = wf_params[n]
        if isinstance(wf_param, str) and ":" in wf_param:
            wf_param = json.loads(wf_param)
        paramlist += format_param_desc(wf_param, wf_inputs, tp_d, level, wf_steps)
    return paramlist


def get_handson_box(step_id, steps, tools):
    """Get the string for an hands-on box based on a step in a workflow"""
    step = steps[step_id]

    # get tool
    tool_name = step['name']
    if len(step['input_connections']) == 0:
        return ''
    tp_desc = tools[tool_name]

    # add description
    wf_inputs = step['input_connections']
    wf_params = json.loads(step['tool_state'])
    paramlist = get_param_desc(wf_params, wf_inputs, tp_desc, 1, steps)

    context = {"tool_name": tool_name, "paramlist": paramlist}
    return templates.render(HANDS_ON_TOOL_BOX_TEMPLATE, **context)


def get_wf_from_running_galaxy(kwds, ctx):
    """Get the workflow dictionary from a running Galaxy instance with the workflow installed there"""
    gi = galaxy.GalaxyInstance(kwds['galaxy_url'], key=kwds['galaxy_api_key'])
    wf = gi.workflows.export_workflow_dict(kwds['workflow_id'])
    tools = get_wf_tool_description(wf, gi)
    return wf, tools


def get_wf_tool_description(wf, gi):
    """Get a dictionary with description of inputs of all tools in a workflow"""
    tools = {}
    for s in wf['steps']:
        step = wf['steps'][s]
        if len(step['input_connections']) == 0:
            continue
        try:
            tool_desc = gi.tools.show_tool(step['tool_id'], io_details = True)
        except:
            tool_desc = {'inputs': []}
        tools.setdefault(step['name'], get_tool_input(tool_desc))
    return tools


def serve_wf_locally(kwds, wf_filepath, ctx):
    """Server local Galaxy and get the workflow dictionary"""
    assert is_galaxy_engine(**kwds)
    runnable = for_path(wf_filepath)
    with engine_context(ctx, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([runnable]) as config:
            workflow_id = config.workflow_id(wf_filepath)
            wf = config.gi.workflows.export_workflow_dict(workflow_id)
            tools = get_wf_tool_description(wf, config.gi) 
    return wf, tools


def create_tutorial_from_workflow(kwds, z_file_links, tuto_dir, ctx):
    """Create tutorial structure from the workflow file"""
    # load workflow
    if kwds['workflow_id']:
        if not kwds['galaxy_url']:
            raise ValueError("No Galaxy URL given")
        if not kwds['galaxy_api_key']:
            raise ValueError("No API key to access Galaxy given")
        wf, tools = get_wf_from_running_galaxy(kwds, ctx)           
    else:
        wf, tools = serve_wf_locally(kwds, kwds["workflow"], ctx)

    body = ''
    for step in range(len(wf['steps'].keys())):
        body += get_handson_box(str(step), wf['steps'], tools)

    context = {
        "topic_name": kwds["topic_name"],
        "tutorial_name": kwds["tutorial_name"],
        "zenodo_link": kwds["zenodo"] if kwds["zenodo"] else '',
        "z_file_links": "\n>    ".join(z_file_links),
        "body": body
    }
    template = templates.render(TUTORIAL_TEMPLATE, **context)
    
    # create the tutorial markdown file
    md_path = os.path.join(tuto_dir, "tutorial.md")
    with open(md_path, 'w') as md:
        md.write(template)


def create_tutorial(kwds, tuto_dir, topic_dir, template_dir, ctx):
    """Create the skeleton of a new tutorial"""
    # copy or rename templates
    template_tuto_path = os.path.join(topic_dir, "tutorials", "tutorial1")
    if os.path.isdir(template_tuto_path):
        os.rename(template_tuto_path, tuto_dir)
    else:
        shutil.copytree(template_dir, tuto_dir)

    # extract the data library from Zenodo and the links for the tutorial
    z_file_links = ''
    if kwds["zenodo"]:
        info("Create the data library from Zenodo")
        z_file_links = extract_from_zenodo(kwds, tuto_dir)

    # create tutorial skeleton from workflow
    if kwds["workflow"] or kwds['workflow_id']:
        info("Create tutorial skeleton from workflow")
        create_tutorial_from_workflow(kwds, z_file_links, tuto_dir, ctx)

    # fill the metadata of the new tutorial
    update_tutorial(kwds, tuto_dir, topic_dir)


def init(ctx, kwds):
    """Create/update a topic/tutorial"""
    topic_template_dir = get_template_dir()

    topic_dir = os.path.join("topics", kwds['topic_name'])
    if not os.path.isdir(topic_dir):
        info("The topic %s does not exist. It will be created" % kwds['topic_name'])
        create_topic(kwds, topic_dir, topic_template_dir)

    if not kwds['tutorial_name']:
        if kwds['workflow'] or kwds['workflow_id']:
            raise Exception("A tutorial name is needed to create the skeleton of the tutorial from a workflow")
        if kwds['zenodo']:
            raise Exception("A tutorial name is needed to add Zenodo information")
    else:
        tuto_dir = os.path.join(topic_dir, "tutorials", kwds['tutorial_name'])
        if not os.path.isdir(tuto_dir):
            tuto_template_dir = os.path.join(topic_template_dir, "tutorials", "tutorial1")
            info("The tutorial %s in topic %s does not exist. It will be created." % (kwds['tutorial_name'], kwds['topic_name']))
            create_tutorial(kwds, tuto_dir, topic_dir, tuto_template_dir, ctx)
        else:
            info("The tutorial %s in topic %s already exists. It will be updated with the other arguments" % (kwds['tutorial_name'], kwds['topic_name']))
            update_tutorial(kwds, tuto_dir, topic_dir)


def fill_data_library(ctx, kwds):
    """Fill a data library for a tutorial"""
    topic_dir = os.path.join("topics", kwds['topic_name'])
    if not os.path.isdir(topic_dir):
        raise Exception("The topic %s does not exists. It should be created" % kwds['topic_name'])

    tuto_dir = os.path.join(topic_dir, "tutorials", kwds['tutorial_name'])
    if not os.path.isdir(tuto_dir):
        raise Exception("The tutorial %s does not exists. It should be created" % kwds['tutorial_name'])

    # get metadata
    metadata_path = os.path.join(topic_dir, "metadata.yaml")
    metadata = load_yaml(metadata_path)
    tuto_metadata = collections.OrderedDict()

    for mat in metadata['material']:
        if mat['name'] == kwds['tutorial_name']:
            tuto_metadata = mat

    # get the zenodo link
    z_link = ''
    if 'zenodo_link' in tuto_metadata and tuto_metadata['zenodo_link'] != '':
        if kwds['zenodo']:
            info("The data library and the metadata will be updated with the new Zenodo link")
            tuto_metadata['zenodo_link'] = z_link
            z_link = kwds['zenodo']
        else:
            info("The data library will be extracted using the Zenodo link in the metadata")
            z_link = tuto_metadata['zenodo_link']
    else:
        info("The data library will be created and the metadata will be filled with the new Zenodo link")
        tuto_metadata['zenodo_link'] = z_link
        z_link = kwds['zenodo']
    
    if z_link == '' or z_link is None:
        raise Exception("A Zenodo link should be provided either in the metadata file or as argument of the command")

    # extract the data library from Zenodo
    topic_kwds = {
        'topic_title': metadata['title'],
        'topic_summary': metadata['summary'],
        'tutorial_title': tuto_metadata['title'],
        'zenodo': z_link,
        'datatypes': kwds['datatypes']
    }
    extract_from_zenodo(topic_kwds, tuto_dir)

    # update the metadata
    save_to_yaml(metadata, metadata_path)