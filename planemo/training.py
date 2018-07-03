""" gtdk: Galaxy training development kit """

import collections
import json
import os
import requests
import shutil
import time
import oyaml as yaml


from planemo import templates
from planemo.io import info
from planemo.runnable import for_path
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)


INPUT_TEMPLATE = """
>   - icon {{icon}} *"{{input_name}}"*: {{input_value}}
"""

INPUT_TEMPLATE_2 = """
>   - {{ '{%' }} icon {{icon}} {{ '%}' }} *"{{input_name}}"*: {{input_value}}
"""


HANDS_ON_TOOL_BOX_TEMPLATE = """
> ### {{ '{%' }} icon hands_on {{ '%}' }} Hands-on: TODO: task description
>
> 1. **{{tool_name}}** {{ '{%' }} icon tool {{ '%}' }} with the following parameters:{{inputlist}}{{paramlist}}
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->
> ### {{ '{%' }} icon question {{ '%}' }} Questions
>
> 1. Question1?
> 2. Question2?
>
>    > ### {{ '{%' }} icon solution {{ '%}' }} Solution
>    >
>    > 1. Answer for question1
>    > 2. Answer for question2
>    >
>    {: .solution}
>
{: .question}
"""


TUTORIAL_TEMPLATE = """
---
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

# Different steps

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


def get_zenodo_file_url(zenodo_link):
    """Get the list of URLs of the files on Zenodo"""
    links = []
    if not zenodo_link:
        return links

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

    # extract the URLs from the JSON
    if 'files' not in req_res:
        return links
    
    for f in req_res['files']:
        links.append(f['links']['self'])

    return links


def get_input_tool_name(step_id, steps):
    """Get the string with the name of the tool that generated an input"""
    inp_provenance = ''
    inp_prov_id = str(step_id)
    if inp_prov_id in steps:
        inp_provenance = '(output of `%s`)' % steps[inp_prov_id]['name']
    return inp_provenance


def get_input_label(inp_n, inputs):
    """Get the label of an input"""
    #for inp in inputs:
    #    if inp["name"] == inp_n:
    #        return inp["label"]
    return inp_n


def get_handson_box(step_id, steps, tools):
    """Get the string for an hands-on box based on a step in a workflow"""
    step = steps[step_id]

    # get tool
    tool_name = step['name']
    if len(step['input_connections']) == 0:
        return ''
    tool = {}#tools[tool_name]

    # add input description
    input_conn = step['input_connections']
    inputlist = ''
    for inp_n, inp in input_conn.items():
        inps = []
        if isinstance(inp, list): # multiple input (not collection)
            icon = 'param-files'
            for i in inp:
                inps.append('`%s` %s' % (
                    i['output_name'],
                    get_input_tool_name(i['id'], steps)))
        else: # sinle input
            icon = 'param-file'
            inps = ['`%s` %s' % (
                inp['output_name'],
                get_input_tool_name(inp['id'], steps))]

        context = {
            "icon": icon,
            "input_name": get_input_label(inp_n, tool["inputs"]),
            "input_value": ', '.join(inps)
        }
        inputlist += templates.render(INPUT_TEMPLATE, **context)

    # add parameters
    parameters = step['tool_state']
    print(parameters)

    #g = nested_dict_iter(json.loads(parameters))
    #print(g)
    
    paramlist = ''

    # while True:
    #    try:
    #        (k, v) = next(g)
    #        print("param: ", k, v)
    #    except StopIteration:
    #        break

    #    if not v or v == 'null' or v == '[]':
    #        pass
    #    elif 'RuntimeValue' in str(v):
    #        pass
            # print("myinputs:", v, inputs)
            # print(inputs)
    #    elif '__' not in k and k != 'chromInfo':
    #        paramlist += '\n>   - *"' + k + '"*: `' + str(v).strip('"[]') + '`'

    # print(paramlist)

    context = {
        "tool_name": tool_name,
        "inputlist": inputlist,
        "paramlist": paramlist
    }
    return templates.render(HANDS_ON_TOOL_BOX_TEMPLATE, **context)


def get_wf_from_running_galaxy(kwds, ctx):
    """Get the workflow dictionary from a running Galaxy instance with the workflow installed there"""
    return {}


def get_wf_tool_description(wf, gi):
    """Get a dictionary with description of all tools in a workflow"""
    tools = {}
    for s in wf['steps']:
        step = wf['steps'][s]
        if len(step['input_connections']) == 0:
            continue
        print()
        print(step)
        tools.setdefault(step['name'],
                         gi.tools.show_tool(step['tool_id'], io_details = True))
    return tools


def serve_wf_locally(kwds, wf_filepath, ctx):
    """Server local Galaxy and get the workflow dictionary"""
    assert is_galaxy_engine(**kwds)
    runnable = for_path(wf_filepath)
    with engine_context(ctx, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([runnable]) as config:
            workflow_id = config.workflow_id(wf_filepath)
            wf = config.gi.workflows.export_workflow_dict(workflow_id)
            print(wf)
            tools = {} # get_wf_tool_description(wf, config.gi) 
    return wf, tools


def create_tutorial_from_workflow(kwds, tuto_dir, ctx):
    """Create tutorial structure from the workflow file"""
    # load workflow
    if kwds['workflow_id']:
        if kwds['galaxy_url']:
            wf = get_wf_from_running_galaxy(kwds, ctx)
    else:
        wf, tools = serve_wf_locally(kwds, kwds["workflow"], ctx)

    # get 
    z_file_links = get_zenodo_file_url(kwds['zenodo'])

    body = ''
    for step in wf['steps']:
        body += get_handson_box(step, wf['steps'], tools)

    context = {
        "topic_name": kwds["topic_name"],
        "tutorial_name": kwds["tutorial_name"],
        "zenodo_link": kwds["zenodo"] if kwds["zenodo"] else '',
        "z_file_links": "\n>    ".join(z_file_links),
        "hands_on_boxes": body
    }
    template = templates.render(TUTORIAL_TEMPLATE, **context)
    
    # create the tutorial markdown file
    md_path = os.path.join(tuto_dir, "tutorial.md")
    with open(md_path, 'w') as md:
        md.write(template)


def extract_tools_from_workflow(kwds, tuto_dir):
    """Create and fill tools.yaml file from workflow"""
    info("Test")


def extract_data_library_from_zenodo(zenodo_link, tuto_dir):
    """Create the data_library from Zenodo"""
    info("Test")


def create_tutorial(kwds, tuto_dir, topic_dir, template_dir, ctx):
    """Create the skeleton of a new tutorial"""
    # copy or rename templates
    template_tuto_path = os.path.join(topic_dir, "tutorials", "tutorial1")
    if os.path.isdir(template_tuto_path):
        os.rename(template_tuto_path, tuto_dir)
    else:
        shutil.copytree(template_dir, tuto_dir)

    print(kwds)
    # create tutorial skeleton from workflow
    if kwds["workflow"] or kwds['workflow_id']:
        info("Create tutorial skeleton from workflow")
        create_tutorial_from_workflow(kwds, tuto_dir, ctx)
        extract_tools_from_workflow(kwds, tuto_dir)

    # fill the metadata of the new tutorial
    update_tutorial(kwds, tuto_dir, topic_dir)

    # extract the data library from Zenodo
    if kwds["zenodo"]:
        extract_data_library_from_zenodo(kwds["zenodo"], tuto_dir)


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
