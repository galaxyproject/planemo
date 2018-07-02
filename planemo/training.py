""" gtdk: Galaxy training development kit """

import collections
import os
import shutil
import oyaml as yaml
from planemo.io import info


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

    elif filepath.is_file():
        filepath.unlink()


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
            mat["hands_on"] = kwds["tutorial_hands_on"]
            mat["slides"] = kwds["tutorial_slides"]
            found = True

    if not found:
        new_mat = collections.OrderedDict()
        new_mat["title"] = kwds["tutorial_title"]
        new_mat["name"] = kwds["tutorial_name"]
        new_mat["type"] = 'tutorial'
        new_mat["zenodo_link"] = ''
        new_mat["hands_on"] = kwds["tutorial_hands_on"]
        new_mat["slides"] = kwds["tutorial_slides"]
        new_mat["workflows"] = False
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
    update_tuto_file(hand_on_path, kwds["tutorial_hands_on"], kwds["topic_name"], kwds["tutorial_name"])

    # update the metadata in top of slides.md or remove it if not needed
    slides_path = os.path.join(tuto_dir, "slides.html")
    update_tuto_file(slides_path, kwds["tutorial_slides"], kwds["topic_name"], kwds["tutorial_name"])


def create_tutorial(args, tuto_dir, topic_dir, template_dir):
    """Create the skeleton of a new tutorial"""
    # copy or rename templates
    template_tuto_path = os.path.join(topic_dir, "tutorials", "tutorial1")
    if template_tuto_path.isdir():
        template_tuto_path.rename(tuto_dir)
    else:
        shutil.copytree(template_dir, tuto_dir)

    # fill the metadata of the new tutorial
    update_tutorial(args, tuto_dir, topic_dir)


def init(kwds):
    """Create/update a topic/tutorial"""
    topic_template_dir = get_template_dir()

    topic_dir = os.path.join("topics", kwds['topic_name'])
    if not os.path.isdir(topic_dir):
        info("The topic %s does not exist. It will be created" % kwds['topic_name'])
        create_topic(kwds, topic_dir, topic_template_dir)

    if kwds['tutorial_name']:
        tuto_dir = os.path.join("topic_dir", "tutorials", kwds['tutorial_name'])
        if not os.path.is_dir(tuto_dir):
            tuto_template_dir = os.path.join(topic_dir, "tutorials", "tutorial1")
            info("The tutorial %s in topic %s does not exist. It will be created." % (args.tutorial_name, args.topic_name))
            create_tutorial(kwds, tuto_dir, topic_dir, tuto_template_dir)
        else:
            info("The tutorial %s in topic %s already exists. It will be updated with the other arguments" % (args.tutorial_name, args.topic_name))
            update_tutorial(kwds, tuto_dir, topic_dir)
#
#
#
    #if output is None:
    #    output = os.path.splitext(workflow_path)[0] + ".ga"
#
    #runnable = for_path(workflow_path)
    #with engine_context(ctx, **kwds) as galaxy_engine:
    #    with galaxy_engine.ensure_runnables_served([runnable]) as config:
    #        workflow_id = config.workflow_id(workflow_path)
    #        output_dict = config.gi.workflows.export_workflow_dict(workflow_id)
    #        print(output_dict)
    #        import time
    #        time.sleep(10000)
    #        write_file(output, "Test File Here...", force=force)
#