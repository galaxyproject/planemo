import json
import os
import shutil

from nose.tools import assert_raises_regexp

from planemo import cli
from planemo import training
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.runnable import for_path
from .test_utils import (
    PROJECT_TEMPLATES_DIR,
    TEST_DATA_DIR
)


METADATA_FP = os.path.join(TEST_DATA_DIR, "training_metadata_w_zenodo.yaml")
TRAINING_TEMPLATE_DIR = os.path.join(PROJECT_TEMPLATES_DIR, "training")
TUTORIAL_FP = os.path.join("tutorials", "tutorial1", "tutorial.md")
DATATYPE_FP = os.path.join(TEST_DATA_DIR, "training_datatypes.yaml")
ZENODO_LINK = 'https://zenodo.org/record/1321885'
WF_FP = os.path.join(TEST_DATA_DIR, "training_workflow.ga")
RUNNABLE = for_path(WF_FP)
CTX = cli.Context()
CTX.planemo_directory = "/tmp/planemo-test-workspace"


def prepare_test():
    topic_name = 'my_new_topic'
    topic_dir = topic_name
    tuto_name = "new_tuto"
    tuto_dir = os.path.join(topic_dir, "tutorials", tuto_name)
    kwds = {
        'topic_name': topic_name,
        'topic_title': "New topic",
        'topic_target': "use",
        'topic_summary': "Topic summary",
        'tutorial_name': tuto_name,
        'tutorial_title': "Title of tuto",
        'hands_on': True,
        'slides': True,
        'workflow': None,
        'workflow_id': None,
        'zenodo': None,
        'datatypes': DATATYPE_FP,
        'templates': None,
        # planemo configuation
        'conda_auto_init': True,
        'conda_auto_install': True,
        'conda_copy_dependencies': False,
        'conda_debug': False,
        'conda_dependency_resolution': False,
        'conda_ensure_channels': 'iuc,bioconda,conda-forge,defaults',
        'conda_exec': None,
        'conda_prefix': None,
        'conda_use_local': False,
        'brew_dependency_resolution': False,
        'daemon': False,
        'database_connection': None,
        'database_type': 'auto',
        'dependency_resolvers_config_file': None,
        'docker': False,
        'docker_cmd': 'docker',
        'docker_extra_volume': None,
        'docker_galaxy_image': 'quay.io/bgruening/galaxy',
        'docker_host': None,
        'docker_sudo': False,
        'docker_sudo_cmd': 'sudo',
        'engine': 'galaxy',
        'extra_tools': (),
        'file_path': None,
        'galaxy_api_key': None,
        'galaxy_branch': None,
        'galaxy_database_seed': None,
        'galaxy_email': 'planemo@galaxyproject.org',
        'galaxy_root': None,
        'galaxy_single_user': True,
        'galaxy_source': None,
        'galaxy_url': None,
        'host': '127.0.0.1',
        'ignore_dependency_problems': False,
        'install_galaxy': False,
        'job_config_file': None,
        'mulled_containers': False,
        'no_cleanup': False,
        'no_cache_galaxy': False,
        'no_dependency_resolution': True,
        'non_strict_cwl': False,
        'pid_file': None,
        'port': '9090',
        'postgres_database_host': None,
        'postgres_database_port': None,
        'postgres_database_user': 'postgres',
        'postgres_psql_path': 'psql',
        'profile': None,
        'shed_dependency_resolution': False,
        'shed_install': True,
        'shed_tool_conf': None,
        'shed_tool_path': None,
        'skip_venv': False,
        'test_data': None,
        'tool_data_table': None,
        'tool_dependency_dir': None
    }
    return (kwds, topic_dir, tuto_dir)


def test_load_yaml():
    """Test :func:`planemo.training.load_yaml`."""
    metadata = training.load_yaml(METADATA_FP)
    # test if name there
    assert metadata["name"] == "test"
    # test if order of material is conserved
    assert metadata["material"][1]["name"] == "test"


def test_save_to_yaml():
    """Test :func:`planemo.training.save_to_yaml`."""
    metadata = training.load_yaml(METADATA_FP)
    new_metadata_fp = "metadata.yaml"
    training.save_to_yaml(metadata, new_metadata_fp)
    assert os.path.exists(new_metadata_fp)
    assert open(new_metadata_fp, 'r').read().find('material') != -1
    os.remove(new_metadata_fp)


def test_get_template_dir_1():
    """Test :func:`planemo.training.get_template_dir`: test exception raising"""
    kwds = {"templates": None}
    exp_exception = "This script needs to be run in the training material repository"
    with assert_raises_regexp(Exception, exp_exception):
        training.get_template_dir(kwds)


def test_get_template_dir_2():
    """Test :func:`planemo.training.get_template_dir`: test default return value"""
    kwds = {"templates": None}
    os.makedirs("templates")
    assert training.get_template_dir(kwds) == "templates"
    shutil.rmtree("templates")


def test_get_template_dir_3():
    """Test :func:`planemo.training.get_template_dir`: test return value"""
    template_path = "temp"
    kwds = {"templates": template_path}
    assert training.get_template_dir(kwds) == template_path


def test_update_top_metadata_file_1():
    """Test :func:`planemo.training.update_top_metadata_file`: test topic change."""
    new_index_fp = "index.md"
    topic_name = 'my_new_topic'
    template_index_fp = os.path.join(TRAINING_TEMPLATE_DIR, "index.md")
    shutil.copyfile(template_index_fp, new_index_fp)
    training.update_top_metadata_file(new_index_fp, topic_name)
    assert open(new_index_fp, 'r').read().find(topic_name) != -1
    os.remove(new_index_fp)


def test_update_top_metadata_file_2():
    """Test :func:`planemo.training.update_top_metadata_file`: test tutorial change."""
    new_tuto_fp = "tutorial.md"
    topic_name = 'my_new_topic'
    tuto_name = 'my_new_tuto'
    template_tuto_fp = os.path.join(TRAINING_TEMPLATE_DIR, TUTORIAL_FP)
    shutil.copyfile(template_tuto_fp, new_tuto_fp)
    training.update_top_metadata_file(new_tuto_fp, topic_name, tuto_name=tuto_name)
    assert open(new_tuto_fp, 'r').read().find(tuto_name) != -1
    os.remove(new_tuto_fp)


def test_update_top_metadata_file_3():
    """Test :func:`planemo.training.update_top_metadata_file`: test tutorial change."""
    new_tuto_fp = "tutorial.md"
    topic_name = 'my_new_topic'
    template_tuto_fp = os.path.join(TRAINING_TEMPLATE_DIR, TUTORIAL_FP)
    shutil.copyfile(template_tuto_fp, new_tuto_fp)
    training.update_top_metadata_file(new_tuto_fp, topic_name, keep=False)
    assert not os.path.exists(new_tuto_fp)


def test_create_topic():
    """Test :func:`planemo.training.create_topic`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    topic_name = kwds['topic_name']
    training.create_topic(kwds, topic_dir, TRAINING_TEMPLATE_DIR)
    # check if files has been moved and updated with topic name
    index_fp = os.path.join(topic_dir, "index.md")
    assert os.path.exists(index_fp)
    assert open(index_fp, 'r').read().find(topic_name) != -1
    tuto_fp = os.path.join(topic_dir, TUTORIAL_FP)
    assert os.path.exists(tuto_fp)
    assert open(tuto_fp, 'r').read().find(topic_name) != -1
    # check metadata content
    metadata = training.load_yaml(os.path.join(topic_dir, "metadata.yaml"))
    assert metadata['name'] == topic_name
    # check in metadata directory
    assert os.path.exists(os.path.join("metadata", "%s.yaml" % topic_name))
    # clean
    shutil.rmtree(topic_dir)
    shutil.rmtree("metadata")


def test_update_tutorial():
    """Test :func:`planemo.training.update_tutorial`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    tuto_title = kwds['tutorial_title']
    metadata_fp = os.path.join(topic_dir, "metadata.yaml")
    tuto_fp = os.path.join(tuto_dir, "tutorial.md")
    slides_fp = os.path.join(tuto_dir, "slides.html")
    # create a topic and prepare the tutorial
    training.create_topic(kwds, topic_dir, TRAINING_TEMPLATE_DIR)
    template_tuto_path = os.path.join(topic_dir, "tutorials", "tutorial1")
    os.rename(template_tuto_path, tuto_dir)
    assert open(metadata_fp, 'r').read().find("tutorial1") != -1
    # test update a new tutorial
    training.update_tutorial(kwds, tuto_dir, topic_dir)
    assert open(metadata_fp, 'r').read().find("tutorial1") == -1
    assert open(metadata_fp, 'r').read().find(tuto_title) != -1
    assert os.path.exists(tuto_fp)
    assert os.path.exists(slides_fp)
    # test update an existing tutorial
    new_tuto_title = "A totally new title"
    kwds['tutorial_title'] = new_tuto_title
    kwds['slides'] = False
    training.update_tutorial(kwds, tuto_dir, topic_dir)
    assert open(metadata_fp, 'r').read().find(tuto_title) == -1
    assert open(metadata_fp, 'r').read().find(new_tuto_title) != -1
    assert not os.path.exists(slides_fp)
    # clean
    shutil.rmtree(topic_dir)
    shutil.rmtree("metadata")


def test_get_zenodo_record():
    """Test :func:`planemo.training.get_zenodo_record`."""
    z_record, req_res = training.get_zenodo_record(ZENODO_LINK)
    file_link_prefix = "https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f"
    assert z_record == "1321885"
    assert 'files' in req_res
    assert req_res['files'][0]['type'] in ['rdata', 'csv']
    assert req_res['files'][0]['links']['self'].find(file_link_prefix) != -1
    # check with wrong zenodo link
    z_record, req_res = training.get_zenodo_record('https://zenodo.org/api/records/zenodooo')
    assert z_record is None
    assert 'files' in req_res
    assert len(req_res['files']) == 0


def test_get_zenodo_record_with_doi():
    """Test :func:`planemo.training.get_zenodo_record`: link with DOI."""
    z_link = 'https://doi.org/10.5281/zenodo.1321885'
    z_record, req_res = training.get_zenodo_record(z_link)
    file_link_prefix = "https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f"
    assert z_record == "1321885"
    assert 'files' in req_res
    assert req_res['files'][0]['type'] in ['rdata', 'csv']
    assert req_res['files'][0]['links']['self'].find(file_link_prefix) != -1


def test_get_galaxy_datatype():
    """Test :func:`planemo.training.get_galaxy_datatype`."""
    assert training.get_galaxy_datatype("csv", DATATYPE_FP) == "csv"
    assert training.get_galaxy_datatype("test", DATATYPE_FP) == "strange_datatype"
    assert training.get_galaxy_datatype("unknown", DATATYPE_FP).find("# Please add") != -1


def test_get_files_from_zenodo():
    """Test :func:`planemo.training.get_files_from_zenodo`."""
    files, links, z_record = training.get_files_from_zenodo(ZENODO_LINK, DATATYPE_FP)
    assert z_record == "1321885"
    # test links
    file_link_prefix = "https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f"
    assert links[0].find(file_link_prefix) != -1
    # test files dict
    assert files[0]['url'].find(file_link_prefix) != -1
    assert files[0]['src'] == 'url'
    assert files[0]['info'] == ZENODO_LINK
    assert files[0]['ext'].find("# Please add") != -1
    assert files[1]['ext'] == 'csv'


def test_init_data_lib():
    """Test :func:`planemo.training.init_data_lib`."""
    data_lib_filepath = 'data-library.yaml'
    datalib = training.init_data_lib(data_lib_filepath)
    assert datalib['destination']['name'] == 'GTN - Material'


def test_prepare_data_library():
    """Test :func:`planemo.training.prepare_data_library`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    os.makedirs(tuto_dir)
    files, links, z_record = training.get_files_from_zenodo(ZENODO_LINK, DATATYPE_FP)
    datalib_fp = os.path.join(tuto_dir, "data-library.yaml")
    # test default prepare_data_library
    training.prepare_data_library(files, kwds, z_record, tuto_dir)
    assert os.path.exists(datalib_fp)
    datalib = training.load_yaml(datalib_fp)
    assert datalib['items'][0]['name'] == kwds['topic_title']
    assert datalib['items'][0]['items'][0]['name'] == kwds['tutorial_title']
    assert datalib['items'][0]['items'][0]['items'][0]['name'] == "DOI: 10.5281/zenodo.%s" % z_record
    assert datalib['items'][0]['items'][0]['items'][0]['description'] == "latest"
    assert datalib['items'][0]['items'][0]['items'][0]['items'][0]['url'] == files[0]['url']
    # test adding a new collection for same tutorial
    new_z_record = '124'
    training.prepare_data_library(files, kwds, new_z_record, tuto_dir)
    datalib = training.load_yaml(datalib_fp)
    assert datalib['items'][0]['items'][0]['items'][0]['name'] == "DOI: 10.5281/zenodo.%s" % new_z_record
    assert datalib['items'][0]['items'][0]['items'][0]['description'] == "latest"
    assert datalib['items'][0]['items'][0]['items'][1]['name'] == "DOI: 10.5281/zenodo.%s" % z_record
    assert datalib['items'][0]['items'][0]['items'][1]['description'] == ""
    # test adding a new tutorial
    new_tuto_title = "New title"
    kwds['tutorial_title'] = new_tuto_title
    training.prepare_data_library(files, kwds, z_record, tuto_dir)
    datalib = training.load_yaml(datalib_fp)
    assert datalib['items'][0]['items'][1]['name'] == new_tuto_title
    assert datalib['items'][0]['items'][1]['items'][0]['name'] == "DOI: 10.5281/zenodo.%s" % z_record
    # test adding a new topic
    new_topic_title = "New title"
    kwds['topic_title'] = new_topic_title
    training.prepare_data_library(files, kwds, z_record, tuto_dir)
    datalib = training.load_yaml(datalib_fp)
    assert datalib['items'][1]['name'] == new_topic_title
    assert datalib['items'][1]['items'][0]['name'] == new_tuto_title
    assert datalib['items'][1]['items'][0]['items'][0]['name'] == "DOI: 10.5281/zenodo.%s" % z_record
    # clean
    shutil.rmtree(topic_dir)


def test_prepare_data_library_from_zenodo():
    """Test :func:`planemo.training.prepare_data_library_from_zenodo`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    os.makedirs(tuto_dir)
    datalib_fp = os.path.join(tuto_dir, "data-library.yaml")
    # test prepare_data_library_from_zenodo with no zenodo
    links = training.prepare_data_library_from_zenodo(kwds, tuto_dir)
    assert len(links) == 0
    assert not os.path.exists(datalib_fp)
    # test prepare_data_library_from_zenodo with a zenodo link
    kwds['zenodo'] = ZENODO_LINK
    links = training.prepare_data_library_from_zenodo(kwds, tuto_dir)
    file_link_prefix = "https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f"
    assert links[0].find(file_link_prefix) != -1
    assert os.path.exists(datalib_fp)
    # clean
    shutil.rmtree(topic_dir)


def test_get_tool_input():
    """Test :func:`planemo.training.get_tool_input`."""
    tool_desc = {
        'inputs': [
            {'name': "name1", 'content': 'c'},
            {'name': "name2", 'content': 'c'}
        ]
    }
    tool_inp = training.get_tool_input(tool_desc)
    assert "name1" in tool_inp
    assert 'content' in tool_inp["name1"]
    assert tool_inp["name1"]['content'] == 'c'


def check_tools(tools):
    """Test the tool return from get_wf_tool_description"""
    assert 'FastQC' in tools
    assert 'input_file' in tools['FastQC']


def test_get_wf_tool_description():
    """Test :func:`planemo.training.get_wf_tool_description`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    assert is_galaxy_engine(**kwds)
    with engine_context(CTX, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([RUNNABLE]) as config:
            workflow_id = config.workflow_id(WF_FP)
            wf = config.gi.workflows.export_workflow_dict(workflow_id)
            wf['steps']['10'] = {
                'input_connections': [],
                'tool_id': 'no_input',
                'name': 'with_no_input'
            }
            wf['steps']['11'] = {
                'input_connections': [1],
                'tool_id': 'no_tool',
                'name': 'with_no_tool'
            }
            tools = training.get_wf_tool_description(wf, config.gi)
    check_tools(tools)
    assert 'with_no_input' not in tools
    assert 'with_no_tool' in tools


def check_workflow(wf):
    """Test the worflow return"""
    assert 'steps' in wf
    assert '1' in wf['steps']
    assert 'name' in wf['steps']['1']


def test_get_wf_tool_from_local_galaxy():
    """Test :func:`planemo.training.get_wf_tool_from_local_galaxy`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    wf, tools = training.get_wf_tool_from_local_galaxy(kwds, WF_FP, CTX)
    check_tools(tools)
    check_workflow(wf)


def test_get_wf_tools_from_running_galaxy():
    """Test :func:`planemo.training.get_wf_tools_from_running_galaxy`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    assert is_galaxy_engine(**kwds)
    kwds['galaxy_url'] = 'http://%s:%s' % (kwds['host'], kwds['port'])
    with engine_context(CTX, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([RUNNABLE]) as config:
            workflow_id = config.workflow_id(WF_FP)
            kwds['workflow_id'] = workflow_id
            kwds['galaxy_api_key'] = config.user_api_key
            wf = config.gi.workflows.export_workflow_dict(workflow_id)
            wf, tools = training.get_wf_tools_from_running_galaxy(kwds)
    check_tools(tools)
    check_workflow(wf)


def test_get_input_tool_name():
    """Test :func:`planemo.training.get_input_tool_name`."""
    steps = {'1': {'name': 'Input dataset'}}
    # test if step not found
    tool_name = training.get_input_tool_name(2, steps)
    assert tool_name == ''
    # test if tool is input
    assert training.get_input_tool_name(1, steps) == '(Input dataset)'
    # test if tool is input
    steps['1']['name'] = 'Input dataset collection'
    assert training.get_input_tool_name(1, steps) == '(Input dataset collection)'
    # test if other case
    steps['1']['name'] = 'Tool name'
    assert training.get_input_tool_name(1, steps) == '(output of **Tool name** {% icon tool %})'


def get_wf_a_tools():
    """Get workflow and tool of a workflow"""
    kwds, topic_dir, tuto_dir = prepare_test()
    assert is_galaxy_engine(**kwds)
    with engine_context(CTX, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([RUNNABLE]) as config:
            workflow_id = config.workflow_id(WF_FP)
            wf = config.gi.workflows.export_workflow_dict(workflow_id)
            tools = training.get_wf_tool_description(wf, config.gi)
    return (wf, tools)


def test_format_inputs():
    """Test :func:`planemo.training.format_inputs`."""
    wf, tools = get_wf_a_tools()
    step = wf['steps']['3']
    step_inputs = step['input_connections']
    tool = tools[step['name']]
    inputlist = training.format_inputs(step_inputs, tool['input_file'], wf['steps'], 1)
    assert inputlist.find('param-collection ') != -1
    assert inputlist.find('Input dataset collection') != -1
    inputlist = training.format_inputs(step_inputs, tool['contaminants'], wf['steps'], 1)
    assert inputlist.find('param-file ') != -1


def test_get_wf_step_inputs():
    """Test :func:`planemo.training.get_wf_step_inputs`."""
    step_inp = {
        'tables_1|table': {'output_name': 'output', 'id': 2},
        'add_to_database|withdb': {'output_name': 'output', 'id': 0},
        'tables_0|table': {'output_name': 'output', 'id': 1},
        'add_to_database|tab_0|tt': {'output_name': 'output', 'id': 0},
        'tables_2|section|sect': {'output_name': 'output', 'id': 1},
        'tables_3|tables_0|sect': {'output_name': 'output', 'id': 1}
    }
    step_inputs = training.get_wf_step_inputs(step_inp)
    assert 'tables' in step_inputs
    assert '0' in step_inputs['tables']
    assert 'table' in step_inputs['tables']['0']
    assert '2' in step_inputs['tables']
    assert 'section' in step_inputs['tables']['2']
    assert 'sect' in step_inputs['tables']['2']['section']
    assert 'output_name' in step_inputs['tables']['2']['section']['sect']
    assert 'add_to_database' in step_inputs
    assert 'withdb' in step_inputs['add_to_database']
    assert 'tab' in step_inputs['add_to_database']
    assert '0' in step_inputs['add_to_database']['tab']
    assert 'tt' in step_inputs['add_to_database']['tab']['0']


def test_json_load():
    """Test :func:`planemo.training.json_load`."""
    assert isinstance(training.json_load('{"name": "1"}'), dict)
    assert isinstance(training.json_load("name"), str)


def test_get_lower_params():
    """Test :func:`planemo.training.get_lower_params`."""
    step_params = {'name': '1'}
    assert 'name' in training.get_lower_params(step_params, 'n1')
    assert training.get_lower_params(step_params, 'name') == '1'
    print(training.get_lower_params('{"name": "1"}', 'n1'))
    assert 'name' in training.get_lower_params('{"name": "1"}', 'n1')


def test_get_lower_inputs():
    """Test :func:`planemo.training.get_lower_inputs`."""
    step_inputs = {'name': '1'}
    assert 'name' in training.get_lower_inputs(step_inputs, 'n1')
    assert training.get_lower_inputs(step_inputs, 'name') == '1'


def test_format_section_param_desc():
    """Test :func:`planemo.training.format_section_param_desc`."""
    wf, tools = get_wf_a_tools()
    step = wf['steps']['4']
    step_inputs = training.get_wf_step_inputs(step['input_connections'])
    step_params = json.loads(step['tool_state'])
    tp_desc = tools[step['name']]['add_to_database']
    section_paramlist = training.format_section_param_desc(
        step_params,
        step_inputs,
        tp_desc,
        0,
        wf['steps'])
    assert section_paramlist.find('In *"Add tables to an existing database"*') != -1
    assert section_paramlist.find('icon param-collection') != -1
    assert section_paramlist.find('Input dataset collection') != -1


def test_format_conditional_param_desc():
    """Test :func:`planemo.training.format_conditional_param_desc`."""
    wf, tools = get_wf_a_tools()
    step = wf['steps']['4']
    step_inputs = training.get_wf_step_inputs(step['input_connections'])
    step_params = json.loads(step['tool_state'])
    tp_desc = tools[step['name']]['query_result']
    section_paramlist = training.format_conditional_param_desc(
        step_params,
        step_inputs,
        tp_desc,
        0,
        wf['steps'])
    print(section_paramlist)
    assert section_paramlist.find('column headers') != -1
    assert section_paramlist.find('`Yes`') != -1
    assert section_paramlist.find('column_header line') != -1


def test_format_repeat_param_desc():
    """Test :func:`planemo.training.format_repeat_param_desc`."""
    wf, tools = get_wf_a_tools()
    step = wf['steps']['4']
    step_inputs = training.get_wf_step_inputs(step['input_connections'])
    step_params = json.loads(step['tool_state'])
    tp_desc = tools[step['name']]['tables']
    repeat_paramlist = training.format_repeat_param_desc(
        step_params,
        step_inputs,
        tp_desc,
        0,
        wf['steps'])
    print(repeat_paramlist)
    assert repeat_paramlist.find('Click on *"Insert Database Table"*') != -1
    assert repeat_paramlist.find('In *"1: Database Table"*') != -1
    assert repeat_paramlist.find('In *"1: Database Table"*') != -1
    assert repeat_paramlist.find('Click on *"Insert Filter Tabular Input Lines"*') != -1
    assert repeat_paramlist.find('In *"1: Filter Tabular Input Lines"*:') != -1
    assert repeat_paramlist.find('In *"2: Database Table"*:') != -1


def test_get_param_value():
    """Test :func:`planemo.training.get_param_value`."""
    # test same value
    tp_desc = {'type': 'boolean', 'value': 'same'}
    assert training.get_param_value('same', tp_desc) is None
    # test boolean
    tp_desc = {'type': 'boolean', 'value': 'True'}
    assert training.get_param_value(True, tp_desc) is None
    assert training.get_param_value(False, tp_desc) == 'No'
    # test select
    tp_desc = {'type': 'select', 'options': [['opt1', 'val1'], ['opt2', 'val2']], 'value': ''}
    assert training.get_param_value('val1', tp_desc) == 'opt1'
    # test data_column
    tp_desc = {'type': 'data_column', 'value': ''}
    assert training.get_param_value('1', tp_desc) == 'c1'
    # test integer
    tp_desc = {'type': 'integer', 'value': ''}
    assert training.get_param_value('1', tp_desc) == '1'


def test_format_param_desc():
    """Test :func:`planemo.training.format_param_desc`."""
    wf, tools = get_wf_a_tools()
    step = wf['steps']['4']
    step_inputs = training.get_wf_step_inputs(step['input_connections'])
    step_params = json.loads(step['tool_state'])
    # test section (add_to_database)
    n = 'add_to_database'
    tp_desc = tools[step['name']][n]
    step_param = training.get_lower_params(step_params, n)
    paramlist = training.format_param_desc(
        step_param,
        step_inputs,
        tp_desc,
        0,
        wf['steps'],
        force_default=False)
    assert paramlist.find('In *"Add tables to an existing database"*') != -1
    # test repeat (tables)
    n = 'tables'
    tp_desc = tools[step['name']][n]
    step_param = training.get_lower_params(step_params, n)
    paramlist = training.format_param_desc(
        step_param,
        step_inputs,
        tp_desc,
        0,
        wf['steps'],
        force_default=False)
    assert paramlist.find('In *"1: Filter Tabular Input Lines"*:') != -1
    # test boolean (save_db)
    n = 'save_db'
    tp_desc = tools[step['name']][n]
    step_param = training.get_lower_params(step_params, n)
    paramlist = training.format_param_desc(
        step_param,
        step_inputs,
        tp_desc,
        0,
        wf['steps'],
        force_default=False)
    assert '`Yes`' in paramlist
    # test conditional (query_result)
    n = 'query_result'
    tp_desc = tools[step['name']][n]
    step_param = training.get_lower_params(step_params, n)
    paramlist = training.format_param_desc(
        step_param,
        step_inputs,
        tp_desc,
        0,
        wf['steps'],
        force_default=False)
    assert paramlist.find('Prefix character') != -1
    # no type
    exp_exception = "No type for the paramater name"
    with assert_raises_regexp(Exception, exp_exception):
        training.format_param_desc(
            step_params,
            step_inputs,
            {'name': 'name'},
            0,
            wf['steps'],
            force_default=False)


def test_get_param_desc():
    """Test :func:`planemo.training.get_param_desc`."""
    wf, tools = get_wf_a_tools()
    step_3 = wf['steps']['3']
    step_inputs = training.get_wf_step_inputs(step_3['input_connections'])
    step_params = json.loads(step_3['tool_state'])
    # not in workflow and should be there
    step_4 = wf['steps']['4']
    tp_desc = tools[step_4['name']]
    exp_exception = "workdb not in workflow"
    with assert_raises_regexp(Exception, exp_exception):
        training.get_param_desc(
            step_params,
            step_inputs,
            tp_desc,
            0,
            wf['steps'],
            should_be_there=True)
    # not in workflow
    step_4 = wf['steps']['4']
    tp_desc = tools[step_4['name']]
    paramlist = training.get_param_desc(
            step_params,
            step_inputs,
            tp_desc,
            0,
            wf['steps'])
    assert paramlist == ''
    # correct one
    tp_desc = tools[step_3['name']]
    paramlist = training.get_param_desc(
        step_params,
        step_inputs,
        tp_desc,
        0,
        wf['steps'])
    assert 'param-collection' in paramlist
    assert 'param-file' in paramlist


def test_get_handson_box():
    """Test :func:`planemo.training.get_handson_box`."""
    wf, tools = get_wf_a_tools()
    # test normal step
    hand_boxes = training.get_handson_box(wf['steps']['3'], wf['steps'], tools)
    assert '### {% icon hands_on %}' in hand_boxes
    assert '{% icon tool %} with the following parameters:' in hand_boxes
    assert ': .hands_on' in hand_boxes
    # test input step
    hand_boxes = training.get_handson_box(wf['steps']['1'], wf['steps'], tools)
    assert hand_boxes == ''


def test_create_tutorial_from_workflow():
    """Test :func:`planemo.training.create_tutorial_from_workflow`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    os.makedirs(tuto_dir)
    assert is_galaxy_engine(**kwds)
    kwds['galaxy_url'] = 'http://%s:%s' % (kwds['host'], kwds['port'])
    with engine_context(CTX, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([RUNNABLE]) as config:
            workflow_id = config.workflow_id(WF_FP)
            kwds['workflow_id'] = workflow_id
            kwds['galaxy_api_key'] = config.user_api_key
            training.create_tutorial_from_workflow(kwds, '', tuto_dir, CTX)
    # tests
    tuto_path = os.path.join(tuto_dir, "tutorial.md")
    assert os.path.exists(tuto_path)
    with open(tuto_path, 'r') as tuto:
        tuto_content = tuto.read()
        assert 'topic_name: my_new_topic' in tuto_content
        assert 'tutorial_name: new_tuto' in tuto_content
        assert '> ### Agenda' in tuto_content
        assert '## Get data' in tuto_content
        assert '{% icon tool %} with the following parameters:' in tuto_content
        assert 'no_toc' in tuto_content
        assert '# Conclusion' in tuto_content
    # clean after
    shutil.rmtree(topic_dir)


def test_add_workflow_file():
    """Test :func:`planemo.training.add_workflow_file`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    wf_dir = os.path.join(tuto_dir, "workflows")
    os.makedirs(wf_dir)
    wf_path = os.path.join(wf_dir, "init_workflow.ga")
    # test with workflow on a running instance
    assert is_galaxy_engine(**kwds)
    kwds['galaxy_url'] = 'http://%s:%s' % (kwds['host'], kwds['port'])
    with engine_context(CTX, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([RUNNABLE]) as config:
            workflow_id = config.workflow_id(WF_FP)
            kwds['workflow_id'] = workflow_id
            kwds['galaxy_api_key'] = config.user_api_key
            training.add_workflow_file(kwds, tuto_dir)
    assert os.path.exists(wf_path)
    os.remove(wf_path)
    # test with local workflow
    kwds["workflow"] = WF_FP
    training.add_workflow_file(kwds, tuto_dir)
    assert os.path.exists(wf_path)
    # clean after
    shutil.rmtree(topic_dir)


def test_create_tutorial():
    """Test :func:`planemo.training.create_tutorial`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    kwds["templates"] = TRAINING_TEMPLATE_DIR
    topic_template_dir = training.get_template_dir(kwds)
    metadata_fp = os.path.join(topic_dir, 'metadata.yaml')
    tuto_fp = os.path.join(tuto_dir, 'tutorial.md')
    data_library_fp = os.path.join(tuto_dir, 'data-library.yaml')
    # wo zenodo and wo workflow
    training.create_topic(kwds, topic_dir, topic_template_dir)
    tuto_template_dir = os.path.join(topic_template_dir, "tutorials", "tutorial1")
    training.create_tutorial(kwds, tuto_dir, topic_dir, tuto_template_dir, CTX)
    with open(metadata_fp, 'r') as metadata:
        metadata_content = metadata.read()
        assert 'name: new_tuto' in metadata_content
        assert "zenodo_link: ''" in metadata_content
        assert 'workflows: false' in metadata_content
    assert '**My Tool** {% icon tool %}' in open(tuto_fp, 'r').read()
    assert 'name: "Small test files"' in open(data_library_fp, 'r').read()
    shutil.rmtree(topic_dir)
    shutil.rmtree("metadata")
    # w zenodo and wo workflow
    kwds["zenodo"] = ZENODO_LINK
    training.create_topic(kwds, topic_dir, topic_template_dir)
    tuto_template_dir = os.path.join(topic_template_dir, "tutorials", "tutorial1")
    training.create_tutorial(kwds, tuto_dir, topic_dir, tuto_template_dir, CTX)
    with open(metadata_fp, 'r') as metadata:
        metadata_content = metadata.read()
        assert 'name: new_tuto' in metadata_content
        assert 'zenodo_link: %s' % ZENODO_LINK in metadata_content
        assert 'workflows: false' in metadata_content
    assert '**My Tool** {% icon tool %}' in open(tuto_fp, 'r').read()
    assert 'DOI: 10.5281/zenodo.1321885' in open(data_library_fp, 'r').read()
    shutil.rmtree(topic_dir)
    shutil.rmtree("metadata")
    # w zenodo and w workflow
    kwds["workflow"] = WF_FP
    training.create_topic(kwds, topic_dir, topic_template_dir)
    tuto_template_dir = os.path.join(topic_template_dir, "tutorials", "tutorial1")
    training.create_tutorial(kwds, tuto_dir, topic_dir, tuto_template_dir, CTX)
    with open(metadata_fp, 'r') as metadata:
        metadata_content = metadata.read()
        assert 'name: new_tuto' in metadata_content
        assert 'zenodo_link: %s' % ZENODO_LINK in metadata_content
        assert 'workflows: true' in metadata_content
    assert '**FastQC** {% icon tool %} with the following parameters:' in open(tuto_fp, 'r').read()
    assert 'DOI: 10.5281/zenodo.1321885' in open(data_library_fp, 'r').read()
    assert os.path.exists(os.path.join(tuto_dir, 'workflows', 'init_workflow.ga'))
    shutil.rmtree(topic_dir)
    shutil.rmtree("metadata")


def test_init():
    """Test :func:`planemo.training.init`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    kwds["templates"] = TRAINING_TEMPLATE_DIR
    topic_dir = os.path.join('topics', topic_dir)
    tuto_dir = os.path.join('topics', tuto_dir)
    metadata_fp = os.path.join(topic_dir, 'metadata.yaml')
    tuto_name = kwds['tutorial_name']
    # new topic, no tutorial name but workflow
    kwds['tutorial_name'] = None
    kwds['workflow'] = WF_FP
    exp_exception = "A tutorial name is needed to create the skeleton of the tutorial from a workflow"
    with assert_raises_regexp(Exception, exp_exception):
        training.init(CTX, kwds)
    # no new topic, no tutorial name but zenodo
    kwds['workflow'] = None
    kwds['zenodo'] = ZENODO_LINK
    exp_exception = "A tutorial name is needed to add Zenodo information"
    with assert_raises_regexp(Exception, exp_exception):
        training.init(CTX, kwds)
    # no new topic, new tutorial
    kwds['tutorial_name'] = tuto_name
    kwds['workflow'] = None
    kwds['zenodo'] = None
    training.init(CTX, kwds)
    assert kwds['tutorial_title'] in open(metadata_fp, 'r').read()
    # no new topic, update tutorial
    kwds['tutorial_title'] = 'Totally new tutorial title'
    training.init(CTX, kwds)
    assert 'Totally new tutorial title' in open(metadata_fp, 'r').read()
    # clean after
    shutil.rmtree(topic_dir)
    shutil.rmtree("metadata")


def test_prepare_tuto_update():
    """Test :func:`planemo.training.prepare_tuto_update`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    new_topic_name = 'a_topic'
    topic_dir = os.path.join('topics', new_topic_name)
    # non existing topic
    kwds['topic_name'] = new_topic_name
    exp_exception = "The topic %s does not exists. It should be created" % new_topic_name
    with assert_raises_regexp(Exception, exp_exception):
        training.prepare_tuto_update(kwds)
    # non existing tutorial
    kwds["templates"] = TRAINING_TEMPLATE_DIR
    topic_template_dir = training.get_template_dir(kwds)
    training.create_topic(kwds, topic_dir, topic_template_dir)
    exp_exception = "The tutorial new_tuto does not exists. It should be created"
    with assert_raises_regexp(Exception, exp_exception):
        training.prepare_tuto_update(kwds)


def test_fill_data_library():
    """Test :func:`planemo.training.fill_data_library`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    kwds["templates"] = TRAINING_TEMPLATE_DIR
    topic_dir = os.path.join('topics', topic_dir)
    tuto_dir = os.path.join(topic_dir, 'tutorials', 'new_tuto')
    training.init(CTX, kwds)
    metadata_fp = os.path.join(topic_dir, 'metadata.yaml')
    data_library_fp = os.path.join(tuto_dir, 'data-library.yaml')
    # no Zenodo link
    kwds['zenodo'] = None
    kwds['workflow'] = None
    exp_exception = "A Zenodo link should be provided either in the metadata file or as argument of the command"
    with assert_raises_regexp(Exception, exp_exception):
        training.fill_data_library(CTX, kwds)
    # with a given Zenodo link and no Zenodo in metadata
    kwds['zenodo'] = ZENODO_LINK
    training.fill_data_library(CTX, kwds)
    assert 'DOI: 10.5281/zenodo.1321885' in open(data_library_fp, 'r').read()
    assert 'zenodo_link: %s' % ZENODO_LINK in open(metadata_fp, 'r').read()
    # with a given Zenodo link and Zenodo in metadata
    new_z_link = 'https://zenodo.org/record/1324204'
    kwds['zenodo'] = new_z_link
    training.fill_data_library(CTX, kwds)
    assert 'DOI: 10.5281/zenodo.1324204' in open(data_library_fp, 'r').read()
    assert 'zenodo_link: %s' % new_z_link in open(metadata_fp, 'r').read()
    # with no given Zenodo link
    kwds['zenodo'] = None
    training.fill_data_library(CTX, kwds)
    assert 'DOI: 10.5281/zenodo.1324204' in open(data_library_fp, 'r').read()
    assert 'zenodo_link: %s' % new_z_link in open(metadata_fp, 'r').read()
    # clean after
    shutil.rmtree(topic_dir)
    shutil.rmtree("metadata")


def test_generate_tuto_from_wf():
    """Test :func:`planemo.training.generate_tuto_from_wf`."""
    kwds, topic_dir, tuto_dir = prepare_test()
    topic_dir = os.path.join('topics', topic_dir)
    tuto_dir = os.path.join(topic_dir, 'tutorials', 'new_tuto')
    kwds["templates"] = TRAINING_TEMPLATE_DIR
    training.init(CTX, kwds)
    # no workflow
    kwds['workflow'] = None
    exp_exception = "A path to a local workflow or the id of a workflow on a running Galaxy instance should be provided"
    with assert_raises_regexp(Exception, exp_exception):
        training.generate_tuto_from_wf(CTX, kwds)
    # with workflow
    kwds['workflow'] = WF_FP
    training.generate_tuto_from_wf(CTX, kwds)
    tuto_fp = os.path.join(tuto_dir, 'tutorial.md')
    metadata_fp = os.path.join(topic_dir, 'metadata.yaml')
    assert 'workflows: true' in open(metadata_fp, 'r').read()
    assert '**FastQC** {% icon tool %} with the following parameters:' in open(tuto_fp, 'r').read()
    assert os.path.exists(os.path.join(tuto_dir, 'workflows', 'init_workflow.ga'))
