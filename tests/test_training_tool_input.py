"""Training:tool_input functions."""
import json
import os

import pytest

from planemo.training.tool_input import (
    get_empty_input,
    get_empty_param,
    get_input_tool_name,
    ToolInput,
)
from .test_training import (
    wf,
    wf_param_values,
)
from .test_utils import TEST_DATA_DIR

wf_steps = wf["steps"]
# load the output from
# gi.tools.show_tool('toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/2.0.0', io_details=True)
with open(os.path.join(TEST_DATA_DIR, "training_query_tabular.json")) as tool_desc_f:
    tool_desc = json.load(tool_desc_f)
tool_inp_desc = tool_desc["inputs"]


def test_get_input_tool_name() -> None:
    """Test :func:`planemo.training.tool_input.get_input_tool_name`."""
    assert "Input dataset" in get_input_tool_name("1", wf_steps)
    assert "output of" in get_input_tool_name("4", wf_steps)
    assert get_input_tool_name("10", wf_steps) == ""


def test_get_empty_input() -> None:
    """Test :func:`planemo.training.tool_input.get_empty_input`."""
    assert '{% icon param-file %} *"Input file"*: File' in get_empty_input()


def test_get_empty_param() -> None:
    """Test :func:`planemo.training.tool_input.get_empty_param`."""
    assert '*"Parameter"*: `a value`' in get_empty_param()


def test_ToolInput_init() -> None:
    """Test :func:`planemo.training.tool_input.ToolInput.init`."""
    # test type exception
    with pytest.raises(Exception, match="No type for the parameter t"):
        ToolInput(
            tool_inp_desc={"name": "t"},
            wf_param_values=wf_param_values,
            wf_steps=wf_steps,
            level=1,
            should_be_there=False,
            force_default=False,
        )
    # test with param not in workflow and exception
    with pytest.raises(Exception, match="t not in workflow"):
        ToolInput(
            tool_inp_desc={"name": "t", "type": ""},
            wf_param_values=wf_param_values,
            wf_steps=wf_steps,
            level=1,
            should_be_there=True,
            force_default=False,
        )
    # test with param not in workflow but no exception
    tool_input = ToolInput(
        tool_inp_desc={"name": "t", "type": ""},
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=False,
        force_default=False,
    )
    assert "save_db" in tool_input.wf_param_values
    # test with param in workflow
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[0],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=False,
        force_default=False,
    )
    assert "save_db" not in tool_input.wf_param_values
    assert tool_input.wf_param_values == "workdb.sqlite"


def test_ToolInput_get_formatted_inputs() -> None:
    """Test :func:`planemo.training.tool_input.ToolInput.get_formatted_inputs`."""
    # test no input
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[1]["inputs"][0],
        wf_param_values={},
        wf_steps=wf_steps,
        level=1,
        should_be_there=False,
        force_default=False,
    )
    inputlist = tool_input.get_formatted_inputs()
    assert inputlist == ""
    # test collection
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[1]["inputs"][0],
        wf_param_values=wf_param_values["add_to_database"],
        wf_steps=wf_steps,
        level=1,
        should_be_there=False,
        force_default=False,
    )
    inputlist = tool_input.get_formatted_inputs()
    assert "param-collection" in inputlist
    assert "(Input dataset collection)" in inputlist
    # test single input
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[2]["inputs"][0],
        wf_param_values=wf_param_values["tables"][0],
        wf_steps=wf_steps,
        level=1,
        should_be_there=False,
        force_default=False,
    )
    inputlist = tool_input.get_formatted_inputs()
    assert "param-file" in inputlist
    assert "(Input dataset)" in inputlist


def test_ToolInput_get_lower_param_desc() -> None:
    """Test :func:`planemo.training.tool_input.ToolInput.get_lower_param_desc`."""
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[1],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    sub_param_desc = tool_input.get_lower_param_desc()
    assert ">        - {% icon param-collection %}" in sub_param_desc


def test_ToolInput_get_formatted_section_desc() -> None:
    """Test :func:`planemo.training.tool_input.ToolInput.get_formatted_section_desc`."""
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[1],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    section_paramlist = tool_input.get_formatted_section_desc()
    assert '>    - In *"' in section_paramlist
    assert ">        - {%" in section_paramlist


def test_ToolInput_get_formatted_conditional_desc() -> None:
    """Test :func:`planemo.training.tool_input.ToolInput.get_formatted_conditional_desc`."""
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[5],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    conditional_paramlist = tool_input.get_formatted_conditional_desc()
    assert '>    - *"' in conditional_paramlist
    assert '"*: `Yes`' in conditional_paramlist
    assert '>        - *"' in conditional_paramlist


def test_ToolInput_get_formatted_repeat_desc() -> None:
    """Test :func:`planemo.training.tool_input.ToolInput.get_formatted_repeat_desc`."""
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[2],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    repeat_desc = tool_input.get_formatted_repeat_desc()
    assert '>    - In *"' in repeat_desc
    assert '>        - {% icon param-repeat %} *"Insert' in repeat_desc
    assert ">            -" in repeat_desc


def test_ToolInput_get_formatted_other_param_desc() -> None:
    """Test :func:`planemo.training.tool_input.ToolInput.get_formatted_other_param_desc`."""
    # test default value of the tool
    tool_input = ToolInput(
        tool_inp_desc={"value": 10, "name": "t", "type": ""},
        wf_param_values={"t": 10},
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    assert tool_input.get_formatted_other_param_desc() == ""
    # test boolean parameter
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[3],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    assert tool_input.get_formatted_other_param_desc() == ""
    tool_input.wf_param_values = "true"
    assert "*: `Yes`" in tool_input.get_formatted_other_param_desc()
    # test select parameter
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[5]["cases"][0]["inputs"][0],
        wf_param_values=wf_param_values["query_result"],
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    assert "*: `&`" in tool_input.get_formatted_other_param_desc()
    # test other parameter
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[4],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=True,
    )
    assert "*: ``" in tool_input.get_formatted_other_param_desc()


def test_ToolInput_get_formatted_desc() -> None:
    """Test :func:`planemo.training.tool_input.ToolInput.get_formatted_desc`."""
    # test no param values
    tool_input = ToolInput(
        tool_inp_desc={"value": 10, "name": "t", "type": ""},
        wf_param_values={},
        wf_steps=wf_steps,
        level=1,
        should_be_there=False,
        force_default=False,
    )
    assert tool_input.get_formatted_desc() == ""
    # test data
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[2]["inputs"][0],
        wf_param_values=wf_param_values["tables"][0],
        wf_steps=wf_steps,
        level=1,
        should_be_there=False,
        force_default=False,
    )
    inputlist = tool_input.get_formatted_inputs()
    formatted_desc = tool_input.get_formatted_desc()
    assert inputlist == formatted_desc
    # test section
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[1],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    section_paramlist = tool_input.get_formatted_section_desc()
    formatted_desc = tool_input.get_formatted_desc()
    assert section_paramlist == formatted_desc
    # test conditional
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[5],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    conditional_paramlist = tool_input.get_formatted_conditional_desc()
    formatted_desc = tool_input.get_formatted_desc()
    assert conditional_paramlist == formatted_desc
    # test repeat
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[2],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    repeat_desc = tool_input.get_formatted_repeat_desc()
    formatted_desc = tool_input.get_formatted_desc()
    assert repeat_desc == formatted_desc
    # test other
    tool_input = ToolInput(
        tool_inp_desc=tool_inp_desc[3],
        wf_param_values=wf_param_values,
        wf_steps=wf_steps,
        level=1,
        should_be_there=True,
        force_default=False,
    )
    param_desc = tool_input.get_formatted_other_param_desc()
    formatted_desc = tool_input.get_formatted_desc()
    assert param_desc == formatted_desc
