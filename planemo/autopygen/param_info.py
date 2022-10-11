import dataclasses
from typing import Any, List, Union, Optional, Dict


@dataclasses.dataclass
class ParamTypeFlags:
    # No flags => store
    # ===============================================
    # is_flag => store_const, store_true, store_false
    is_flag: bool = False
    # is_repeat => append, append_const, count
    is_repeat: bool = False
    # is_selection => has choices
    is_selection: bool = False
    # is_extend => extend
    is_extend: bool = False
    # is_version => version action
    is_version: bool = False
    # is_help => help action
    is_help: bool = False
    # is_output => could be parameter for result output file
    is_output: bool = False


@dataclasses.dataclass
class ParamInfo:
    """
    Class containing data of a single extracted parameter
    """
    is_positional: bool
    param_type: ParamTypeFlags
    type: str
    name: str
    argument: str
    label: str
    section: str
    section_label: str
    default_val: Any
    custom_attributes: Dict[str, str]

    nargs: Union[float, int] = 0
    help: Optional[str] = None
    optional: bool = False
    choices: Optional[List[Any]] = None
    format: Optional[str] = None
