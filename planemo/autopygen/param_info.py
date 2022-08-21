import dataclasses
from typing import Any, List, Tuple, Union, Optional


@dataclasses.dataclass
class ParamInfo:
    """
    Class containing data of a single extracted parameter
    """
    type: str
    name: str
    argument: str
    label: str
    section: str
    section_label: str
    default_val: Any
    custom_attributes: List[Tuple[str, str]]

    nargs: Union[float, int] = 0
    help: Optional[str] = None
    optional: bool = False
    is_repeat: bool = False
    is_select: bool = False
    choices: Optional[List[Any]] = None
    is_flag: bool = False
    format: Optional[str] = None