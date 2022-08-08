"""
Functions used to extract initialized ArgumentParser from target source code,
and to transform it into easily usable ParamInfo class containing
 all the necessary info for 'inputs' element initialization
"""
import ast
import logging
import math
import uuid



from typing import Optional, Set, Any, List, Dict, Tuple, Union
from argparse import ArgumentParser
import dataclasses

from planemo.autopygen.source_file_parsing.constants import LINTER_MAGIC
from planemo.autopygen.source_file_parsing.decoy_parser import DecoyParser
from planemo.autopygen.source_file_parsing.local_module_parsing import \
    handle_local_module_names
from planemo.autopygen.source_file_parsing.parser_discovery_and_init import \
    get_parser_init_and_actions
from planemo.autopygen.source_file_parsing.parsing_commons import \
    create_module_tree_from_path, create_module_tree_from_str
from planemo.autopygen.source_file_parsing.parsing_exceptions import \
    ArgumentParsingDiscoveryError
from planemo.autopygen.source_file_parsing.unknown_names_discovery import \
    initialize_variables_in_module


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


def obtain_and_convert_parser(path: str) -> Optional[ArgumentParser]:
    """
    Function parses python source code located at 'path',
    extracts initialized argument parser from the source file and returns it

    Parameters
    ----------
    path : str
     path to the source file containing argument parser init

    Returns
    -------
    Initialized argument parser, or None, in case error happened
    """
    try:
        tree = create_module_tree_from_path(path)
    except FileNotFoundError:
        logging.error("Input file not found")
        return None

    return obtain_parser(tree)


def obtain_and_convert_parser_from_str(text: str) -> Optional[ArgumentParser]:
    """
        Function parses python source code stored in text and
        extracts argument parser

        Parameters
        ----------
        text : str
         variable containing target source code

        Returns
        -------
        Initialized argument parser, or None, in case error happened
        """
    tree = create_module_tree_from_str(text)

    return obtain_parser(tree)


def obtain_parser(tree: ast.Module) -> Optional[ArgumentParser]:
    try:
        actions, name, section_names, imported_names = \
            get_parser_init_and_actions(tree)

        actions, unknown_names = \
            initialize_variables_in_module(tree, name,
                                           section_names, actions,
                                           imported_names)

        result_module = handle_local_module_names(actions, unknown_names)
    except ArgumentParsingDiscoveryError as e:
        logging.error(e)
        return None

    # result_module.body.pop(2)
    ast.fix_missing_locations(result_module)
    compiled_module = compile(result_module, filename="<parser>", mode="exec")
    namespace = {}
    try:
        print(ast.unparse(result_module))

        exec(compiled_module, namespace)
    except Exception as e:
        print(e)
        print(ast.unparse(result_module))
        logging.error("Parser couldn't be extracted")
        return None
    return namespace[name]


def extract_useful_info_from_parser(parser: DecoyParser,
                                    data_inputs: Dict[str, str],
                                    reserved_names: Set[str]) \
        -> Tuple[List[ParamInfo], Dict[str, str]]:
    """
    Converts extracted argument parser object into tuple of Param info objects

    It also renames parameters whose names are equal to
    name in 'reserved_names' set

    containing parsed out data of arguments, and dictionary mapping new names
    of these arguments, compatible with galaxy wrappers, to the old names

    Parameters
    ----------
    parser : ArgumentParser
     extracted argument parser, initialized with possible command
     line arguments of target tool
    data_inputs : Dict[str, str]
     dictionary mapping names of arguments which should be interpreted
     as paths to input datasets
    reserved_names : set of names reserved by Galaxy, parameters whose name
     are equal to one of the names in reserved names must be renamed

    Returns
    -------
    Parsed out data of arguments, and dictionary mapping new names
    of these arguments, compatible with galaxy wrappers, to the old names
    """
    params = []
    name_map = dict()
    section_map = dict()

    for action in parser.report_arguments_and_groups():
        name = action.argument.lstrip("-").replace("-", "_")

        argument = action.argument
        type_ = _determine_type(data_inputs, name,
                                action.kwargs.get("type", str))
        nargs = _determine_nargs(action.kwargs.get("nargs", None))

        update_name_map(name, name_map, reserved_names)
        # these actions are of hidden type, and they contain the container
        # field. This field contains the information about their groups.
        # currently, only two level hierarchies are supported
        section = action.scope.name
        update_name_map(section, section_map, reserved_names)
        default_val = action.kwargs.get("default", None)

        help_ = action.kwargs.get("help", None)
        optional = action.kwargs.get("required", False)
        is_repeat = action.action == "APPEND" or nargs > 1
        is_select = action.kwargs.get("choices", None) is not None
        choices = action.kwargs.get("choices", None)
        is_flag = action.action == "STORE_TRUE"

        custom_attributes = _determine_custom_attributes(type_, nargs)

        params.append(ParamInfo(type_, name_map[name], argument, name,
                                section_map[section], section,
                                default_val, custom_attributes, nargs, help_,
                                optional, is_repeat, is_select, choices,
                                is_flag))

    return params, {new: old for old, new in name_map.items()}


def update_name_map(name: str, name_map: Dict[str, str],
                    reserved_names: Set[str]):
    """
    Updates names that are equal to one of the values in 'reserved_names'

    Parameters
    ----------
    name :
    name_map :
    reserved_names :
    """
    old_name = name
    while name.lower() in reserved_names:
        name = name + str(uuid.uuid4())[:4]

    name_map[old_name] = name


def _determine_type(data_inputs: Dict[str, str], name: str, type_):
    if type_ is None or type_ == bool:
        type_ = "boolean"
    elif type_ == int:
        type_ = "integer"
    elif type_ == float:
        type_ = "float"
    elif type_ == str:
        if name in data_inputs:
            type_ = "data"
        else:
            type_ = "text"
    else:
        type_ = f"{LINTER_MAGIC} argument uses complex type," \
                f" it's type cannot be determined"
    return type_


def _determine_nargs(nargs: Union[str, int, None]) -> Union[float, int]:
    if type(nargs) == str:
        if nargs == "?":
            return 1
        return math.inf
    if nargs is None:
        return 0
    return int(nargs)


def _determine_custom_attributes(type_: str, nargs: int) -> List[
    Tuple[str, str]]:
    flags: List[Tuple[str, str]] = []

    if type_ == "data" and nargs > 1:
        flags.append(("multiple", "true"))

    return flags
