"""
Functions used to extract initialized ArgumentParser from target source code,
and to transform it into easily usable ParamInfo class containing
 all the necessary info for 'inputs' element initialization
"""
import ast
import logging
import math
import re
import uuid

from typing import Optional, Set, List, Dict, Tuple, Union
from argparse import ArgumentParser

from planemo.autopygen.commands.command_utils import transform_param_info
from planemo.autopygen.param_info import ParamInfo
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


def obtain_and_convert_parser(path: str) -> Optional[DecoyParser]:
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


def _options(param_info: ParamInfo):
    options = []
    for option in param_info.choices:
        options.append(
            f'<option value="{option}">{option.capitalize()}</option>\n')

    return _default(param_info, "\n".join(options).rstrip())


def _repeat(param_info: ParamInfo):
    attributes = {
        "name": param_info.name + "_repeat",
        "title": param_info.name + "_repeat",
        "min_reps": None,  # TODO add support for reps
        "max_reps": None,
        "default_reps": None,
    }

    names = list(attributes.keys())
    for name in names:
        if attributes[name] is None:
            attributes.pop(name)

    attr_str = attributes_to_str(attributes)

    if param_info.is_select:
        inner = _options(param_info)
    else:
        inner = _default(param_info)

    return f'<repeat {attr_str}>{inner}</repeat>\n'


def _default(param_info: ParamInfo, body: Optional[str] = None):
    attributes = {
        "argument": param_info.argument,
        "type": param_info.type,
        "format": param_info.format,
        "optional": str(param_info.optional).lower(),
        "label": param_info.label
    }

    if param_info.help is not None:
        attributes["help"] = param_info.help
    # this might seem weird but it is done like this for correct order
    # of attributes
    if param_info.format is None:
        attributes.pop("format")

    if body:
        return f'<param {attributes_to_str(attributes)}>\n{body}\n</param>\n'

    return f'<param {attributes_to_str(attributes)}/>\n'


def attributes_to_str(attributes):
    result = []
    for key, value in attributes.items():
        result.append(f'{key}="{value}"')
    return " ".join(result)


def _action_to_param(action: DecoyParser.Action,
                     data_inputs: Dict[str, str],
                     reserved_names: Set[str],
                     name_map: Dict[str, str],
                     section_map: Dict[str, str]):
    param_info = obtain_param_info(action, data_inputs, reserved_names,
                                   name_map, section_map)

    if param_info.is_repeat:
        return _repeat(param_info)
    elif param_info.is_select:
        return _options(param_info)

    return _default(param_info)


def generate_inputs_from_section(section: DecoyParser.Section,
                                 data_inputs: Dict[str, str],
                                 reserved_names: Set[str],
                                 name_map: Dict[str, str],
                                 section_map: Dict[str, str]):
    sub_actions = [_action_to_param(action,
                                    data_inputs,
                                    reserved_names,
                                    name_map,
                                    section_map) for action in section.actions]
    sub_sections = [
        generate_inputs_from_section(subsection,
                                     data_inputs,
                                     reserved_names,
                                     name_map,
                                     section_map) for subsection in
        section.subsections if subsection.actions]

    transformed_name = re.sub("[/\\-* ()]", "_", section_map.get(section.name, section.name)).lower()
    return f'<section name="{transformed_name}" title="{section.name}" expanded="true">\n' \
           f'{"".join(sub_actions)}' \
           f'{"".join(sub_sections)}' \
           f'</section>\n'


def _command_recursion(section: DecoyParser.Section,
                       data_inputs: Dict[str, str],
                       reserved_names: Set[str],
                       name_map: Dict[str, str],
                       section_map: Dict[str, str], nesting: str, depth: int):
    for action in section.actions:
        param_info = obtain_param_info(action, data_inputs, reserved_names,
                                       name_map, section_map)
        yield transform_param_info(param_info, nesting,
                                   depth)

    for subsection in section.subsections:
        sec_name = re.sub("[/\\-* ()]", "_", section_map[subsection.name]).lower()

        yield from _command_recursion(subsection,
                                      data_inputs,
                                      reserved_names,
                                      name_map,
                                      section_map,
                                      nesting + f".{sec_name}",
                                      depth + 1)


def inputs_from_decoy(parser: DecoyParser, data_inputs: Dict[str, str],
                      reserved_names: Set[str],
                      name_map: Dict[str, str],
                      section_map: Dict[str, str]):
    result = []

    for subsection in parser.default_section.subsections:
        if not subsection.actions and not subsection.subsections:
            continue
        result.append(generate_inputs_from_section(subsection, data_inputs,
                                                   reserved_names,
                                                   name_map,
                                                   section_map))

    for sub_parsers in parser.sub_parsers:
        for parser in sub_parsers.parsers:
            result.append(inputs_from_decoy(parser, data_inputs,
                                            reserved_names, name_map, section_map))

    return "\n".join(result)


def command_from_decoy(parser: DecoyParser, data_inputs: Dict[str, str],
                       reserved_names: Set[str],
                       name_map: Dict[str, str],
                       section_map: Dict[str, str]):
    """
    The function generates commands from decoy parser. It requires name and section maps that have already been
    initialised and contain correct mapping
    """
    result = []
    for subsection in parser.default_section.subsections:
        if not subsection.actions:
            continue
        sec_name = re.sub("[/\\-* ()]", "_", section_map[subsection.name]).lower()
        result.extend([command for command in
                       _command_recursion(subsection,
                                          data_inputs,
                                          reserved_names,
                                          name_map,
                                          section_map,
                                          sec_name,
                                          depth=0)])

        for sub_parsers in parser.sub_parsers:
            for parser in sub_parsers.parsers:
                result.append(command_from_decoy(parser, data_inputs,
                                                 reserved_names, name_map, section_map))

    return "\n".join(result)


def obtain_param_info(action: DecoyParser.Action,
                      data_inputs: Dict[str, str],
                      reserved_names: Set[str],
                      name_map: Dict[str, str],
                      section_map: Dict[str, str]):
    name = action.argument.lstrip("-").replace("-", "_")
    argument = action.argument
    type_ = _determine_type(action, data_inputs, name,
                            action.kwargs.get("type", str))
    nargs = _determine_nargs(action.kwargs.get("nargs", None))

    update_name_map(name, name_map, reserved_names)

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

    return ParamInfo(type_, name_map[name], argument, name,
                     section_map[section], section,
                     default_val, custom_attributes, nargs, help_,
                     optional, is_repeat, is_select, choices,
                     is_flag)


def extract_useful_info_from_parser(parser: DecoyParser,
                                    data_inputs: Dict[str, str],
                                    reserved_names: Set[str]) \
    -> Tuple[List[ParamInfo], Dict[str, str]]:
    """
    !!! SOON TO BE DEPRECATED !!!
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
        param_info = obtain_param_info(action, data_inputs, reserved_names,
                                       name_map, section_map)
        params.append(param_info)

    return params, {new: old for old, new in name_map.items()}


def update_name_map(name: str, name_map: Dict[str, str],
                    reserved_names: Set[str]):
    """
    Updates names that are equal to one of the values in 'reserved_names'
    Doesn't change names already present in the mapping

    Parameters
    ----------
    name :
    name_map :
    reserved_names :
    """
    old_name = name

    if old_name in name_map:
        return

    while name.lower() in reserved_names:
        name = name + str(uuid.uuid4())[:4]

    name_map[old_name] = name


def _determine_type(action: DecoyParser.Action, data_inputs: Dict[str, str],
                    name: str, type_):
    if action.kwargs.get("choices", []):
        return "select"

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


def _determine_custom_attributes(type_: str, nargs: int) -> List[Tuple[str, str]]:
    flags: List[Tuple[str, str]] = []

    if type_ == "data" and nargs > 1:
        flags.append(("multiple", "true"))

    return flags
