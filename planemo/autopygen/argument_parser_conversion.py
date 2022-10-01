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
from planemo.autopygen.param_info import ParamInfo, ParamTypeFlags
from planemo.autopygen.source_file_parsing.constants import LINTER_MAGIC
from planemo.autopygen.source_file_parsing.decoy_parser import DecoyParser
from planemo.autopygen.source_file_parsing.local_module_parsing import handle_local_module_names
from planemo.autopygen.source_file_parsing.parser_discovery_and_init import \
    get_parser_init_and_actions
from planemo.autopygen.source_file_parsing.parsing_commons import \
    create_module_tree_from_path, create_module_tree_from_str
from planemo.autopygen.source_file_parsing.parsing_exceptions import \
    ArgumentParsingDiscoveryError
from planemo.autopygen.source_file_parsing.unknown_names_discovery import initialize_variables_in_module


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
        actions, name, known_names = \
            get_parser_init_and_actions(tree)

        actions, unknown_names = \
            initialize_variables_in_module(tree, name, actions,
                                           known_names)

        result_module = handle_local_module_names(actions, unknown_names)
    except ArgumentParsingDiscoveryError as e:
        logging.error(e)
        return None

    ast.fix_missing_locations(result_module)
    compiled_module = compile(result_module, filename="<parser>", mode="exec")
    namespace = {}
    try:
        exec(compiled_module, namespace)
    except Exception as e:
        print(e)
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

    if param_info.param_type.is_selection:
        inner = _options(param_info)
    else:
        inner = _default(param_info)

    return f'<repeat {attr_str}>{inner}</repeat>\n'


def _default(param_info: ParamInfo, body: Optional[str] = None):
    attributes = {
        "argument": param_info.argument,
        "type": param_info.type,
        "format": param_info.format,
    }

    if param_info.param_type.is_flag:
        attributes["truevalue"] = param_info.argument
        attributes["falsevalue"] = ""
        attributes["checked"] = "false"

    attributes["optional"] = str(param_info.optional).lower()
    attributes["label"] = param_info.label

    if param_info.help is not None:
        attributes["help"] = param_info.help
    # this might seem weird, but it is done like this for correct order
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


def _action_to_param(param_info: ParamInfo):
    if param_info.param_type.is_repeat:
        return _repeat(param_info)
    elif param_info.param_type.is_selection:
        return _options(param_info)

    return _default(param_info)


def generate_inputs_from_section(section: DecoyParser.Section,
                                 data_inputs: Dict[str, str],
                                 reserved_names: Set[str],
                                 name_map: Dict[str, str],
                                 section_map: Dict[str, str]):
    sub_actions = []
    for action in section.actions:
        param_info = obtain_param_info(action, data_inputs, reserved_names,
                                       name_map, section_map)
        if param_info.param_type.is_help or param_info.param_type.is_version:
            continue

        sub_actions.append(_action_to_param(param_info))

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

        if param_info.param_type.is_help or param_info.param_type.is_version:
            continue

        yield transform_param_info(param_info, nesting,
                                   depth)

    for subsection in section.subsections:
        sec_name = re.sub("[/\\-* ()]", "_", section_map.get(subsection.name, subsection.name)).lower()

        yield from _command_recursion(subsection,
                                      data_inputs,
                                      reserved_names,
                                      name_map,
                                      section_map,
                                      nesting + f".{sec_name}",
                                      depth + 1)


def _sub_parsers_conditionals(parsers: List[DecoyParser], index: int,
                              data_inputs: Dict[str, str],
                              reserved_names: Set[str],
                              name_map: Dict[str, str],
                              section_map: Dict[str, str]):
    result = []
    conditional_options = []
    is_first = True
    for parser in parsers:
        conditional_options.append(f'<option value="{parser.name}" selected="{str(is_first).lower()}">'
                                   f'{parser.name}'
                                   f'</option>')
        if is_first:
            is_first = False

        inner = inputs_from_decoy(parser, data_inputs,
                                  reserved_names, name_map, section_map)
        result.append(f'<when value="{parser.name}">\n{inner}\n</when>')

    conditional_options_joined = "\n".join(conditional_options)
    result_joined = "\n".join(result)

    return f'<conditional name="subparsers{index}">' \
           f'<param name="subparser_selector" type="select" label="">' \
           f'{conditional_options_joined}' \
           f'</param>' \
           f'{result_joined}' \
           f'</conditional>'


def inputs_from_decoy(parser: DecoyParser, data_inputs: Dict[str, str],
                      reserved_names: Set[str],
                      name_map: Dict[str, str],
                      section_map: Dict[str, str]) -> str:
    result = [generate_inputs_from_section(parser.default_section, data_inputs,
                                           reserved_names,
                                           name_map,
                                           section_map)]

    for index, sub_parsers in enumerate(parser.sub_parsers):
        result.append(_sub_parsers_conditionals(sub_parsers.parsers,
                                                index,
                                                data_inputs,
                                                reserved_names,
                                                name_map,
                                                section_map))

    return "\n".join(result)


def command_from_decoy(parser: DecoyParser, data_inputs: Dict[str, str],
                       reserved_names: Set[str],
                       name_map: Dict[str, str],
                       section_map: Dict[str, str]) -> str:
    """
    The function generates commands from decoy parser. It requires name and section maps that have already been
    initialised and contain correct mapping
    """
    sec_name = re.sub("[/\\-* ()]", "_", section_map[parser.default_section.name]).lower()
    result = ["\n".join(_command_recursion(parser.default_section, data_inputs,
                                           reserved_names,
                                           name_map,
                                           section_map, sec_name, 0))]

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

    nargs = _determine_nargs(action.kwargs.get("nargs", None))

    update_name_map(name, name_map, reserved_names)

    section = action.scope.name
    update_name_map(section, section_map, reserved_names)
    default_val = action.kwargs.get("default", None)

    help_ = action.kwargs.get("help", None)
    optional = action.kwargs.get("required", False)

    param_type = _determine_param_type(action, nargs)
    type_ = _determine_type(action, data_inputs, name,
                            action.kwargs.get("type", str), param_type)
    choices = action.kwargs.get("choices", None)
    custom_attributes = _determine_custom_attributes(type_, nargs)

    return ParamInfo(param_type, type_, name_map[name], argument, name,
                     section_map[section], section,
                     default_val, custom_attributes, nargs, help_,
                     optional, choices)


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


def _determine_param_type(action: DecoyParser.Action, nargs: Union[int, float]):
    param_flags = ParamTypeFlags()

    param_flags.is_selection = action.kwargs.get("choices", None) is not None
    param_flags.is_repeat = action.action in ["APPEND", "APPEND_CONST", "COUNT"] or nargs > 1
    param_flags.is_flag = action.action in ["STORE_TRUE", "STORE_CONST", "STORE_FALSE"]
    param_flags.is_extend = action.action == "EXTEND"
    param_flags.is_help = action.action == "HELP"
    param_flags.is_version = action.action == "VERSION"

    return param_flags


def _determine_type(action: DecoyParser.Action, data_inputs: Dict[str, str],
                    name: str, type_, param_type_flags: ParamTypeFlags):
    if action.kwargs.get("choices", []):
        return "select"

    if type_ is None or type_ == bool or param_type_flags.is_flag:
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
