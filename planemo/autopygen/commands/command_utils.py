"""
Module used for creation and manipulation of template elements
are parts of template, separated from the rest by comments
with specific structure
Example of the comments:
    ## foo definition
    ... block itself ...
    ## end foo definition

Elements can be nested
"""
from typing import List

from planemo.autopygen.param_info import ParamInfo

SPACE = " "


class DefinitionNotFoundException(Exception):
    """
    Exception raised if part of template definition cannot be found
    """
    pass


def create_flag(variable: str, comment: str, depth: int, indent=3) -> str:
    """
    Function used to create a flag definition, wrapped in a comment

    Parameters
    ----------
    variable : str
        name of variable, containing $ at the beginning
    comment : str
        wrapping comment
    depth : int
     integer, used to set the depth of the current element.
     This value is used to indent the block properly
    indent : int
      default value for size of the block indent
    """
    return f"{depth * indent * SPACE}## FLAG {comment}\n" \
           f"{depth * indent * SPACE}{variable}\n" \
           f"{depth * indent * SPACE}## end FLAG {comment}\n"


def create_element_with_body(kind: str, head: str,
                             body: List[str], comment: str,
                             depth: int,
                             indent: int = 3,
                             body_indented: bool = True) -> str:
    """
    Function used to create block of template, like if or loop

    Parameters
    ----------
    kind : str
     string defining what kind of element is created, for example if or for
     (loop)
    head : str
     body of block header, for example predicate of condition, or the
     body of loop
    body : str
     body of the block, can be another element
    comment : str
     comment, used to set the start and end of the block
    depth : int
     integer, used to set the depth of the current element.
     This value is used to indent the block properly
    indent : int
      default value for size of the block indent
    body_indented : bool
      option that defines whether body is already correctly indented

    Returns
    -------
    string containing the created template element
    """

    result = (f"{depth * indent * SPACE}## {comment}\n"
              f"{depth * indent * SPACE}#{kind} {head}:\n")

    body_indent = ""
    if body:
        if not body_indented:
            body_indent = (depth + 1) * indent * SPACE
            body[0] = f"{body_indent}{body[0]}"
        result += ("\n" + body_indent).join(body) + "\n"

    result += f"{depth * indent * SPACE}#end {kind}\n"
    result += f"{depth * indent * SPACE}## end {comment}\n"
    return result


def transform_param_info(info: ParamInfo, namespace: str, depth: int):
    if info.param_type.is_help or info.param_type.is_version:
        raise ParamTypeNotSupported("Transformation for these param types are not supported")

    name = info.name
    separator = "." if namespace else ""
    variable = f"${namespace}{separator}{name}"
    if not info.param_type.is_repeat:
        if info.param_type.is_flag:
            return create_flag(variable, f"{name} definition", depth)
        else:
            body_expression = create_body_expression(info, variable, depth + 1)
            return create_element_with_body("if", variable, [body_expression],
                                            f"{name} definition", depth)

    iteration_var = "$item"
    if info.param_type.is_flag:
        param = create_flag(iteration_var, f"{name} definition", depth + 1)
    else:
        body_expression = create_body_expression(info, iteration_var, depth + 2)
        param = create_element_with_body("if", iteration_var, [body_expression],
                                         f"{name} definition", depth + 1)

    head_expression = f"{iteration_var} in ${namespace}{separator}{info.name}"

    return create_element_with_body("for", head_expression,
                                    [param],
                                    f"{info.name} definition",
                                    depth)


# FIXME generating command like this assumes that parameters are added to argparse in the right order.
# Separating arguments into positional and non-positional is necessary,
# otherwise the command generator will not work correctly
def create_body_expression(info: ParamInfo, variable: str, depth: int, indentation: int = 3) -> str:
    stripped_arg = info.argument.lstrip("-")
    indentation = SPACE * depth * indentation
    if stripped_arg == info.argument:
        return f"{indentation}{variable}"
    return f"{indentation}{info.argument} {variable}"


class ParamTypeNotSupported(Exception):
    pass
