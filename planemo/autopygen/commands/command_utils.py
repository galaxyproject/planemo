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


def create_element_with_body(kind: str, head: str,
                             body: List[str], comment: str,
                             depth: int,
                             indent: int = 3) -> str:
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

    Returns
    -------
    string containing the created template element
    """
    result = (f"{depth * indent * SPACE}## {comment}\n"
              f"{depth * indent * SPACE}#{kind} {head}:\n{indent * SPACE}")

    result += ("\n" + indent * SPACE).join(body) + "\n"

    result += f"{depth * indent * SPACE}#end {kind}\n"
    result += f"{depth * indent * SPACE}## end {comment}\n"
    return result


def transform_param_info(info: ParamInfo, sections: str, depth: int):
    body_expression = info.argument
    name = info.name
    variable = f"${sections}.{name}"

    if info.type != "boolean":
        body_expression += f" {variable}"

    if info.is_repeat:
        iteration_var = "$item"
        param = create_element_with_body("if", iteration_var, [f"{info.argument} {iteration_var}"],
                                         f"{name} definition", depth + 1)
        head_expression = (f"{iteration_var}"
                           f" in ${sections}."
                           f"{info.name}")

        return create_element_with_body("for", head_expression,
                                        [param], f"{info.name} definition",
                                        depth)

    return create_element_with_body("if", variable, [body_expression],
                                    f"{name} definition", depth)
