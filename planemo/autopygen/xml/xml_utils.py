from typing import Optional

from planemo.autopygen.param_info import ParamInfo

SPACE = " "


def options(param_info: ParamInfo, depth: int):
    opts = []
    for option in param_info.choices:
        opts.append(
            formatted_xml_elem("option", {"value": option}, depth + 1, option.capitalize(), inline=True))

    return param(param_info, depth, "".join(opts).rstrip())


def repeat(param_info: ParamInfo, depth: int):
    attributes = {
        "name": param_info.name + "_repeat",
        "title": param_info.name + "_repeat",
        "min": param_info.custom_attributes.get("min", None),
        "max": param_info.custom_attributes.get("max", None),
        "default": param_info.custom_attributes.get("default", None),
    }

    names = list(attributes.keys())
    for name in names:
        if attributes[name] is None:
            attributes.pop(name)

    if param_info.param_type.is_selection:
        inner = options(param_info, depth + 1)
    else:
        inner = param(param_info, depth + 1)

    return formatted_xml_elem("repeat", attributes, depth, inner)


def param(param_info: ParamInfo, depth: int, body: Optional[str] = None):
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

    return formatted_xml_elem("param", attributes, depth, body)


def formatted_xml_elem(name, attributes, depth, body: Optional[str] = None, indent: int = 4,
                       body_indented: bool = True, inline: bool = False):
    left_indent: str = depth * indent * SPACE
    right_indent: str = left_indent
    if body:
        if not body_indented:
            body = f"{left_indent}{body}"

        body_newlines = "\n"
        if inline:
            body_newlines = ""
            right_indent = ""

        return f'{left_indent}<{name} {attributes_to_str(attributes)}>' \
               f'{body_newlines}{body.rstrip()}{body_newlines}' \
               f'{right_indent}</{name}>\n'

    return f'{left_indent}<{name} {attributes_to_str(attributes)} />\n'


def attributes_to_str(attributes):
    result = []
    for key, value in attributes.items():
        result.append(f'{key}="{value}"')
    return " ".join(result)
