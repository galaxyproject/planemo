import json

with open("tool_index.json", "r") as f:
    tool_index = json.load(f)

tool_map = {}

for tool in tool_index:
    tool_id = tool["id"]
    if "/" in tool_id:
        continue
    config_file = tool["config_file"]
    path = None
    if "lib/galaxy/tools/" in config_file:
        _, path = config_file.split("lib/galaxy/tools/", 1)
        path = "${model_tools_path}/%s" % path
    elif "tools/" in config_file:
        _, path = config_file.split("tools/", 1)

    if path:
        tool_map[tool_id] = path

as_python = "DISTRO_TOOLS_ID_TO_PATH = %s\n" % json.dumps(tool_map, indent=4)
as_python = as_python.replace(" \n", "\n")
with open("../planemo/galaxy/distro_tools.py", "w") as f:
    f.write("# file auto generated with scripts/tool_index_to_id_map.py\n")
    f.write(as_python)
