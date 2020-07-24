import copy
import os

try:
    import schema_salad
except (ImportError, SyntaxError):
    schema_salad = None

try:
    import cwltool
except (ImportError, SyntaxError):
    cwltool = None

try:
    from cwltool.load_tool import load_tool
except (ImportError, SyntaxError):
    load_tool = None

try:
    from cwltool import process
except (ImportError, SyntaxError):
    process = None

try:
    from .cwl2script import cwl2script
except (ImportError, SyntaxError):
    cwl2script = None


def to_script(ctx, path, job_path, **kwds):
    if schema_salad is None:
        raise Exception("This functionality requires schema_salad.")

    if cwltool is None:
        raise Exception("This functionality requires cwltool.")

    uri = "file://" + os.path.abspath(job_path)

    loader = schema_salad.ref_resolver.Loader({
        "@base": uri,
        "path": {
            "@type": "@id"
        }
    })

    job, _ = loader.resolve_ref(uri)

    t = load_tool(path)

    if type(t) == int:
        return t

    process.checkRequirements(t.tool, cwl2script.supportedProcessRequirements)

    for inp in t.tool["inputs"]:
        if process.shortname(inp["id"]) in job:
            pass
        elif process.shortname(inp["id"]) not in job and "default" in inp:
            job[process.shortname(inp["id"])] = copy.copy(inp["default"])
        elif process.shortname(inp["id"]) not in job and inp["type"][0] == "null":
            pass
        else:
            raise Exception("Missing inputs `%s`" % process.shortname(inp["id"]))

    if not kwds.get("basedir", None):
        kwds["basedir"] = os.path.dirname(os.path.abspath(job_path))

    outdir = kwds.get("outdir")

    if t.tool["class"] == "Workflow":
        print(cwl2script.generateScriptForWorkflow(t, job, outdir))
    elif t.tool["class"] == "CommandLineTool":
        print(cwl2script.generateScriptForTool(t, job, outdir))

    return 0
