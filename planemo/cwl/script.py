import copy
import os

import schema_salad
import cwltool
from cwltool.main import load_tool
from cwltool.process import (
    checkRequirements,
    shortname,
)
from cwl2script.cwl2script import (
    supportedProcessRequirements,
    generateScriptForWorkflow,
    generateScriptForTool,
)


def to_script(ctx, path, job_path, **kwds):
    uri = "file://" + os.path.abspath(job_path)

    loader = schema_salad.ref_resolver.Loader({
        "@base": uri,
        "path": {
            "@type": "@id"
        }
    })

    job, _ = loader.resolve_ref(uri)

    t = load_tool(path, False, False, cwltool.workflow.defaultMakeTool, True)

    if type(t) == int:
        return t

    checkRequirements(t.tool, supportedProcessRequirements)

    for inp in t.tool["inputs"]:
        if shortname(inp["id"]) in job:
            pass
        elif shortname(inp["id"]) not in job and "default" in inp:
            job[shortname(inp["id"])] = copy.copy(inp["default"])
        elif shortname(inp["id"]) not in job and inp["type"][0] == "null":
            pass
        else:
            raise Exception("Missing inputs `%s`" % shortname(inp["id"]))

    if not kwds.get("basedir", None):
        kwds["basedir"] = os.path.dirname(os.path.abspath(job_path))

    outdir = kwds.get("outdir")

    if t.tool["class"] == "Workflow":
        print(generateScriptForWorkflow(t, job, outdir))
    elif t.tool["class"] == "CommandLineTool":
        print(generateScriptForTool(t, job, outdir))

    return 0
