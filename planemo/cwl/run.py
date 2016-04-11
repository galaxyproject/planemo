from planemo.io import conditionally_captured_io
from planemo.galaxy.serve import serve_daemon
from .client import run_cwl_tool
from planemo import io

try:
    from cwltool.main import main
except (ImportError, SyntaxError):
    # TODO: Add Python 2.6 and 3.X support to cwltool so this SyntaxError
    # check isn't needed.
    main = None


def run_galaxy(ctx, path, job_path, **kwds):
    kwds["cwl"] = True
    conformance_test = kwds.get("conformance_test", False)
    with conditionally_captured_io(conformance_test):
        with serve_daemon(ctx, [path], **kwds) as config:
            try:
                cwl_run = run_cwl_tool(path, job_path, config, **kwds)
            except Exception:
                io.warn("Problem running cwl tool...")
                print(config.log_contents)
                raise

    print(cwl_run.cwl_command_state)
    return 0


def run_cwltool(ctx, path, job_path, **kwds):
    if main is None:
        raise Exception("cwltool dependency not found.")

    args = []
    if kwds.get("conformance_test", False):
        args.append("--conformance-test")
    if ctx.verbose:
        args.append("--verbose")

    args.extend([path, job_path])
    ctx.vlog("Calling cwltool with arguments %s" % args)
    return main(args)
