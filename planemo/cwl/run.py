"""Module defines a planemo abstraction around running cwltool.

cwltool is an executable Python script and library mostly maintained by
Peter Amstutz and serves the reference implementation for the CWL.
It can be found at https://github.com/common-workflow-language/cwltool,
"""
import json
import tempfile

from galaxy.tools.cwl.cwltool_deps import (
    ensure_cwltool_available,
    main,
)

from planemo.io import error, real_io
from planemo.runnable import (
    ErrorRunResponse,
    SuccessfulRunResponse,
)

JSON_PARSE_ERROR_MESSAGE = ("Failed to parse JSON from cwltool output [%s] "
                            "in file [%s]. cwltool logs [%s].")


class CwlToolRunResponse(SuccessfulRunResponse):
    """Describe the resut of a cwltool invocation."""

    def __init__(self, log, cwl_command_state=None, outputs=None):
        self._log = log
        self._cwl_command_state = cwl_command_state
        self._outputs = outputs

    @property
    def log(self):
        return self._log

    @property
    def job_info(self):
        return None

    @property
    def cwl_command_state(self):
        if self._cwl_command_state is None:
            message = "Can only call cwl_command_state if running conformance_test."
            raise NotImplementedError(message)

        return self._cwl_command_state

    @property
    def outputs_dict(self):
        if self._outputs is None:
            message = "Can not call outputs if running conformance_test."
            raise NotImplementedError(message)
        return self._outputs


def run_cwltool(ctx, path, job_path, **kwds):
    """Translate planemo kwds to cwltool kwds and run cwltool main function."""
    ensure_cwltool_available()

    args = []
    conformance_test = kwds.get("conformance_test", False)
    if conformance_test:
        args.append("--conformance-test")
    if ctx.verbose:
        args.append("--verbose")
    output_directory = kwds.get("output_directory", None)
    if output_directory:
        args.append("--outdir")
        args.append(output_directory)
    if kwds.get("no_container", False):
        args.append("--no-container")

    args.extend([path, job_path])
    ctx.vlog("Calling cwltool with arguments %s" % args)
    with tempfile.NamedTemporaryFile() as tmp_stdout, \
            tempfile.NamedTemporaryFile() as tmp_stderr:
        # cwltool passes sys.stderr to subprocess.Popen - ensure it has
        # and actual fileno.
        with real_io():
            ret_code = main.main(
                args,
                stdout=tmp_stdout,
                stderr=tmp_stderr
            )
        tmp_stdout.flush()
        tmp_stderr.flush()
        with open(tmp_stderr.name, "r") as stderr_f:
            log = stderr_f.read()
            ctx.vlog("cwltool log output [%s]" % log)
        with open(tmp_stdout.name, "r") as stdout_f:
            try:
                result = json.load(stdout_f)
            except ValueError:
                message = JSON_PARSE_ERROR_MESSAGE % (
                    open(tmp_stdout.name, "r").read(),
                    tmp_stdout.name,
                    log
                )
                error(message)
                raise Exception(message)

        if ret_code != 0:
            return ErrorRunResponse("Error running cwltool", log=log)
        if conformance_test:
            cwl_command_state = result
            outputs = None
        else:
            cwl_command_state = None
            outputs = result
    return CwlToolRunResponse(
        log,
        cwl_command_state=cwl_command_state,
        outputs=outputs,
    )


__all__ = [
    "run_cwltool",
]
