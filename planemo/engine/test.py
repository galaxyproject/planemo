import os

from planemo.engine import engine_context
from planemo.galaxy import galaxy_config
from planemo.galaxy import galaxy_serve
from planemo.galaxy.api import (
    DEFAULT_ADMIN_API_KEY
)
from planemo.galaxy.config import _find_test_data
from planemo.galaxy.ephemeris_sleep import sleep
from planemo.galaxy.test import (
    handle_reports_and_summary,
    run_in_config,
)
from planemo.runnable import (
    flatten_to_single_artifacts,
    for_paths,
    RunnableType,
)


def test_runnables(ctx, runnables, original_paths=None, **kwds):
    """Return exit code indicating test or failure."""
    if kwds.get("serve"):
        if "galaxy" not in kwds["engine"]:
            raise ValueError("The serve option is only supported by Galaxy-based engines.")
        kwds["galaxy_url"] = kwds["galaxy_url"] or ''.join(("http://", kwds["host"], ":", str(kwds["port"])))
        kwds["galaxy_admin_key"] = kwds["galaxy_admin_key"] or DEFAULT_ADMIN_API_KEY
        pid = os.fork()
        if pid == 0:
            # wait for served Galaxy instance to start
            sleep(kwds["galaxy_url"], verbose=ctx.verbose, timeout=500)
            # then proceed to test against it
            kwds["engine"] = "external_galaxy"
        else:
            # serve Galaxy instance
            galaxy_serve(ctx, runnables, **kwds)
            exit(1)

    engine_type = kwds["engine"]
    test_engine_testable = {RunnableType.galaxy_tool, RunnableType.galaxy_datamanager, RunnableType.directory}
    enable_test_engines = any(r.type not in test_engine_testable for r in runnables)
    enable_test_engines = enable_test_engines or engine_type != "galaxy"

    if enable_test_engines:
        runnables = flatten_to_single_artifacts(runnables)  # the test engines cannot deal with directories
        ctx.vlog("Using test engine type %s" % engine_type)
        with engine_context(ctx, **kwds) as engine:
            test_data = engine.test(runnables)
            ctx.vlog("engine.test returning [%s]" % test_data)
            return_value = handle_reports_and_summary(ctx, test_data.structured_data, kwds=kwds)
    else:
        ctx.vlog("Running traditional Galaxy tool tests using run_tests.sh in Galaxy root %s" % engine_type)
        kwds["for_tests"] = True
        if kwds.get("update_test_data"):
            non_copied_runnables = for_paths(original_paths)
            kwds["test_data_target_dir"] = _find_test_data(non_copied_runnables, **kwds)
        with galaxy_config(ctx, runnables, **kwds) as config:
            return_value = run_in_config(ctx, config, **kwds)
    return return_value
