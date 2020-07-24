from planemo.engine import (
    engine_context,
)
from planemo.galaxy import galaxy_config
from planemo.galaxy.config import _find_test_data
from planemo.galaxy.test import (
    handle_reports_and_summary,
    run_in_config,
)
from planemo.runnable import (
    for_paths,
    RunnableType,
)


def test_runnables(ctx, runnables, original_paths=None, **kwds):
    """Return exit code indicating test or failure."""
    engine_type = kwds["engine"]
    test_engine_testable = {RunnableType.galaxy_tool, RunnableType.galaxy_datamanager, RunnableType.directory}
    enable_test_engines = any(r.type not in test_engine_testable for r in runnables)
    enable_test_engines = enable_test_engines or engine_type != "galaxy"
    if enable_test_engines:
        ctx.vlog("Using test engine type %s" % engine_type)
        with engine_context(ctx, **kwds) as engine:
            test_data = engine.test(runnables)
            ctx.vlog("engine.test returning [%s]" % test_data)
            return_value = handle_reports_and_summary(ctx, test_data.structured_data, kwds=kwds)
    else:
        ctx.vlog("Running traditional Galaxy tool tests using run_tests.sh in Galaxy root %s" % engine_type)
        kwds["for_tests"] = True
        if kwds.get('update_test_data'):
            non_copied_runnables = for_paths(original_paths)
            kwds['test_data_target_dir'] = _find_test_data(non_copied_runnables, **kwds)
        with galaxy_config(ctx, runnables, **kwds) as config:
            return_value = run_in_config(ctx, config, **kwds)
    return return_value
