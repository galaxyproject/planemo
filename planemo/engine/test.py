from planemo.engine import engine_context
from planemo.galaxy.test import handle_reports_and_summary
from planemo.runnable import RunnableType


def test_runnables(ctx, runnables, original_paths=None, **kwds):
    """Return exit code indicating test or failure."""
    engine_type = kwds["engine"]
    test_engine_testable = {RunnableType.galaxy_tool, RunnableType.galaxy_datamanager, RunnableType.directory}
    enable_test_engines = any(r.type not in test_engine_testable for r in runnables)
    enable_test_engines = enable_test_engines or engine_type != "galaxy"
    ctx.vlog(f"Using test engine type {engine_type}")
    with engine_context(ctx, **kwds) as engine:
        test_data = engine.test(runnables, test_timeout=kwds.get("test_timeout"))
        ctx.vlog(f"engine.test returning [{test_data}]")
        return handle_reports_and_summary(ctx, test_data.structured_data, kwds=kwds)
