"""Module contianing the :class:`GalaxyEngine` implementation of :class:`Engine`."""
import contextlib

from .interface import BaseEngine
from planemo.runnable import RunnableType

from planemo.galaxy.serve import serve_daemon
from planemo.galaxy.activity import execute


class GalaxyEngine(BaseEngine):
    """An :class:`Engine` implementation backed by Galaxy.

    More information on Galaxy can be found at http://galaxyproject.org/.
    """

    handled_runnable_types = [
        RunnableType.cwl_tool,
        RunnableType.galaxy_workflow,
        RunnableType.galaxy_tool
    ]

    def _run(self, runnable, job_path):
        """Run CWL job in Galaxy."""
        self._ctx.vlog("Serving artifact [%s] with Galaxy." % (runnable,))
        with self._serve([runnable]) as config:
            self._ctx.vlog("Running job path [%s]" % job_path)
            run_response = execute(config, runnable, job_path, **self._kwds)

        return run_response

    @contextlib.contextmanager
    def _serve(self, runnables):
        with serve_daemon(self._ctx, runnables, **self._kwds) as config:
            yield config


__all__ = ["GalaxyEngine"]
