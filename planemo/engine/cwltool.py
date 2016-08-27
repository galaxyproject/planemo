"""Module contianing the :class:`CwlToolEngine` implementation of :class:`Engine`."""

from planemo import cwl
from planemo.runnable import RunnableType

from .interface import BaseEngine


class CwlToolEngine(BaseEngine):
    """An :class:`Engine` implementation backed by cwltool.

    More information on cwltool can be found at https://github.com/common-workflow-language/cwltool.
    """

    handled_runnable_types = [RunnableType.cwl_tool, RunnableType.cwl_workflow]

    def _run(self, runnable, job_path):
        """Run CWL job using cwltool."""
        path = runnable.path
        return cwl.run_cwltool(self._ctx, path, job_path, **self._kwds)


__all__ = ["CwlToolEngine"]
