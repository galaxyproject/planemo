"""Module contianing the :class:`GalaxyEngine` implementation of :class:`Engine`."""

import abc
import contextlib

from galaxy.tool_util.verify import interactor

from planemo.galaxy.activity import (
    execute,
    execute_rerun,
)
from planemo.galaxy.config import external_galaxy_config
from planemo.galaxy.serve import serve_daemon
from planemo.runnable import RunnableType
from .interface import BaseEngine


class GalaxyEngine(BaseEngine, metaclass=abc.ABCMeta):
    """An :class:`Engine` implementation backed by a managed Galaxy.

    More information on Galaxy can be found at http://galaxyproject.org/.
    """

    handled_runnable_types = [
        RunnableType.cwl_tool,
        RunnableType.cwl_workflow,
        RunnableType.galaxy_workflow,
        RunnableType.galaxy_tool,
        RunnableType.galaxy_datamanager,
    ]

    def _run(self, runnables, job_paths):
        """Run CWL job in Galaxy."""
        results = []
        for runnable, job_path in zip(runnables, job_paths):
            self._ctx.vlog(f"Serving artifact [{runnable}] with Galaxy.")
            with self.ensure_runnables_served([runnable]) as config:
                self._ctx.vlog("Running job path [%s]" % job_path)
                if self._ctx.verbose:
                    self._ctx.log("Running Galaxy with API configuration [%s]" % config.user_api_config)
                run_response = execute(self._ctx, config, runnable, job_path, **self._kwds)
                results.append(run_response)

        return results

    @abc.abstractmethod
    def ensure_runnables_served(self, runnables):
        """Use a context manager and describe Galaxy instance with runnables being served."""

    def _run_test_cases(self, test_cases):
        test_results = []
        file_based_test_cases = []
        embedded_test_cases = []
        # TODO: unify interface so we don't need to split test cases
        for test_case in test_cases:
            if hasattr(test_cases, "job_path"):
                # Simple file-based job path.
                file_based_test_cases.append(test_case)
            else:
                embedded_test_cases.append(test_case)
        if file_based_test_cases:
            test_results.extend(super()._run_test_cases(file_based_test_cases))
        if embedded_test_cases:
            runnables = [test_case.runnable for test_case in embedded_test_cases]
            with self.ensure_runnables_served(runnables) as config:
                for test_case in embedded_test_cases:
                    galaxy_interactor_kwds = {
                        "galaxy_url": config.galaxy_url,
                        "master_api_key": config.master_api_key,
                        "api_key": config.user_api_key,
                        "keep_outputs_dir": "",  # TODO: this...
                    }
                    tool_id = test_case.tool_id
                    test_index = test_case.test_index
                    tool_version = test_case.tool_version
                    galaxy_interactor = interactor.GalaxyInteractorApi(**galaxy_interactor_kwds)

                    def _register_job_data(job_data):
                        test_results.append(
                            {
                                "id": tool_id + "-" + str(test_index),
                                "has_data": True,
                                "data": job_data,
                            }
                        )

                    verbose = self._ctx.verbose
                    try:
                        if verbose:
                            # TODO: this is pretty hacky, it'd be better to send a stream
                            # and capture the output information somehow.
                            interactor.VERBOSE_GALAXY_ERRORS = True

                        interactor.verify_tool(
                            tool_id,
                            galaxy_interactor,
                            test_index=test_index,
                            tool_version=tool_version,
                            register_job_data=_register_job_data,
                            quiet=not verbose,
                        )
                    except Exception:
                        pass

        return test_results


class LocalManagedGalaxyEngine(GalaxyEngine):
    """An :class:`Engine` implementation backed by a managed Galaxy.

    More information on Galaxy can be found at http://galaxyproject.org/.
    """

    @contextlib.contextmanager
    def ensure_runnables_served(self, runnables):
        # TODO: define an interface for this - not everything in config would make sense for a
        # pre-existing Galaxy interface.
        with serve_daemon(self._ctx, runnables, **self._serve_kwds()) as config:
            yield config

    def _serve_kwds(self):
        return self._kwds.copy()


class DockerizedManagedGalaxyEngine(LocalManagedGalaxyEngine):
    """An :class:`Engine` implementation backed by Galaxy running in Docker.

    More information on Galaxy can be found at http://galaxyproject.org/.
    """

    def _serve_kwds(self):
        serve_kwds = self._kwds.copy()
        serve_kwds["dockerize"] = True
        return serve_kwds


class ExternalGalaxyEngine(GalaxyEngine):
    """An :class:`Engine` implementation backed by an external Galaxy instance."""

    @contextlib.contextmanager
    def ensure_runnables_served(self, runnables):
        # TODO: ensure tools are available
        with external_galaxy_config(self._ctx, runnables, **self._kwds) as config:
            config.install_workflows()
            yield config

    def rerun(self, ctx, rerunnable, **kwds):
        with self.ensure_runnables_served([]) as config:
            rerun_response = execute_rerun(ctx, config, rerunnable, **kwds)
            return rerun_response


__all__ = (
    "DockerizedManagedGalaxyEngine",
    "ExternalGalaxyEngine",
    "LocalManagedGalaxyEngine",
)
