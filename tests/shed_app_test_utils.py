import contextlib
import multiprocessing
import shutil
from tempfile import mkdtemp
from typing import NamedTuple

from werkzeug.serving import run_simple

from planemo import network_util
from .shed_app import (
    app,
    InMemoryShedDataModel,
)

DEFAULT_OP_TIMEOUT = 2


def mock_model(directory):
    return (
        InMemoryShedDataModel(directory)
        .add_category("c1", "Text Manipulation")
        .add_category("c2", "Sequence Analysis")
        .add_category("c3", "Tool Dependency Packages")
        .add_repository(
            "r1",
            name="test_repo_1",
            owner="iuc",
        )
    )


def setup_mock_shed():
    port = network_util.get_free_port()
    directory = mkdtemp()
    model = mock_model(directory)

    def run():
        app.debug = True
        app.config["model"] = model
        run_simple("localhost", port, app, use_reloader=False, use_debugger=True)

    p = multiprocessing.Process(target=run)
    p.start()
    network_util.wait_net_service("localhost", port, DEFAULT_OP_TIMEOUT)
    return MockShed(f"http://localhost:{port}", directory, p, model)


@contextlib.contextmanager
def mock_shed():
    mock_shed_obj = None
    try:
        mock_shed_obj = setup_mock_shed()
        yield mock_shed_obj
    finally:
        if mock_shed_obj is not None:
            mock_shed_obj.shutdown()


class MockShed(NamedTuple):
    url: str
    directory: str
    process: multiprocessing.Process
    model: InMemoryShedDataModel

    def shutdown(self):
        self.process.terminate()
        shutil.rmtree(self.directory)


__all__ = (
    "setup_mock_shed",
    "mock_shed",
)
