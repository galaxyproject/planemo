"""Tool Shed repository installs must give up on a Galaxy that never answers."""

import contextlib
import json
import threading
from http.server import (
    BaseHTTPRequestHandler,
    ThreadingHTTPServer,
)

import pytest
import requests

from planemo.galaxy.api import gi
from planemo.galaxy.config import BaseGalaxyConfig
from planemo.galaxy.workflows import (
    _install_shed_repos_from_tools_info,
    FAILED_REPOSITORIES_MESSAGE,
)

REPOSITORY = {
    "name": "rgi",
    "owner": "card",
    "changeset_revision": "84bd24ac33fd",
    "tool_shed_url": "https://toolshed.g2.bx.psu.edu",
}


@contextlib.contextmanager
def _fake_galaxy(hang_get=False, hang_names=("rgi",)):
    """A Galaxy whose install requests for ``hang_names`` (and optionally every GET) never return."""
    release = threading.Event()

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            if hang_get:
                release.wait()
            else:
                self._json([])

        def do_POST(self):
            payload = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            if payload["name"] in hang_names:
                release.wait()
            else:
                self._json([])

        def _json(self, value):
            body = json.dumps(value).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield gi(url=f"http://127.0.0.1:{server.server_port}", key="admin")
    finally:
        release.set()
        server.shutdown()
        server.server_close()


@pytest.fixture
def hanging_galaxy():
    with _fake_galaxy() as admin_gi:
        yield admin_gi


def _shed_test_config(admin_gi):
    config = BaseGalaxyConfig.__new__(BaseGalaxyConfig)
    config.galaxy_url = admin_gi.base_url
    config.master_api_key = "admin"
    config._kwds = {"shed_install_timeout": 1}
    return config


def test_hung_install_times_out(hanging_galaxy, capsys):
    with pytest.raises(Exception, match=FAILED_REPOSITORIES_MESSAGE):
        _install_shed_repos_from_tools_info([dict(REPOSITORY)], hanging_galaxy, False, install_timeout=1)
    err = capsys.readouterr().err
    assert "card/rgi" in err
    assert "84bd24ac33fd" in err
    assert hanging_galaxy.timeout is None


def test_hung_install_timeout_ignored_with_dependency_problems(hanging_galaxy):
    installed, _ = _install_shed_repos_from_tools_info([dict(REPOSITORY)], hanging_galaxy, True, install_timeout=1)
    assert installed == []


def test_hung_install_does_not_block_remaining_repositories(hanging_galaxy):
    other = dict(REPOSITORY, name="other", changeset_revision="0123456789ab")
    installed, _ = _install_shed_repos_from_tools_info(
        [dict(REPOSITORY), other], hanging_galaxy, True, install_timeout=1
    )
    assert [r["name"] for r in installed] == ["other"]


def test_shed_test_install_repo_times_out(hanging_galaxy):
    with pytest.raises(Exception, match="card/rgi at revision 84bd24ac33fd"):
        _shed_test_config(hanging_galaxy).install_repo(**REPOSITORY)


def test_shed_test_install_status_poll_times_out():
    with _fake_galaxy(hang_get=True) as admin_gi:
        with pytest.raises(requests.exceptions.Timeout):
            _shed_test_config(admin_gi).wait_for_all_installed()
