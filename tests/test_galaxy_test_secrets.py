"""Tests for :mod:`planemo.galaxy.test.secrets` encoded-secret matching."""

import base64

import pytest

from planemo.galaxy.test.secrets import (
    _leak_problems,
    Secret,
)

SECRET = Secret("Secret 'TOKEN' of credentials 'svc'", "fake/token+4c1d2e")


@pytest.mark.parametrize("prefix", [b"", b"a", b"ab"])
def test_finds_base64_secret_at_any_alignment(prefix):
    content = base64.b64encode(prefix + SECRET.value.encode() + b" trailing")
    assert list(_leak_problems([SECRET], "job stdout", content)) == [
        "Secret 'TOKEN' of credentials 'svc' was found (base64-encoded) in job stdout."
    ]


def test_finds_url_encoded_secret():
    content = b"GET /?token=fake%2Ftoken%2B4c1d2e"
    assert list(_leak_problems([SECRET], "job stderr", content)) == [
        "Secret 'TOKEN' of credentials 'svc' was found (URL-encoded) in job stderr."
    ]


def test_reports_plain_match_once():
    content = SECRET.value.encode() * 2
    assert list(_leak_problems([SECRET], "output 'out'", content)) == [
        "Secret 'TOKEN' of credentials 'svc' was found in output 'out'."
    ]


def test_ignores_content_without_secret():
    assert list(_leak_problems([SECRET], "job stdout", b"token length 17")) == []
