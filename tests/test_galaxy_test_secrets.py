"""Tests for :mod:`planemo.galaxy.test.secrets` encoded-secret matching."""

import base64

import pytest

from planemo.galaxy.test.secrets import (
    _leak_problems,
    check_for_secret_leaks,
    collect_secrets,
    Secret,
)

SECRET = Secret("Secret 'TOKEN' of credentials 'svc'", "fake/token+4c1d2e")


@pytest.mark.parametrize("prefix", [b"", b"a", b"ab"])
def test_finds_base64_secret_at_any_alignment(prefix):
    content = base64.b64encode(prefix + SECRET.value.encode() + b" trailing")
    assert list(_leak_problems([SECRET], "job stdout", content)) == [
        "Secret 'TOKEN' of credentials 'svc' was found (base64-encoded) in job stdout."
    ]


def test_finds_wrapped_and_urlsafe_base64_secret():
    encoded = base64.urlsafe_b64encode(b"x" * 60 + SECRET.value.encode())
    wrapped = encoded[:76] + b"\n" + encoded[76:]
    assert list(_leak_problems([SECRET], "job stdout", wrapped)) == [
        "Secret 'TOKEN' of credentials 'svc' was found (base64-encoded) in job stdout."
    ]


@pytest.mark.parametrize("prefix", [b"", b"a", b"ab"])
def test_finds_unpadded_base64_secret_ending_the_stream(prefix):
    content = b"jwt=" + base64.urlsafe_b64encode(prefix + SECRET.value.encode()).rstrip(b"=")
    assert list(_leak_problems([SECRET], "job stdout", content)) == [
        "Secret 'TOKEN' of credentials 'svc' was found (base64-encoded) in job stdout."
    ]


def test_ignores_base64_of_partial_secret():
    secret = Secret("Secret 'TOKEN' of credentials 'svc'", "abcdefgh")
    assert list(_leak_problems([secret], "job stdout", base64.b64encode(b"abcdefZZ"))) == []


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


def test_skips_short_secrets():
    test_dict = {"credentials": [{"name": "svc", "secrets": [{"name": "TOKEN", "value": "secret"}]}]}
    assert collect_secrets(test_dict, None) == []


class _FailingInteractor:
    def job_outputs(self, job_id):
        raise Exception("Galaxy is down")


def test_scan_errors_fail_closed():
    job_data = {"status": "success", "job": {"id": "1", "stdout": ""}}
    check_for_secret_leaks(_FailingInteractor(), job_data, [SECRET])
    assert job_data["status"] == "error"
    assert job_data["output_problems"] == ["Checking the test for leaked secrets failed: Galaxy is down"]
