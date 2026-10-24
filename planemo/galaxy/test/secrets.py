"""Fail Galaxy tool tests whose stdio or outputs contain one of the test's secrets.

Best effort: catches plain, base64 and URL-encoded copies, not otherwise transformed ones.
"""

import base64
import binascii
from typing import (
    Any,
    Dict,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Tuple,
)
from urllib.parse import quote

from galaxy.tool_util.verify.interactor import GalaxyInteractorApi

from planemo.io import warn

_JOB_STREAMS = ("stdout", "stderr", "job_messages")
# Shorter values would match unrelated text and binary outputs.
MIN_SECRET_LENGTH = 8
_URLSAFE_BASE64 = bytes.maketrans(b"+/", b"-_")


class Secret(NamedTuple):
    """A secret value a test must not print, with how to name it in a failure."""

    description: str
    value: str


def collect_secrets(test_dict: Optional[Dict[str, Any]], api_key: Optional[str]) -> List[Secret]:
    """Collect the ``<secret>`` values of a test's ``<credentials>`` plus the test user's API key."""
    secrets = []
    for credential in (test_dict or {}).get("credentials") or []:
        for secret in credential.get("secrets") or []:
            description = f"Secret '{secret['name']}' of credentials '{credential['name']}'"
            value = secret.get("value") or ""
            if len(value) < MIN_SECRET_LENGTH:
                warn(f"{description} is shorter than {MIN_SECRET_LENGTH} characters, not checking it for leaks.")
                continue
            secrets.append(Secret(description, value))
    if api_key:
        secrets.append(Secret("The test user's API key", api_key))
    return secrets


def _base64_needles(raw: bytes) -> Iterator[bytes]:
    # One needle per alignment of the secret in the encoded stream, keeping only the
    # characters encoded purely from secret bytes.
    for offset in range(3):
        encoded = base64.b64encode(b"\0" * offset + raw)
        start = 4 if offset else 0
        end = (offset + len(raw)) // 3 * 4
        if end - start >= MIN_SECRET_LENGTH:
            yield encoded[start:end]


def _decodes_to_secret(raw: bytes, encoded: bytes) -> bool:
    # A needle drops up to two secret bytes at each end, so decode it with the
    # neighbouring groups to confirm the whole secret is there.
    for altchars in (None, b"-_"):
        for needle in _base64_needles(raw):
            if altchars:
                needle = needle.translate(_URLSAFE_BASE64)
            index = encoded.find(needle)
            while index != -1:
                for start in (max(index - 4, 0), index):
                    for end in (index + len(needle) + 4, index + len(needle)):
                        try:
                            if raw in base64.b64decode(encoded[start:end], altchars=altchars):
                                return True
                        except binascii.Error:
                            pass
                index = encoded.find(needle, index + 1)
    return False


def _leak_problems(secrets: List[Secret], location: str, content: bytes) -> Iterator[str]:
    unwrapped = b"".join(content.split())  # base64 output is often wrapped at 76 columns
    for secret in secrets:
        raw = secret.value.encode()
        url_encoded = quote(secret.value, safe="").encode()
        if raw in content:
            encoding = ""
        elif _decodes_to_secret(raw, unwrapped):
            encoding = " (base64-encoded)"
        elif url_encoded != raw and url_encoded in content:
            encoding = " (URL-encoded)"
        else:
            continue
        yield f"{secret.description} was found{encoding} in {location}."


def _job_texts(job: Dict[str, Any]) -> Iterator[Tuple[str, bytes]]:
    for stream in _JOB_STREAMS:
        if job.get(stream):
            yield f"job {stream}", str(job[stream]).encode()


def _output_texts(galaxy_interactor: GalaxyInteractorApi, job_id: str) -> Iterator[Tuple[str, bytes]]:
    for output in galaxy_interactor.job_outputs(job_id):
        if "dataset" not in output:
            continue  # a collection; its elements are listed as datasets
        dataset_id = output["dataset"]["id"]
        name = output["name"]
        dataset = galaxy_interactor._get(f"datasets/{dataset_id}").json()
        if dataset.get("misc_info"):
            yield f"the info of output '{name}'", dataset["misc_info"].encode()
        response = galaxy_interactor._get(f"datasets/{dataset_id}/display?raw=true")
        response.raise_for_status()
        yield f"output '{name}'", response.content


def check_for_secret_leaks(
    galaxy_interactor: GalaxyInteractorApi, job_data: Dict[str, Any], secrets: List[Secret]
) -> None:
    """Fail ``job_data`` if a secret shows up in its job's stdio or outputs."""
    job = job_data.get("job")
    if not secrets or not job:
        return
    problems: List[str] = []
    status = "failure"
    try:
        for location, content in _job_texts(job):
            problems.extend(_leak_problems(secrets, location, content))
        for location, content in _output_texts(galaxy_interactor, job["id"]):
            problems.extend(_leak_problems(secrets, location, content))
    except Exception as e:
        problems.append(f"Checking the test for leaked secrets failed: {e}")
        status = "error"
    if problems:
        job_data.setdefault("output_problems", []).extend(problems)
        if job_data.get("status") == "success":
            job_data["status"] = status
