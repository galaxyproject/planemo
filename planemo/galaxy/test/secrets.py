"""Fail Galaxy tool tests whose stdio or outputs contain one of the test's secrets.

Best effort: catches plain, base64 and URL-encoded copies, not otherwise transformed ones.
"""

import base64
import re
from itertools import chain
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
_MIN_SECRET_LENGTH = 8
_TO_STANDARD_BASE64 = bytes.maketrans(b"-_", b"+/")
# Shorter runs can't encode a _MIN_SECRET_LENGTH secret.
_BASE64_RUN = re.compile(rb"[A-Za-z0-9+/]{%d,}" % _MIN_SECRET_LENGTH)


class _Secret(NamedTuple):
    description: str
    value: str


def _collect_secrets(test_dict: Dict[str, Any], api_key: Optional[str]) -> List[_Secret]:
    secrets = []
    for credential in test_dict.get("credentials") or []:
        for secret in credential.get("secrets") or []:
            description = f"Secret '{secret['name']}' of credentials '{credential['name']}'"
            value = secret.get("value") or ""
            if len(value) < _MIN_SECRET_LENGTH:
                warn(f"{description} is shorter than {_MIN_SECRET_LENGTH} characters, not checking it for leaks.")
                continue
            secrets.append(_Secret(description, value))
    if api_key:
        secrets.append(_Secret("The test user's API key", api_key))
    return secrets


def _base64_decodings(content: bytes) -> List[bytes]:
    # Decode each base64-alphabet run at all four alignments; one of them lines up with any encoded secret.
    normalized = b"".join(content.split()).translate(_TO_STANDARD_BASE64)  # output is often wrapped
    decodings = []
    for run in _BASE64_RUN.findall(normalized):
        for offset in range(4):
            chunk = run[offset:]
            if len(chunk) % 4 == 1:  # a lone trailing character holds no whole byte
                chunk = chunk[:-1]
            decodings.append(base64.b64decode(chunk + b"=" * (-len(chunk) % 4)))
    return decodings


def _leak_problems(secrets: List[_Secret], location: str, content: bytes) -> Iterator[str]:
    decodings = _base64_decodings(content)
    for secret in secrets:
        raw = secret.value.encode()
        url_encoded = quote(secret.value, safe="").encode()
        if raw in content:
            encoding = ""
        elif any(raw in decoding for decoding in decodings):
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
        response = galaxy_interactor._get(f"datasets/{dataset_id}")
        response.raise_for_status()
        dataset = response.json()
        if dataset.get("misc_info"):
            yield f"the info of output '{name}'", dataset["misc_info"].encode()
        response = galaxy_interactor._get(f"datasets/{dataset_id}/display?raw=true")
        response.raise_for_status()
        yield f"output '{name}'", response.content


def check_for_secret_leaks(
    galaxy_interactor: GalaxyInteractorApi, test_dict: Dict[str, Any], job_data: Dict[str, Any]
) -> None:
    """Fail ``job_data`` if a secret of ``test_dict`` or the API key shows up in its job's stdio or outputs."""
    job = job_data.get("job")
    if not job:
        return
    problems: List[str] = []
    status = "failure"
    try:
        secrets = _collect_secrets(test_dict, galaxy_interactor.api_key)
        if not secrets:
            return
        for location, content in chain(_job_texts(job), _output_texts(galaxy_interactor, job["id"])):
            problems.extend(_leak_problems(secrets, location, content))
    except Exception as e:
        problems.append(f"Checking the test for leaked secrets failed: {e}")
        status = "error"
    if problems:
        job_data.setdefault("output_problems", []).extend(problems)
        if job_data.get("status") == "success":
            job_data["status"] = status
