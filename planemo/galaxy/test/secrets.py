"""Fail Galaxy tool tests whose stdio or outputs contain one of the test's secrets."""

import base64
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

JOB_STREAMS = ("stdout", "stderr", "job_messages")
MIN_BASE64_NEEDLE = 8


class Secret(NamedTuple):
    """A secret value a test must not print, with how to name it in a failure."""

    description: str
    value: str


def collect_secrets(test_dict: Optional[Dict[str, Any]], api_key: Optional[str]) -> List[Secret]:
    """Collect the ``<secret>`` values of a test's ``<credentials>`` plus the test user's API key."""
    secrets = []
    for credential in (test_dict or {}).get("credentials") or []:
        for secret in credential.get("secrets") or []:
            if secret.get("value"):
                secrets.append(
                    Secret(f"Secret '{secret['name']}' of credentials '{credential['name']}'", secret["value"])
                )
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
        if end - start >= MIN_BASE64_NEEDLE:
            yield encoded[start:end]


def _encodings(value: str) -> Iterator[Tuple[str, bytes]]:
    raw = value.encode()
    yield "", raw
    for needle in _base64_needles(raw):
        yield " (base64-encoded)", needle
    url_encoded = quote(value, safe="")
    if url_encoded != value:
        yield " (URL-encoded)", url_encoded.encode()


def _leak_problems(secrets: List[Secret], location: str, content: bytes) -> Iterator[str]:
    for secret in secrets:
        for encoding, needle in _encodings(secret.value):
            if needle in content:
                yield f"{secret.description} was found{encoding} in {location}."
                break


def _job_texts(job: Dict[str, Any]) -> Iterator[Tuple[str, bytes]]:
    for stream in JOB_STREAMS:
        if job.get(stream):
            yield f"job {stream}", str(job[stream]).encode()


def _output_texts(galaxy_interactor: GalaxyInteractorApi, job_id: str) -> Iterator[Tuple[str, bytes]]:
    for output in galaxy_interactor.job_outputs(job_id):
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
    for location, content in _job_texts(job):
        problems.extend(_leak_problems(secrets, location, content))
    for location, content in _output_texts(galaxy_interactor, job["id"]):
        problems.extend(_leak_problems(secrets, location, content))
    if problems:
        job_data.setdefault("output_problems", []).extend(problems)
        if job_data.get("status") == "success":
            job_data["status"] = "failure"
