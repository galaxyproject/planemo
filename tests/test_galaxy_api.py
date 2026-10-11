"""Unit tests for ``planemo.galaxy.api``."""

from planemo.galaxy.api import (
    DEFAULT_ADMIN_API_KEY,
    gi,
)


def test_gi_defaults_to_anonymous_access() -> None:
    galaxy_instance = gi(url="https://example.org")

    assert galaxy_instance.key is None


def test_gi_uses_explicit_api_key() -> None:
    galaxy_instance = gi(port=9090, key=DEFAULT_ADMIN_API_KEY)

    assert galaxy_instance.url == "http://localhost:9090/api"
    assert galaxy_instance.key == DEFAULT_ADMIN_API_KEY
