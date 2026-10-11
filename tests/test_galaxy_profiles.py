"""Unit tests for ``planemo.galaxy.profiles``."""

import json
import os
from types import SimpleNamespace

import pytest

from planemo.galaxy import profiles

GALAXY_URL = "https://example.org"
ANONYMOUS_WARNING = "will access the external Galaxy instance anonymously"


def _create_external_profile(tmp_path, monkeypatch, **kwds):
    monkeypatch.setattr(profiles, "test_credentials_valid", lambda **_: True)
    ctx = SimpleNamespace(galaxy_profiles_directory=str(tmp_path))
    profiles.create_profile(ctx, "external", engine="external_galaxy", galaxy_url=GALAXY_URL, **kwds)
    with open(os.path.join(tmp_path, "external", profiles.PROFILE_OPTIONS_JSON_NAME)) as f:
        return json.load(f)


def test_external_profile_without_keys_warns_about_anonymous_access(tmp_path, monkeypatch, capsys) -> None:
    options = _create_external_profile(tmp_path, monkeypatch)

    assert ANONYMOUS_WARNING in capsys.readouterr().err
    assert options["galaxy_url"] == GALAXY_URL
    assert options["galaxy_user_key"] is None
    assert options["galaxy_admin_key"] is None


@pytest.mark.parametrize("key_option", ["galaxy_user_key", "galaxy_admin_key"])
def test_external_profile_with_key_does_not_warn(tmp_path, monkeypatch, capsys, key_option) -> None:
    options = _create_external_profile(tmp_path, monkeypatch, **{key_option: "secret"})

    assert ANONYMOUS_WARNING not in capsys.readouterr().err
    assert options[key_option] == "secret"
