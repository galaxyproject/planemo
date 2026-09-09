"""Unit tests for how a database backend gets picked.

Nothing here probes ``PATH`` - the presence of a ``psql`` client says nothing about
whether a postgres server is reachable, so a backend is only used when it is named.
"""

import contextlib
import json
import os
from types import SimpleNamespace
from unittest import mock

import pytest

from planemo.database import (
    create_database_source,
    database_source_context,
)
from planemo.galaxy.config import DATABASE_LOCATION_TEMPLATE
from planemo.galaxy.profiles import (
    _create_profile_local,
    _profile_to_database_identifier,
    create_profile,
    delete_profile,
)
from .test_utils import TempDirectoryContext


def test_create_database_source_requires_a_named_backend():
    for kwds in ({}, {"database_type": "auto"}):
        with pytest.raises(Exception) as exc_info:
            create_database_source(**kwds)
        assert "--database_type" in str(exc_info.value), kwds


def test_create_database_source_dispatches_on_the_named_backend():
    source = create_database_source(database_type="postgres")
    assert type(source).__name__ == "LocalPostgresDatabaseSource"


def test_database_source_context_stops_sources_that_do_not_stay_running(tmp_path):
    source = mock.Mock(keep_running_after_database_commands=False)
    with mock.patch("planemo.database.factory.started_database_source", return_value=source):
        with database_source_context(
            database_type="postgres_singularity",
            postgres_storage_location=str(tmp_path / "postgres"),
        ) as yielded:
            assert yielded is source
            source.stop.assert_not_called()
    source.stop.assert_called_once_with()


def test_profile_defaults_to_its_own_sqlite_file():
    for kwds in ({}, {"database_type": None}, {"database_type": "auto"}, {"database_type": "sqlite"}):
        with TempDirectoryContext() as temp_directory_context:
            profile_directory = temp_directory_context.temp_directory
            with mock.patch("planemo.galaxy.profiles.database_source_context") as database_source_context:
                options = _create_profile_local(None, profile_directory, "profile1234", dict(kwds))
            database_source_context.assert_not_called()
            assert options["database_type"] == "sqlite", kwds
            database_location = os.path.join(profile_directory, "galaxy.sqlite")
            assert options["database_connection"] == DATABASE_LOCATION_TEMPLATE % database_location, kwds


def test_profile_creates_a_database_for_a_named_backend():
    database_source = mock.Mock()
    database_source.store_connection_in_profile = True
    database_source.sqlalchemy_url.return_value = "postgresql://galaxy@localhost/plnmoprof_profile1234"
    with TempDirectoryContext() as temp_directory_context:
        with mock.patch(
            "planemo.galaxy.profiles.database_source_context", return_value=contextlib.nullcontext(database_source)
        ):
            options = _create_profile_local(
                None, temp_directory_context.temp_directory, "profile1234", {"database_type": "postgres"}
            )
    identifier = _profile_to_database_identifier("profile1234")
    database_source.create_database.assert_called_once_with(identifier)
    assert options["database_type"] == "postgres"
    assert options["database_connection"] == "postgresql://galaxy@localhost/plnmoprof_profile1234"


def test_profile_creation_failure_is_not_swallowed_into_sqlite():
    database_source = mock.Mock()
    database_source.create_database.side_effect = RuntimeError("role does not exist")
    with TempDirectoryContext() as temp_directory_context:
        with mock.patch(
            "planemo.galaxy.profiles.database_source_context", return_value=contextlib.nullcontext(database_source)
        ):
            with pytest.raises(RuntimeError):
                _create_profile_local(
                    None, temp_directory_context.temp_directory, "profile1234", {"database_type": "postgres"}
                )


def test_singularity_profile_stores_managed_database_inputs_not_connection():
    database_source = mock.Mock()
    database_source.store_connection_in_profile = False
    database_source.profile_options.return_value = {
        "postgres_storage_location": "/profiles/profile1234/postgres",
        "singularity_cmd": "apptainer",
        "singularity_sudo": False,
    }
    database_source.sqlalchemy_url.return_value = (
        "postgresql://galaxy:mysecretpassword@/plnmoprof_profile1234?host=/profiles/profile1234/postgres/pgrun"
    )
    with TempDirectoryContext() as temp_directory_context:
        with mock.patch(
            "planemo.galaxy.profiles.database_source_context", return_value=contextlib.nullcontext(database_source)
        ):
            options = _create_profile_local(
                None,
                temp_directory_context.temp_directory,
                "profile1234",
                {
                    "database_type": "postgres_singularity",
                    "singularity_cmd": "apptainer",
                    "singularity_sudo": False,
                },
            )

    database_source.create_database.assert_called_once_with("plnmoprof_profile1234")
    assert "database_connection" not in options
    assert options == {
        "database_type": "postgres_singularity",
        "database_identifier": "plnmoprof_profile1234",
        "postgres_storage_location": "/profiles/profile1234/postgres",
        "singularity_cmd": "apptainer",
        "singularity_sudo": False,
        "engine": "galaxy",
    }


def test_failed_profile_creation_removes_partial_profile_directory(tmp_path):
    profiles_directory = tmp_path / "profiles"
    profiles_directory.mkdir()
    database_source = mock.Mock()
    database_source.create_database.side_effect = RuntimeError("database unavailable")
    context = contextlib.nullcontext(database_source)

    with mock.patch("planemo.galaxy.profiles.database_source_context", return_value=context):
        with pytest.raises(RuntimeError, match="database unavailable"):
            create_profile(
                SimpleNamespace(galaxy_profiles_directory=str(profiles_directory)),
                "profile1234",
                database_type="postgres",
            )

    assert not (profiles_directory / "profile1234").exists()


def test_singularity_profile_delete_drops_named_database_using_stored_location(tmp_path):
    profiles_directory = tmp_path / "profiles"
    profile_directory = profiles_directory / "profile1234"
    profile_directory.mkdir(parents=True)
    external_storage = tmp_path / "external-postgres"
    options = {
        "database_type": "postgres_singularity",
        "database_identifier": "plnmoprof_profile1234",
        "postgres_storage_location": str(external_storage),
        "singularity_cmd": "apptainer",
        "engine": "galaxy",
    }
    (profile_directory / "planemo_profile_options.json").write_text(json.dumps(options))
    database_source = mock.Mock()
    context = contextlib.nullcontext(database_source)

    with mock.patch("planemo.galaxy.profiles.database_source_context", return_value=context) as source_context:
        delete_profile(SimpleNamespace(galaxy_profiles_directory=str(profiles_directory)), "profile1234")

    source_context.assert_called_once_with(
        profile_directory=str(profile_directory),
        database_type="postgres_singularity",
        postgres_storage_location=str(external_storage),
        singularity_cmd="apptainer",
    )
    database_source.delete_database.assert_called_once_with("plnmoprof_profile1234")
    assert not profile_directory.exists()
