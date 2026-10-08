import pytest

from planemo.galaxy.api import gi


def test_gi_requires_a_url_or_port():
    with pytest.raises(ValueError, match="Either port or url must be provided"):
        gi()
