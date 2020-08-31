"""Planemo layer for ensuring bioblend available."""
try:
    from bioblend import toolshed
    from bioblend import galaxy
except ImportError:
    toolshed = None
    galaxy = None

BIOBLEND_UNAVAILABLE = ("This functionality requires the bioblend library "
                        " which is unavailable, please install `pip install "
                        "bioblend`")


def ensure_module():
    """Throw an exception if bioblend is not available to Planemo."""
    if toolshed is None:
        raise Exception(BIOBLEND_UNAVAILABLE)
