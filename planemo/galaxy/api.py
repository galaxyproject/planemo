"""A high-level interface to local Galaxy instances using bioblend."""
from planemo.bioblend import ensure_module
from planemo.bioblend import galaxy

DEFAULT_MASTER_API_KEY = "test_key"


def gi(port, key=None):
    """Return a bioblend ``GalaxyInstance`` for Galaxy on this port."""
    ensure_module()
    if key is None:
        key = DEFAULT_MASTER_API_KEY
    return galaxy.GalaxyInstance(
        url="http://localhost:%d" % int(port),
        key=key
    )


def user_api_key(admin_gi):
    """Use an admin authenticated account to generate a user API key."""
    ensure_module()
    # TODO: thread-safe
    users = admin_gi.users
    all_users = users.get_users()

    user_id = None
    for user in all_users:
        if user["email"] == "planemo@galaxyproject.org":
            user_id = user["id"]

    if user_id is None:
        # TODO: Allow override with --user_api_key.
        galaxy_config = admin_gi.config.get_config()
        use_remote_user = bool(galaxy_config["use_remote_user"])
        if not use_remote_user:
            user_response = users.create_local_user(
                "planemo",
                "planemo@galaxyproject.org",
                "planemo",
            )
            user_id = user_response["id"]
        else:
            user_response = users.create_remote_user(
                "planemo@galaxyproject.org",
            )
            user_id = user_response["id"]
    return users.create_user_apikey(user_id)


__all__ = [
    "DEFAULT_MASTER_API_KEY",
    "gi",
    "user_api_key",
]
