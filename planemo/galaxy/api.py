"""A high-level interface to local Galaxy instances using bioblend."""
from six import StringIO

from planemo.bioblend import ensure_module
from planemo.bioblend import galaxy

DEFAULT_MASTER_API_KEY = "test_key"


def gi(port=None, url=None, key=None):
    """Return a bioblend ``GalaxyInstance`` for Galaxy on this port."""
    ensure_module()
    if key is None:
        key = DEFAULT_MASTER_API_KEY
    if port is None:
        url = url
    else:
        url = "http://localhost:%d" % int(port)

    return galaxy.GalaxyInstance(
        url=url,
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


def summarize_history(ctx, gi, history_id):
    """Summarize a history with print() based on similar code in Galaxy for populators.
    """
    if not ctx.verbose:
        return

    if history_id is None:
        raise ValueError("summarize_history passed empty history_id")
    try:
        history_contents = gi.histories.show_history(history_id, contents=True)
    except Exception:
        print("Failed to fetch history contents in summarize_history.")
        return

    for history_content in history_contents:
        history_content_id = history_content.get('id', None)
        print("| %d - %s (HID - NAME) " % (int(history_content['hid']), history_content['name']))
        if history_content['history_content_type'] == 'dataset_collection':
            history_contents_json = gi.histories.show_dataset_collection(history_id, history_content["id"])
            print("| Dataset Collection: %s" % history_contents_json)
            continue
        try:
            dataset_info = gi.histories.show_dataset(history_id, history_content_id)
            print("| Dataset State:")
            print(_format_for_summary(dataset_info.get("state"), "Dataset state is unknown."))
            print("| Dataset Blurb:")
            print(_format_for_summary(dataset_info.get("misc_blurb", ""), "Dataset blurb was empty."))
            print("| Dataset Info:")
            print(_format_for_summary(dataset_info.get("misc_info", ""), "Dataset info is empty."))
            print("| Peek:")
            print(_format_for_summary(dataset_info.get("peek", ""), "Peek unavilable."))
        except Exception:
            print("| *PLANEMO ERROR FETCHING DATASET DETAILS*")
        try:
            provenance_info = _dataset_provenance(gi, history_id, history_content_id)
            print("| Dataset Job Standard Output:")
            print(_format_for_summary(provenance_info.get("stdout", ""), "Standard output was empty."))
            print("| Dataset Job Standard Error:")
            print(_format_for_summary(provenance_info.get("stderr", ""), "Standard error was empty."))
        except Exception:
            print("| *PLANEMO ERROR FETCHING JOB DETAILS*")
        print("|")


def _format_for_summary(blob, empty_message, prefix="|  "):
    contents = "\n".join(["%s%s" % (prefix, line.strip()) for line in StringIO(blob).readlines() if line.rstrip("\n\r")])
    return contents or "%s*%s*" % (prefix, empty_message)


def _dataset_provenance(gi, history_id, id):
    provenance = gi.histories.show_dataset_provenance(history_id, id)
    return provenance


__all__ = (
    "DEFAULT_MASTER_API_KEY",
    "gi",
    "user_api_key",
)
