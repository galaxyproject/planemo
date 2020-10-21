"""Module describing the planemo ``list_invocations`` command."""
import json

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.galaxy.api import get_invocations
from planemo.galaxy.profiles import translate_alias
from planemo.io import error, info

try:
    from tabulate import tabulate
except ImportError:
    tabulate = None  # type: ignore


@click.command('list_invocations')
@click.argument(
    "workflow_id",
    type=click.STRING,
)
@options.profile_option(required=True)
@command_function
def cli(ctx, workflow_id, **kwds):
    """
    Get a list of invocations for a particular workflow ID or alias.
    """
    workflow_id = translate_alias(ctx, workflow_id, kwds.get('profile'))
    info("Looking for invocations for workflow {}...".format(workflow_id))
    workflow_id = profiles.translate_alias(ctx, workflow_id, kwds.get('profile'))
    profile = profiles.ensure_profile(ctx, kwds.get('profile'))

    invocations = get_invocations(url=profile['galaxy_url'], key=profile['galaxy_admin_key'] or profile['galaxy_user_key'], workflow_id=workflow_id)
    if tabulate:
        state_colors = {
            'ok': '\033[92m',           # green
            'running': '\033[93m',      # yellow
            'error': '\033[91m',        # red
            'paused': '\033[96m',       # cyan
            'deleted': '\033[95m',      # magenta
            'deleted_new': '\033[95m',  # magenta
            'new': '\033[96m',          # cyan
            'queued': '\033[93m',       # yellow
        }
        print(tabulate({
                "Invocation ID": invocations.keys(),
                "Jobs status": [', '.join(['{}{} jobs {}\033[0m'.format(state_colors[k], v, k) for k, v in inv['states'].items()]
                                          ) for inv in invocations.values()],
                "Invocation report URL": ['{}/workflows/invocations/report?id={}'.format(profile['galaxy_url'].strip('/'), inv_id
                                                                                         ) for inv_id in invocations],
                "History URL": ['{}/histories/view?id={}'.format(profile['galaxy_url'].strip('/'), invocations[inv_id]['history_id']
                                                                 ) for inv_id in invocations]
            }, headers="keys"))
    else:
        error("The tabulate package is not installed, invocations could not be listed correctly.")
        print(json.dumps(invocations, indent=4, sort_keys=True))
    info("{} invocations found.".format(len(invocations)))

    return
