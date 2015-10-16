"""
"""
import sys
import time

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.io import info, error
from planemo.reports import build_report


@click.command("shed_update")
@click.option(
    "--report_xunit",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output update as a faked XUnit report, useful when you want to "
         "automatically update repositories in bulk.",
    default=None,
)
@options.shed_publish_options()
@options.shed_upload_options()
@options.shed_skip_upload()
@options.shed_skip_metadata()
@pass_context
def cli(ctx, paths, **kwds):
    """Update repository in shed from a ``.shed.yml`` file.

    By default this command will update both repository metadata
    from ``.shed.yml`` and upload new contents from the repository
    directory.

    ::

        % planemo shed_update

    This will update the main tool shed with the repository defined
    by a ``.shed.yml`` file in the current working directory. Both
    the location of the ``.shed.yml`` and the tool shed to upload to
    can be easily configured. For instance, the following command can
    be used if ``.shed.yml`` if contained in ``path/to/repo`` and the
    desire is to update the test tool shed.

    ::

        % planemo shed_update --shed_target testtoolshed path/to/repo

    Another important option is ``--check_diff`` - this doesn't affect the
    updating of shed metadata but it will check for differences before
    uploading new contents to the tool shed. This may important because the
    tool shed will automatically populate certain attributes in tool shed
    artifact files (such as ``tool_dependencies.xml``) and this may
    cause unwanted installable revisions to be created when there are no
    important changes.

    The lower-level ``shed_upload`` command should be used instead if
    the repository doesn't define complete metadata in a ``.shed.yml``.
    """
    # In a little bit of cheating, we're defining this variable here to collect
    # a "report" on the shed_update command
    collected_data = {
        'results': {
            'total': 0,
            'errors': 0,
            'failures': 0,
            'skips': 0,
        },
        'suitename': 'update',
        'tests': [],
    }

    shed_context = shed.get_shed_context(ctx, **kwds)

    def update(realized_repository):
        collected_data['results']['total'] += 1
        upload_ret_code = 0
        upload_ok = True

        # Start the time
        time1 = time.time()
        if not kwds["skip_upload"]:
            upload_ret_code = shed.upload_repository(
                ctx, realized_repository, **kwds
            )
            upload_ok = not upload_ret_code
        time2 = time.time()

        # Now that we've uploaded (or skipped appropriately), collect results.
        if upload_ret_code == 2:
            collected_data['results']['failures'] += 1
            collected_data['tests'].append({
                'classname': realized_repository.name,
                'errorType': 'FailedUpdate',
                'errorMessage': 'Failed to update repository as it does not exist in target ToolShed',
                'time': (time2 - time1),
            })
            error("Failed to update repository it does not exist "
                  "in target ToolShed.")
            return upload_ret_code
        repo_id = realized_repository.find_repository_id(ctx, shed_context)
        metadata_ok = True
        if not kwds["skip_metadata"]:
            metadata_ok = realized_repository.update(ctx, shed_context, repo_id)
        if metadata_ok:
            info("Repository metadata updated.")
        else:
            error("Failed to update repository metadata.")
        if metadata_ok and upload_ok:
            collected_data['tests'].append({
                'classname': realized_repository.name,
                'time': (time2 - time1),
            })
            return 0
        elif upload_ok:
            collected_data['results']['skips'] += 1
            collected_data['tests'].append({
                'classname': realized_repository.name,
                'errorType': 'FailedMetadata',
                'errorMessage': 'Failed to update repository metadata',
                'time': (time2 - time1),
            })
            error("Repo updated but metadata was not.")
            return 1
        else:
            collected_data['results']['failures'] += 1
            collected_data['tests'].append({
                'classname': realized_repository.name,
                'errorType': 'FailedUpdate',
                'errorMessage': 'Failed to update repository',
                'time': (time2 - time1),
            })
            error("Failed to update a repository.")
            return 1

    exit_code = shed.for_each_repository(ctx, update, paths, **kwds)

    if kwds.get('report_xunit', False):
        with open(kwds['report_xunit'], 'w') as handle:
            handle.write(build_report.template_data(
                collected_data, template_name='xunit.tpl'))

    sys.exit(exit_code)
