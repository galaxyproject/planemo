import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.reports import build_report
from planemo.io import captured_io_for_xunit


@click.command("check_urls")
@options.shed_read_options()
@click.option(
    "--report_xunit",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output diff as a faked XUnit report, useful when you want to "
         "automatically diff repositories and be warned when out-of-date.",
    default=None,
)
@pass_context
def cli(ctx, paths, **kwds):
    """Check URLs listed in package files for download.

    By default, this tool requests the first hundred or so bytes of each listed
    URL and validates that a 200OK was received.

        % # diff for this repository and the main Tool Shed
        % planemo shed_diff

    This command will return an exit code of:

    - 0 if all files are OK
    - 1 if any file has errors.
    """

    # In a little bit of cheating, we're defining this variable here to collect
    # a "report" on our shed_diff
    collected_data = {
        'results': {
            'total': 0,
            'errors': 0,
            'failures': 0,
            'skips': 0,
        },
        'suitename': 'check_urls',
        'tests': [],
    }

    def diff(realized_repository):
        captured_io = {}

        with captured_io_for_xunit(kwds, captured_io):
            results = shed.check_urls(ctx, realized_repository, **kwds)

        # Collect data about what happened
        collected_data['results']['total'] += 1
        xunit_case = {
            'name': 'check-urls',
            'classname': realized_repository.name,
            'time': captured_io["time"],
            'stdout': captured_io["stdout"],
            'stderr': captured_io["stderr"],
        }
        ret = 0
        if sum([r[0] for r in results]) > 0:
            inaccessible = [r for r in results if r[0] > 0]

            xunit_case.update({
                'errorType': 'InaccessibleUrls',
                'errorMessage': '%s/%s URLs inaccessible' % (len(inaccessible), len(results)),
                'errorContent': '\n'.join([r[1] for r in inaccessible]),
            })
            ret = 1

        # Append our xunit test case
        collected_data['tests'].append(xunit_case)
        return ret

    exit_code = shed.for_each_repository(ctx, diff, paths, **kwds)

    if kwds.get('report_xunit', False):
        with open(kwds['report_xunit'], 'w') as handle:
            handle.write(build_report.template_data(
                collected_data, template_name='xunit.tpl'))

    sys.exit(exit_code)
