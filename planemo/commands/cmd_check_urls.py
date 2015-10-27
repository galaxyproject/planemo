import sys
import re
from xml.sax.saxutils import escape

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
# @click.option(
    # "--verify",
    # type=bool,
    # help="If an sha256sum is available, download the entire file AND validate it.",
    # default=False,
# )
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

        for res in results:
            safe_test_name = res[1]
            if '/' in safe_test_name:
                safe_test_name = safe_test_name[safe_test_name.rindex('/'):]
            safe_test_name = re.sub(r'[^A-Za-z0-9_.-]', '', safe_test_name)

            xunit_case = {
                'name': safe_test_name,
                'classname': realized_repository.name,
                'time': captured_io["time"],
                'stdout': captured_io["stdout"],
                'stderr': captured_io["stderr"],
            }

            if res[0] > 0:
                xunit_case.update({
                    'errorType': 'InaccessibleUrl',
                    'errorMessage': '%s inaccessible' % escape(res[1]),
                })
                ret = 1
            collected_data['tests'].append(xunit_case)

        return ret

    exit_code = shed.for_each_repository(ctx, diff, paths, **kwds)

    if kwds.get('report_xunit', False):
        with open(kwds['report_xunit'], 'w') as handle:
            handle.write(build_report.template_data(
                collected_data, template_name='xunit.tpl'))

    sys.exit(exit_code)
