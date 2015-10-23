from __future__ import print_function
from __future__ import absolute_import

import contextlib
import os
import shutil
import sys
import tempfile
import time
from cStringIO import StringIO
from xml.sax.saxutils import escape

import click
from galaxy.tools.deps import commands
from galaxy.tools.deps.commands import which


def communicate(cmds, **kwds):
    info(cmds)
    p = commands.shell_process(cmds, **kwds)
    ret_val = p.communicate()
    if p.returncode != 0:
        template = "Problem executing commands {0} - ({1}, {2})"
        msg = template.format(cmds, ret_val[0], ret_val[1])
        raise RuntimeError(msg)
    return ret_val


def shell(cmds, **kwds):
    info(cmds)
    return commands.shell(cmds, **kwds)


def info(message, *args):
    if args:
        message = message % args
    _echo(click.style(message, bold=True, fg='green'))


def can_write_to_path(path, **kwds):
    if not kwds["force"] and os.path.exists(path):
        error("%s already exists, exiting." % path)
        return False
    return True


def error(message, *args):
    if args:
        message = message % args
    _echo(click.style(message, bold=True, fg='red'), err=True)


def warn(message, *args):
    if args:
        message = message % args
    _echo(click.style(message, fg='red'), err=True)


def _echo(message, err=False):
    if sys.version_info[0] == 2:
        click.echo(message, err=err)
    else:
        print(message)


def write_file(path, content):
    with open(path, "w") as f:
        f.write(content)


def untar_to(url, path, tar_args):
    if which("wget"):
        download_cmd = "wget -q --recursive -O - '%s'"
    else:
        download_cmd = "curl '%s'"
    download_cmd = download_cmd % url
    if tar_args:
        if not os.path.exists(path):
            os.makedirs(path)

        untar_cmd = "tar %s" % tar_args
        shell("%s | %s" % (download_cmd, untar_cmd))
    else:
        shell("%s > '%s'" % (download_cmd, path))


@contextlib.contextmanager
def temp_directory(prefix="planemo_tmp_"):
    temp_dir = tempfile.mkdtemp(prefix=prefix)
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir)


def kill_pid_file(pid_file):
    if not os.path.exists(pid_file):
        return

    pid = int(open(pid_file, "r").read())
    kill_posix(pid)


def kill_posix(pid):
    def _check_pid():
        try:
            os.kill(pid, 0)
            return True
        except OSError:
            return False

    if _check_pid():
        for sig in [15, 9]:
            try:
                os.kill(pid, sig)
            except OSError:
                return
            time.sleep(1)
            if not _check_pid():
                return


@contextlib.contextmanager
def captured_io_for_xunit(kwds, captured_io):
    captured_std = []
    with_xunit = kwds.get('report_xunit', False)
    if with_xunit:
        with Capturing() as captured_std:
            time1 = time.time()
            yield
            time2 = time.time()
        tee_captured_output(captured_std)
    else:
        time1 = time.time()
        yield
        time2 = time.time()

    if with_xunit:
        stdout = [escape(m['data']) for m in captured_std
                  if m['logger'] == 'stdout']
        stderr = [escape(m['data']) for m in captured_std
                  if m['logger'] == 'stderr']
        captured_io["stdout"] = stdout
        captured_io["stderr"] = stderr
        captured_io["time"] = (time2 - time1)
    else:
        captured_io["stdout"] = None
        captured_io["stderr"] = None
        captured_io["time"] = None


class Capturing(list):
    """Function context which captures stdout/stderr

    This keeps planemo's codebase clean without requiring planemo to hold onto
    messages, or pass user-facing messages back at all. This could probably be
    solved by swapping planemo entirely to a logger and reading from/writing
    to that, but this is easier.

    This swaps sys.std{out,err} with StringIOs and then makes that output
    available.
    """
    # http://stackoverflow.com/a/16571630

    def __enter__(self):
        self._stdout = sys.stdout
        self._stderr = sys.stderr
        sys.stdout = self._stringio_stdout = StringIO()
        sys.stderr = self._stringio_stderr = StringIO()
        return self

    def __exit__(self, *args):
        self.extend([{'logger': 'stdout', 'data': x} for x in
                     self._stringio_stdout.getvalue().splitlines()])
        self.extend([{'logger': 'stderr', 'data': x} for x in
                     self._stringio_stderr.getvalue().splitlines()])

        sys.stdout = self._stdout
        sys.stderr = self._stderr


def tee_captured_output(output):
    """For messages captured with Capturing, send them to their correct
    locations so as to not interfere with normal user experience.
    """
    for message in output:
        # Append '\n' due to `splitlines()` above
        if message['logger'] == 'stdout':
            sys.stdout.write(message['data'] + '\n')
        if message['logger'] == 'stderr':
            sys.stderr.write(message['data'] + '\n')
