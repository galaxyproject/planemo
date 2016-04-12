"""The module describes a CLI framework extending ``click``."""
import os
import sys
import traceback

import click

from .io import error
from .config import read_global_config
from planemo import __version__

PYTHON_2_7_COMMANDS = ["cwl_run", "cwl_script"]
IS_PYTHON_2_7 = sys.version_info[0] == 2 and sys.version_info[1] >= 7


CONTEXT_SETTINGS = dict(auto_envvar_prefix='PLANEMO')
COMMAND_ALIASES = {
    "l": "lint",
    "t": "test",
    "s": "serve",
}


class Context(object):
    """Describe context of Planemo computation.

    Handles cross cutting concerns for Planemo such as verbose log
    tracking, the definition of the Planemo workspace (``~/.planemo``),
    and the global configuraton defined in ``~/.planemo.yml``.
    """

    def __init__(self):
        """Construct a Context object using execution environment."""
        self.home = os.getcwd()
        self._global_config = None
        # Will be set by planemo CLI driver
        self.verbose = False
        self.planemo_config = None
        self.planemo_directory = None
        self.option_source = {}

    def set_option_source(self, param_name, option_source):
        """Specify how an option was set."""
        assert param_name not in self.option_source
        self.option_source[param_name] = option_source

    def get_option_source(self, param_name):
        """Return OptionSource value indicating how the option was set."""
        assert param_name not in self.option_source
        return self.option_source[param_name]

    @property
    def global_config(self):
        """Read Planemo's global configuration.

        As defined most simply by ~/.planemo.yml.
        """
        if self._global_config is None:
            self._global_config = read_global_config(self.planemo_config)
        return self._global_config

    def log(self, msg, *args):
        """Log a message to stderr."""
        if args:
            msg %= args
        click.echo(msg, file=sys.stderr)

    def vlog(self, msg, *args, **kwds):
        """Log a message to stderr only if verbose is enabled."""
        if self.verbose:
            self.log(msg, *args)
            if kwds.get("exception", False):
                traceback.print_exc(file=sys.stderr)

    @property
    def workspace(self):
        """Create and return Planemo's workspace.

        By default this will be ``~/.planemo``.
        """
        if not self.planemo_directory:
            raise Exception("No planemo workspace defined.")
        workspace = self.planemo_directory
        return self._ensure_directory(workspace, "workspace")

    @property
    def galaxy_profiles_directory(self):
        """Create a return a directory for storing Galaxy profiles."""
        path = os.path.join(self.workspace, "profiles")
        return self._ensure_directory(path, "Galaxy profiles")

    def _ensure_directory(self, path, name):
        if not os.path.exists(path):
            os.makedirs(path)
        if not os.path.isdir(path):
            template = "Planemo %s directory [%s] unavailable."
            message = template % (name, path)
            raise Exception(message)
        return path


pass_context = click.make_pass_decorator(Context, ensure=True)
cmd_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          'commands'))


def list_cmds():
    """List planemo commands from commands folder."""
    rv = []
    for filename in os.listdir(cmd_folder):
        if filename.endswith('.py') and \
           filename.startswith('cmd_'):
            rv.append(filename[len("cmd_"):-len(".py")])
    rv.sort()
    if not IS_PYTHON_2_7:
        for command in PYTHON_2_7_COMMANDS:
            rv.remove(command)
    return rv


def _name_to_command(name):
    try:
        if sys.version_info[0] == 2:
            name = name.encode('ascii', 'replace')
        mod_name = 'planemo.commands.cmd_' + name
        mod = __import__(mod_name, None, None, ['cli'])
    except ImportError as e:
        error("Problem loading command %s, exception %s" % (name, e))
        return
    return mod.cli


class PlanemoCLI(click.MultiCommand):

    def list_commands(self, ctx):
        return list_cmds()

    def get_command(self, ctx, name):
        if name in COMMAND_ALIASES:
            name = COMMAND_ALIASES[name]
        return _name_to_command(name)


def command_function(f):
    """Extension point for processing kwds after click callbacks."""
    def outer(*args, **kwargs):
        # arg_spec = inspect.getargspec(f)
        return f(*args, **kwargs)
    outer.__doc__ = f.__doc__
    return pass_context(outer)


@click.command(cls=PlanemoCLI, context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
@click.option('-v', '--verbose', is_flag=True,
              help='Enables verbose mode.')
@click.option('--config',
              default="~/.planemo.yml",
              envvar="PLANEMO_GLOBAL_CONFIG_PATH",
              help="Planemo configuration YAML file.")
@click.option('--directory',
              default="~/.planemo",
              envvar="PLANEMO_GLOBAL_WORKSPACE",
              help="Workspace for planemo.")
@pass_context
def planemo(ctx, config, directory, verbose):
    """A command-line toolkit for building tools and workflows for Galaxy.

    Check out the full documentation for Planemo online
    http://planemo.readthedocs.org or open with ``planemo docs``.
    """
    ctx.verbose = verbose
    ctx.planemo_config = os.path.expanduser(config)
    ctx.planemo_directory = os.path.expanduser(directory)


__all__ = [
    "Context",
    "list_cmds",
    "command_function",
    "planemo",
]
