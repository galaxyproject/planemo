"""Entry point for modules describing abstractions for dealing with CWL artifacts."""
from .run import run_cwltool
from .script import to_script
from .toil import run_toil


__all__ = (
    'run_cwltool',
    'run_toil',
    'to_script',
)
