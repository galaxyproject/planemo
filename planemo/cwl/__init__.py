"""Entry point for modules describing abstractions for dealing with CWL artifacts."""
from .run import run_cwltool
from .script import to_script


__all__ = [
    'run_cwltool',
    'to_script',
]
