"""Docker utilities for planemo.

Built on Galaxy abstractions in :module:`galaxy.tools.deps.dockerfiles` and
:module:`galaxy.tools.deps.docker_util`.
"""
from __future__ import absolute_import

from galaxy.tools.deps.dockerfiles import docker_host_args


__all__ = [
    'docker_host_args',
]
