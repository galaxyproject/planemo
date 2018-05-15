
Publishing BioContainers_
----------------------------------------------------------------

Building unpublished BioContainers on the fly is great for testing but
for production use and to increase reproducibility such containers should
ideally be published as well.

BioContainers_ maintains a registry of package combinations to be published
using these long mulled hashes. This registry is represented as a Github repository
named `multi-package-containers <https://github.com/biocontainers/multi-package-containers>`__.
The Planemo command ``container_register`` will inspect a tool and open a
Github pull request to add the tool's combination
of packages to the registry. Once merged, this pull request will
result in the corresponding BioContainers image to be published (with the
correct mulled has as its name) - these can be subsequently be picked up by
Galaxy.

Various Github related settings need to be configured in order for Planemo
to be able to open pull requests on your behalf as part of the
``container_register`` command. To simplify all of this - the Planemo community
maintains a list of Github repositories containing Galaxy and/or CWL tools that
are scanned daily by Travis_. For each such repository, the Travis job will run
``container_register`` across the repository on all tools resulting in new registry
pull requests for all new combinations of tools. This list is maintained
in a script named ``monitor.sh`` in the `planemo-monitor
<https://github.com/galaxyproject/planemo-monitor/>`__ repository. The easiest way
to ensure new containers are built for your tools is simply to open open a pull
request to add your tool repositories to this list.


.. _BioContainers: http://biocontainers.pro/
.. _Travis: https://travis-ci.org
