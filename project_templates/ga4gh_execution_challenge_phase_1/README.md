# GA4GH Tool Execution Challenge - Phase 1 - Entry Planemo + Galaxy

This document describes how to use Planemo ([docs](http://planemo.readthedocs.io/en/latest/?badge=latest)
| [github](https://github.com/galaxyproject/planemo)) to complete Phase 1 of the
[GA4GH Tool Execution Challenge](https://www.synapse.org/#!Synapse:syn8080249) - running the example tool
through both a fork of Galaxy implementing initial support for CWL tool and [cwltool](https://github.com/common-workflow-language/cwltool) (the CWL reference implementation).

## Install Planemo

First install Planemo. Planemo is available for Linux and Mac OS X and can be installed using pip, Conda, or Homebrew and is also available as a virtual appliance for Docker, Vagrant, or VirtualBox. If you have virtualenv and Python installed the following is perhaps the easiest way to get going.

```
$ virtualenv planemo-venv; . planemo-venv/bin/activate
$ pip install -U pip
$ pip install planemo
```

Full documentation on installing Planemo is available [here](http://planemo.readthedocs.io/en/latest/installation.html).

## Setup Test Files

Planemo has a project and tool templates to quickly get going with CWL and Galaxy tool development. A project has been created for this phase that just populates a directory with these directions and the required input files.

```
$ planemo project_init --template ga4gh_execution_challenge_phase_1 phase_1
$ cd phase_1
$ ls 
README.md    job.json     md5sum.input
```

These files are:

* ``README.md`` is this document.
* ``md5sum.input`` is the document provided for the challenge to run through the challenge tool.
* ``job.json`` is a CWL job crafted to run ``md5sum.input`` with the supplied tool.

## Run the Test

Planemo can now be used to run the test tool.

```
$ time planemo run --docker --non_strict_cwl --output_directory . dockstore://quay.io/briandoconnor/dockstore-tool-md5sum job.json
** Lots of Galaxy output as a server is launched containing the referenced tool **
{u'output_file': {'path': u'/Users/john/workspace/planemo/project_templates/ga4gh_execution_challenge_phase_1/output_file', 'class': 'File'}}
planemo run --docker --non_strict_cwl --output_directory .  job.json  7.06s user 3.19s system 25% cpu 39.912 total
$ cat output_file
00579a00e3e7fa0674428ac7049423e2
```

The various Planemo flags here used are:

* ``--docker``: Unlike cwltool, Planemo doesn't attempt to use Docker by default but the tool used in the
  challenge requires a Docker container. This flags enables Docker - be sure Docker is running and available
  to the user of Planemo.
* ``--non_strict_cwl``: The test tool (at least one version of it) was non-strict CWL, this flag allows such
  tools to execute when using Planemo.
* ``--output_directory .``: This tells ``planemo run`` to write the tool outputs to the current directory.

As you can see, Planemo copied the output of that tool execution to the current directory as ``output_file`` (the
id of output in the CWL tool definition.

Planemo has the concept of engines - and it uses the Galaxy engine to run CWL tools by default but it can
also use ``cwltool`` (the reference implementation).

Below is an example of that.

```
$ time planemo run --engine cwltool --docker --non_strict_cwl dockstore://quay.io/briandoconnor/dockstore-tool-md5sum job.json
{u'output_file': {u'format': u'http://edamontology.org/data_3671', u'checksum': u'sha1$1250816c19c6f5524c5366b56c7a1eed6f3c3ab3', u'basename': u'md5sum.txt', u'location': u'file:///Users/john/workspace/planemo/project_templates/ga4gh_execution_challenge_phase_1/md5sum.txt', u'path': u'/Users/john/workspace/planemo/project_templates/ga4gh_execution_challenge_phase_1/md5sum.txt', u'class': u'File', u'size': 33}}
planemo run --engine cwltool --docker --non_strict_cwl  job.json  1.30s user 0.22s system 46% cpu 3.262 total
$ cat md5sum.txt
00579a00e3e7fa0674428ac7049423e2
```
