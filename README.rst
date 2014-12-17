===============================
Planemo
===============================

.. image:: https://travis-ci.org/galaxyproject/planemo.png?branch=master
        :target: https://travis-ci.org/galaxyproject/planemo

.. image:: https://readthedocs.org/projects/pip/badge/?version=latest
		:target: https://planemo.readthedocs.org.

.. .. image:: https://badge.fury.io/py/planemo.png
..    :target: http://badge.fury.io/py/planemo
.. .. image:: https://pypip.in/d/planemo/badge.png
..        :target: https://pypi.python.org/pypi/planemo


Command-line utilities to assist in building tools for the Galaxy_ project.

* Free software: Academic Free License version 3.0
* Documentation: https://planemo.readthedocs.org.
* Code: https://github.com/galaxyproject/planemo

Quick Start
-----------

This quick start demonstrates using ``planemo`` commands to help
develop Galaxy tools. Planemo can quickly be installed via
Homebrew_ or as a more traditional Python project.

To install using Homebrew_ or linuxbrew_:

::

   brew tap galaxyproject/tap
   brew install planemo

For a more traditional Python installation simply setup a virtualenv
for ``planemo`` (this example creates a new one in ``.venv``) and then
install with ``pip``.

::

   % virtualenv .venv; . .venv/bin/activate
   % pip install git+https://github.com/galaxyproject/planemo.git

This quick start will assume you will have a directory with one or more
tool XML files. If none is available, one can be quickly create for
demonstrating ``planemo`` as follows ``mkdir mytools; cd mytools; planemo
project_init --template=demo``.

On OSX Mavericks, if gnutar is not installed, please install with homebrew using:

   % brew install gnu-tar;
   
   % cd /usr/bin;
   
   % sudo ln -s /usr/local/opt/gnu-tar/libexec/gnubin/tar gnutar  

Force the demo to initialize in the mytools directory:

   % cd mytools;
   
   % wget -q --recursive -O - 'https://codeload.github.com/galaxyproject/planemo/tar.gz/master' | gnutar -C . -zxvf - --strip-components=3 --wildcards --no-anchored 'demo/**'



Planemo can check tools containing XML for common problems and best 
practices using the ``lint`` `command <http://planemo.readthedocs.org/en/latest/commands.html#lint-command>`_
(also aliased as ``l``). ::

    % planemo lint
    ...

Like many ``planemo`` commands - by default this will search the 
current directory and use all tool files it finds. It can be explicitly
passed other tool files or a directory of tool files. ::

    % planemo l randomlines.xml

The ``lint`` command takes in a additional options related to 
reporting levels, exit code, etc.... These options are descirbed here
or (like all available commands) be accessed by passing it ``--help``.::

    % planemo l --help
    Usage: planemo lint [OPTIONS] TOOL_PATH
    ...

Once tools are syntically correct - it is time to test. The ``test`` 
`command <http://planemo.readthedocs.org/en/latest/commands.html#test-command>`_
can be used to test a tool or directory of tools.::

	% planemo test --galaxy_root=../galaxy-central randomlines.xml

If no ``--galaxy_root`` is defined, ``planemo`` will check for a default in
`~/.planemo.yml
<http://planemo.readthedocs.org/en/latest/configuration.html>`_) and finally
search the tool's parent directories for a Galaxy root directory (developing
tools under Galaxy ``tools`` directory is a common development workflow).
Planemo can also download and configure a disposable Galaxy instance just for
testing by passing it ``-install_galaxy`` instead of a Galaxy root.::

	% planemo t --install_galaxy

**Warning**: The features of Planemo that require a ``--galaxy_root`` will
only work with the latest ``galaxy-central`` ``default`` branch. Planemo will
not work when used with older versions of Galaxy - even the latest stable
``latest_2014.10.06``. Planemo will stablize with the next release of Galaxy
and serious attempts at backward compatibility going forward will be made at
that time.

Planemo will create a HTML an output report in the current directory named
``tool_test_output.html`` (override with ``--test_output``). `Here <http://galaxyproject.github.io/planemo/tool_test_viewer.html?test_data_url=https://gist.githubusercontent.com/jmchilton/9d4351c9545d34209904/raw/9ed285d3cf98e435fc4a743320363275949ad63c/index>`_ is an
example of such a report for Tophat.

Once tools have been linted and tested - the tools can be viewed in a
Galaxy interface using the ``serve`` (``s``) `command
<http://planemo.readthedocs.org/en/latest/commands.html#serve-command>`_.::

	% planemo serve

Like ``test``, ``serve`` requires a Galaxy root and one can be 
explicitly specified with ``--galaxy_root`` or installed dynamically
with ``--install_galaxy``.

Finally, when tools are ready to be published to GitHub_, it may be valuable
to setup contineous integration to test changes committed and new pull
requests. `Travis CI <http://travis-ci.org/>`_ is a service providing free
testing and deep integration with GitHub_.

The ``travis_init`` `command
<http://planemo.readthedocs.org/en/latest/commands.html#travis_init-command>`_
will bootstrap a project with files to ease  contineous inegration of tools
using a Planemo driven approach inspired by this great `blog post
<http://bit.ly/gxtravisci>`_ by `Peter Cock <https://github.com/peterjc>`_.

::

    % planemo travis_init .
    % # setup Ubuntu 12.04 w/tool dependencies
    % vim .travis/setup_custom_dependencies.bash
    % git add .travis.yml .travis
    % git commit -m "Add Travis CI testing infrastructure for tools."
    % git push # and register repository @ http://travis-ci.org/

Experimental Features
---------------------

While Planemo is very experimental itself - it can also be used to explore
some more experimental features related to Galaxy tooling - including support
for Docker and Brew.

-----------
Docker
-----------

Galaxy has `experimental support
<https://wiki.galaxyproject.org/Admin/Tools/Docker>`_ for running jobs in
Docker_ containers. Planemo contains tools to assist in development of Docker
images for Galaxy tools.

A shell can be launched to explore the Docker enviornment referenced by tools 
that are annotated with publically registered Docker images.::

    % $(planemo docker_shell bowtie2.xml)

For Docker containers still in development - a Dockerfile can be associated
with a tool by sticking it in the tool's directory. Planemo can then build
and tag a Docker image for this tool and launch a shell into it using the
following commands.::

    % planemo docker_build bowtie2.xml # asssumes Dockerfile in same dir
    % $(planemo docker_shell --from_tag bowtie2.xml)

For more details see the documentation for the `docker_build
<http://planemo.readthedocs.org/en/latest/commands.html#docker_build-command>`_
and `docker_shell
<http://planemo.readthedocs.org/en/latest/commands.html#docker_shell-command>`_
commands.

-----------
Brew
-----------

The Galaxy development team is exploring different options for integrating
Homebrew_ and linuxbrew_ with Galaxy. One angle is resolving the tool requirements
using ``brew``. An experimental approach for versioning of brew recipes will be
used. See full discussion on the homebrew-science issues page here -
https://github.com/Homebrew/homebrew-science/issues/1191. Information on the
implementation can be found https://github.com/jmchilton/platform-brew until a
more permanent project home is setup.

::

    % planemo brew_init # install linuxbrew (only need if not already installed)
    % planemo brew # install dependencies for all tools in directory.
    % planemo brew bowtie2.xml # install dependencies for one tool
    % which bowtie2
    bowtie2 not found
    % . <(planemo brew_env --shell bowtie2.xml) # shell w/brew deps resolved
    (bowtie2) % which bowtie2
    /home/john/.linuxbrew/Cellar/bowtie2/2.1.0/bin/bowtie2
    (bowtie2) % exit
    % . <(planemo brew_env bowtie2.xml) # or just source deps in cur env
    % which bowtie2
    /home/john/.linuxbrew/Cellar/bowtie2/2.1.0/bin/bowtie2

For more information see the documentation for the `brew
<http://planemo.readthedocs.org/en/latest/commands.html#brew-command>`_
and `brew_env
<http://planemo.readthedocs.org/en/latest/commands.html#brew_env-command>`_ commands.

.. _Galaxy: (http://galaxyproject.org/)
.. _GitHub: https://github.com/
.. _Docker: https://www.docker.com/
.. _Homebrew: http://brew.sh/
.. _linuxbrew: https://github.com/Homebrew/linuxbrew
