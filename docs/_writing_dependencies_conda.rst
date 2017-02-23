Dependencies and Conda
===========================================

----------------------------------------------------------------
Specifying and Using Tool Requirements
----------------------------------------------------------------

.. note:: This document discusses using Conda to satisfy tool dependencies from a tool developer
    perspective. An in depth discussion of using Conda to satisfy dependencies from an
    admistrator's perspective can be found `here <https://docs.galaxyproject.org/en/latest/admin/conda_faq.html>`__.
    That document also serves as good background for this discussion.

.. note:: Planemo requires a Conda installation to target with its various Conda
    related commands. A properly configured Conda installation can be initialized
    with the ``conda_init`` command. This should only need to be executed once
    per development machine.

    ::

        $ planemo conda_init
        wget -q --recursive -O '/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/conda_installLuGDHE.sh' 'https://repo.continuum.io/miniconda/Miniconda3-4.2.12-MacOSX-x86_64.sh' && bash '/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/conda_installLuGDHE.sh' -b -p '/Users/john/miniconda2' && /Users/john/miniconda2/bin/conda install -y -q conda=4.2.13
        PREFIX=/Users/john/miniconda2
        installing: python-3.5.2-0 ...
        installing: conda-env-2.6.0-0 ...
        installing: openssl-1.0.2j-0 ...
        installing: pycosat-0.6.1-py35_1 ...
        installing: readline-6.2-2 ...
        installing: requests-2.11.1-py35_0 ...
        installing: ruamel_yaml-0.11.14-py35_0 ...
        installing: sqlite-3.13.0-0 ...
        installing: tk-8.5.18-0 ...
        installing: xz-5.2.2-0 ...
        installing: yaml-0.1.6-0 ...
        installing: zlib-1.2.8-3 ...
        installing: conda-4.2.12-py35_0 ...
        installing: pycrypto-2.6.1-py35_4 ...
        installing: pip-8.1.2-py35_0 ...
        installing: wheel-0.29.0-py35_0 ...
        installing: setuptools-27.2.0-py35_0 ...
        Python 3.5.2 :: Continuum Analytics, Inc.
        creating default environment...
        installation finished.
        Fetching package metadata .......
        Solving package specifications: ..........

        Package plan for installation in environment /Users/john/miniconda2:

        The following packages will be downloaded:

            package                    |            build
            ---------------------------|-----------------
            conda-4.2.13               |           py35_0         389 KB

        The following packages will be UPDATED:

            conda: 4.2.12-py35_0 --> 4.2.13-py35_0

        Conda installation succeeded - Conda is available at '/Users/john/miniconda2/bin/conda'

While Galaxy can be configured to resolve dependencies various ways, Planemo
is configured with opinionated defaults geared at making building tools that
target Conda_ as easy as possible.

During the introductory tool development tutorial, we called ``planemo tool_init``
with the argument ``--requirement seqtk@1.2`` and the resulting tool contained
the XML::

    <requirements>
        <requirement type="package" version="1.2">seqtk</requirement>
    </requirements>

As configured by Planemo, when Galaxy encounters these ``requirement`` tags it
will attempt to install Conda, check for referenced packages (such as
``seqtk``), and install them as needed for tool testing.

We can check if the requirements on a tool are available in best practice
Conda channels using an extended form of the ``planemo lint`` command. Passing
``--conda_requirements`` flag will ensure all listed requirements are found.

::

    $ planemo lint --conda_requirements seqtk_seq.xml
    Linting tool /Users/john/workspace/planemo/docs/writing/seqtk_seq_v6.xml
      ...
    Applying linter requirements_in_conda... CHECK
    .. INFO: Requirement [seqtk@1.2] matches target in best practice Conda channel [bioconda].


.. note:: You can download the final version of the seqtk from the Planemo tutorial using
    the command::

        $ planemo project_init --template=seqtk_complete seqtk_example
        $ cd seqtk_example

We can verify these tool requirements install with the ``conda_install`` command. With
its default parameters ``conda_install`` processes tools and creates isolated environments
for their declared requirements.

::

    $ planemo conda_install seqtk_seq.xml
    Install conda target CondaTarget[seqtk,version=1.2]
    /home/john/miniconda2/bin/conda create -y --name __seqtk@1.2 seqtk=1.2
    Fetching package metadata ...............
    Solving package specifications: ..........

    Package plan for installation in environment /home/john/miniconda2/envs/__seqtk@1.2:

    The following packages will be downloaded:

        package                    |            build
        ---------------------------|-----------------
        seqtk-1.2                  |                0          29 KB  bioconda

    The following NEW packages will be INSTALLED:

        seqtk: 1.2-0   bioconda
        zlib:  1.2.8-3

    Fetching packages ...
    seqtk-1.2-0.ta 100% |#############################################################| Time: 0:00:00 444.71 kB/s
    Extracting packages ...
    [      COMPLETE      ]|################################################################################| 100%
    Linking packages ...
    [      COMPLETE      ]|################################################################################| 100%
    #
    # To activate this environment, use:
    # > source activate __seqtk@1.2
    #
    # To deactivate this environment, use:
    # > source deactivate __seqtk@1.2
    #
    $ which seqtk
    seqtk not found
    $

The above install worked properly, but seqtk is not on your ``PATH`` because this merely
created an environment within the Conda directory for the seqtk installation. Planemo
will configure Galaxy to exploit this installation. If you wish to interactively explore
the resulting enviornment to explore the installed tool or produce test data the output
of the ``conda_env`` command can be sourced.

::

    $ . <(planemo conda_env seqtk_seq.xml)
    Deactivate environment with conda_env_deactivate
    (seqtk_seq) $ which seqtk
    /home/planemo/miniconda2/envs/jobdepsiJClEUfecc6d406196737781ff4456ec60975c137e04884e4f4b05dc68192f7cec4656/bin/seqtk
    (seqtk_seq) $ seqtk seq

    Usage:   seqtk seq [options] <in.fq>|<in.fa>

    Options: -q INT    mask bases with quality lower than INT [0]
             -X INT    mask bases with quality higher than INT [255]
             -n CHAR   masked bases converted to CHAR; 0 for lowercase [0]
             -l INT    number of residues per line; 0 for 2^32-1 [0]
             -Q INT    quality shift: ASCII-INT gives base quality [33]
             -s INT    random seed (effective with -f) [11]
             -f FLOAT  sample FLOAT fraction of sequences [1]
             -M FILE   mask regions in BED or name list FILE [null]
             -L INT    drop sequences with length shorter than INT [0]
             -c        mask complement region (effective with -M)
             -r        reverse complement
             -A        force FASTA output (discard quality)
             -C        drop comments at the header lines
             -N        drop sequences containing ambiguous bases
             -1        output the 2n-1 reads only
             -2        output the 2n reads only
             -V        shift quality by '(-Q) - 33'
             -U        convert all bases to uppercases
             -S        strip of white spaces in sequences
    (seqtk_seq) $ conda_env_deactivate
    $

As shown above the ``conda_env_deactivate`` will be created in this environment and can
be used to restore your initial shell configuration.

Confident the underlying application works, we can now use ``planemo test`` or
``planemo serve`` and it will reuse this environment and find our dependency (in this
case ``seqtk`` as needed).

Here is a portion of the output from the testing command ``planemo test seqtk_seq.xml``
demonstrating using this tool.

::

    $ planemo test seqtk_seq.xml
    ...
    2017-02-22 10:13:28,902 INFO  [galaxy.tools.actions] Handled output named output1 for tool seqtk_seq (20.136 ms)
    2017-02-22 10:13:28,914 INFO  [galaxy.tools.actions] Added output datasets to history (12.782 ms)
    2017-02-22 10:13:28,935 INFO  [galaxy.tools.actions] Verified access to datasets for Job[unflushed,tool_id=seqtk_seq] (10.954 ms)
    2017-02-22 10:13:28,936 INFO  [galaxy.tools.actions] Setup for job Job[unflushed,tool_id=seqtk_seq] complete, ready to flush (21.053 ms)
    2017-02-22 10:13:28,962 INFO  [galaxy.tools.actions] Flushed transaction for job Job[id=2,tool_id=seqtk_seq] (26.510 ms)
    2017-02-22 10:13:29,064 INFO  [galaxy.jobs.handler] (2) Job dispatched
    2017-02-22 10:13:29,281 DEBUG [galaxy.tools.deps] Using dependency seqtk version 1.2 of type conda
    2017-02-22 10:13:29,282 DEBUG [galaxy.tools.deps] Using dependency seqtk version 1.2 of type conda
    2017-02-22 10:13:29,317 INFO  [galaxy.jobs.command_factory] Built script [/tmp/tmpLvKwta/job_working_directory/000/2/tool_script.sh] for tool command [[ "$CONDA_DEFAULT_ENV" = "/Users/john/miniconda2/envs/__seqtk@1.2" ] || . /Users/john/miniconda2/bin/activate '/Users/john/miniconda2/envs/__seqtk@1.2' >conda_activate.log 2>&1 ; seqtk seq -a '/tmp/tmpLvKwta/files/000/dataset_1.dat' > '/tmp/tmpLvKwta/files/000/dataset_2.dat']
    2017-02-22 10:13:29,516 DEBUG [galaxy.tools.deps] Using dependency samtools version None of type conda
    2017-02-22 10:13:29,516 DEBUG [galaxy.tools.deps] Using dependency samtools version None of type conda
    ok
    
    ----------------------------------------------------------------------
    XML: /private/tmp/tmpLvKwta/xunit.xml
    ----------------------------------------------------------------------
    Ran 1 test in 15.936s
    
    OK
    2017-02-22 10:13:37,014 INFO  [test_driver] Shutting down
    2017-02-22 10:13:37,014 INFO  [test_driver] Shutting down embedded galaxy web server
    2017-02-22 10:13:37,016 INFO  [test_driver] Embedded web server galaxy stopped
    2017-02-22 10:13:37,017 INFO  [test_driver] Stopping application galaxy
    ....
    2017-02-22 10:13:37,018 INFO  [galaxy.jobs.handler] sending stop signal to worker thread
    2017-02-22 10:13:37,018 INFO  [galaxy.jobs.handler] job handler stop queue stopped
    Testing complete. HTML report is in "/Users/john/workspace/planemo/project_templates/seqtk_complete/tool_test_output.html".
    All 1 test(s) executed passed.
    seqtk_seq[0]: passed

In this case the tests passed and the line containing ``[galaxy.tools.deps] Using dependency seqtk version 1.2 of type conda``
indicates Galaxy dependency resolution was successful and it found the environment we previously installed with ``conda_install``.

----------------------------------------------------------------
Finding Existing Conda Packages
----------------------------------------------------------------

How did we know what software name and software version to use? We found the existing
packages available for Conda and referenced them. To do this yourself, you can simply
use the planemo command ``conda_search``. If we do a search for ``seqt`` it will show
all the software and all the versions available matching that search term - including
``seqtk``.

::

    $ planemo conda_search seqt
    Fetching package metadata ...............
    seqtk                        r75                           0  bioconda
                                 r82                           0  bioconda
                                 r93                           0  bioconda
                                 1.2                           0  bioconda

.. note:: The Planemo command ``conda_search`` is a light wrapper around the underlying
   ``conda search`` command but configured to use the same channels and other options as
   Planemo and Galaxy. The following Conda command would also work to search::

       $ $HOME/miniconda3/bin/conda -c iuc -c bioconda -c conda-forge seqt

Alternatively the Anaconda_ website can be used to search for packages. Typing ``seqtk``
into the search form on that page and clicking the top result will bring on to `this page
https://anaconda.org/bioconda/seqtk`__ with information about the Bioconda package.

When using the website to search though, you need to aware of what channel you are using. By
default, Planemo and Galaxy will search a few different Conda channels. While it is possible
to configure a local Planemo or Galaxy to target different channels - the current best practice
is to add tools to the existing channels.

The existing channels include:

* Bioconda (`github <https://github.com/bioconda/bioconda-recipes>`__ | `conda <https://anaconda.org/bioconda>`__) - best practice channel for various bioinformatics packages.
* Conda-Forge (`github <https://github.com/conda-forge/staged-recipes>`__ | `conda <https://anaconda.org/conda-forge>`__) - best practice channel for general purpose and widely useful computing packages and libraries.
* iuc (`github <https://github.com/galaxyproject/conda-iuc>`__ | `conda <https://anaconda.org/iuc>`__) - best practice channel for other more Galaxy specific packages.

----------------------------------------------------------------
Exercise - Leveraging Bioconda
----------------------------------------------------------------

Use the ``project_init`` command to download this exercise.

::

    $ planemo project_init --template conda_exercise conda_exercise
    $ cd conda_exercise
    $ ls 
    pear.xml              test-data

This will download a tool for `PEAR - Paired-End reAd mergeR
<http://sco.h-its.org/exelixis/web/software/pear/>`__. This tool however has
no ``requirement`` tags and so will not work properly.

1. Run ``planemo test pear.xml`` to verify the tool does not function
   without dependencies defined.
1. Use ``--conda_requirements`` flag with ``planemo lint`` to verify it does
   indeed lack requirements.
1. Use ``planemo conda_search`` or the Anaconda_ website to search for the
   correct package and version in a best practice channel.
1. Update ``pear.xml`` with the correct ``requirement`` tags.
1. Re-run the ``lint`` command from above to verify the tool now has the
   correct dependency definition.
1. Re-run the ``test`` command from above to verify the tool test now
   works properly.

----------------------------------------------------------------
Building New Conda Packages
----------------------------------------------------------------

Frequently packages your tool will require are not found in Bioconda_
or conda-forge yet. In these cases, it is likely best to contribute
your package to one of these projects. Unless the tool is exceedingly
general Bioconda_ is usually the correct starting point.

.. note:: Many things that are not strictly or even remotely "bio" have
    been accepted into Bioconda_ - including tools for image analysis,
    natural language processing, and cheminformatics.

At this time, the most relevant source for information on building Conda packages for Galaxy
is probably the Bioconda_ documentation - in particular check out the `contributing documentation
<https://bioconda.github.io/contributing.html>`__.

----------------------------------------------------------------
Exercise - Build a Recipe
----------------------------------------------------------------

1. Build a recipe for Bioconda.
2. Open a pull request to contribute it.


.. _Bioconda: https://github.com/bioconda/bioconda-recipes
.. _Conda: https://conda.io/docs/
.. _Anaconda: https://anaconda.org/
