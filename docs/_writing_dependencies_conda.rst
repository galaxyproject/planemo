.. _dependencies_and_conda:

Dependencies and Conda
===========================================

----------------------------------------------------------------
Specifying and Using Tool Requirements
----------------------------------------------------------------

.. note:: This document discusses using Conda to satisfy tool dependencies from a tool developer
    perspective. An in depth discussion of using Conda to satisfy dependencies from an
    admistrator's perspective can be found `here <https://docs.galaxyproject.org/en/latest/admin/conda_faq.html>`__.
    That document also serves as good background for this discussion.

.. include:: _writing_conda_init.rst

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

.. figure:: images/dependency_resolution.png
   :alt: Diagram describing mapping tool requirements to dependencies.

   Galaxy's dependency resolution maps tool requirement tags to concrete
   applications and libraries setup by the Galaxy deployer (or Planemo). As
   the above diagram indicates the same requirements may be used by multiple
   Galaxy tools and a single Galaxy tool may depend on multiple requirements.
   The document describes working with Conda dependencies from a developer
   perspective but other dependency resolution techniques are covered in
   the `Galaxy docs <https://docs.galaxyproject.org/en/latest/admin/dependency_resolvers.html>`__.

.. include:: _writing_conda_overview.rst

We can check if the requirements on a tool are available in best practice
Conda channels using an extended form of the ``planemo lint`` command. Passing
``--conda_requirements`` flag will ensure all listed requirements are found.

::

    $ planemo lint --conda_requirements seqtk_seq.xml
    Linting tool /Users/john/workspace/planemo/docs/writing/seqtk_seq.xml
      ...
    Applying linter requirements_in_conda... CHECK
    .. INFO: Requirement [seqtk@1.2] matches target in best practice Conda channel [bioconda].


.. note:: You can download the final version of the seqtk seq wrapper from the Planemo
    tutorial using the command::

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

The above install worked properly, but ``seqtk`` is not on your ``PATH`` because this merely
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
    2017-02-22 10:13:29,317 INFO  [galaxy.jobs.command_factory] Built script [/tmp/tmpLvKwta/job_working_directory/000/2/tool_script.sh] for tool command [[ "$CONDA_DEFAULT_ENV" = "/Users/john/miniconda2/envs/__seqtk@1.2" ] || . /Users/john/miniconda2/bin/activate '/Users/john/miniconda2/envs/__seqtk@1.2' >conda_activate.log 2>&1 ; seqtk seq -A '/tmp/tmpLvKwta/files/000/dataset_1.dat' > '/tmp/tmpLvKwta/files/000/dataset_2.dat']
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

.. include:: _writing_conda_search.rst

----------------------------------------------------------------
Exercise - Leveraging Bioconda
----------------------------------------------------------------

Use the ``project_init`` command to download this exercise.

::

    $ planemo project_init --template conda_exercises conda_exercises
    $ cd conda_exercises/exercise_1
    $ ls
    pear.xml              test-data

This project template contains a few exercises. The first uses an adapted
version of an IUC tool for `PEAR - Paired-End reAd mergeR
<http://sco.h-its.org/exelixis/web/software/pear/>`__. This tool however has
no ``requirement`` tags and so will not work properly.

1. Run ``planemo test pear.xml`` to verify the tool does not function
   without dependencies defined.
2. Use ``--conda_requirements`` flag with ``planemo lint`` to verify it does
   indeed lack requirements.
3. Use ``planemo conda_search`` or the Anaconda_ website to search for the
   correct package and version in a best practice channel.
4. Update ``pear.xml`` with the correct ``requirement`` tags.
5. Re-run the ``lint`` command from above to verify the tool now has the
   correct dependency definition.
6. Re-run the ``test`` command from above to verify the tool test now
   works properly.

.. include:: _writing_conda_new.rst

::

    $ planemo project_init --template conda_exercises conda_exercises
    $ cd conda_exercises/exercise_2
    $ ls
    fleeqtk_seq.xml              test-data

.. include:: _writing_conda_fleeqtk.rst

1. Clone and branch Bioconda_.
2. Build a recipe for fleeqtk version 1.3. You may wish to start from scratch
   (``conda skeleton`` is not available for C programs like fleeqtk), or copy
   the recipe of seqtk and modify it for fleeqtk.
3. Use ``conda build`` or Bioconda tooling to build the recipe.
4. Run ``planemo test --conda_use_local fleeqtk_seq.xml`` to verify the resulting package works as expected.

.. note: The planemo flag ``--conda_use_local`` causes Planemo to use locally built
     packages during dependency resolution and related commands.

.. include:: _writing_conda_recipe_complete.rst

.. _fleeqtk: https://github.com/jmchilton/fleeqtk
.. _Bioconda: https://github.com/bioconda/bioconda-recipes
.. _Conda: https://conda.io/docs/
.. _Anaconda: https://anaconda.org/
