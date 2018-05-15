.. _dependencies_and_conda_cwl:

Dependencies and Conda
===========================================

----------------------------------------------------------------
Specifying and Using `Software Requirements`_
----------------------------------------------------------------

.. include:: _writing_conda_init.rst

.. note:: Why not just use containers?

    Containers are great, use containers (be it Docker, Singularity, etc.) whenever possible to
    increase reproducibility and portability - but building ad hoc containers to support CWL
    tools has some limitations that this document describes a process for addressing.

    There are technical reasons to describe `Software Requirements`_ in addition or in lieu
    of just using ad hoc containers - it will allow your tool to be used in environments without
    container runtimes available and the containers built from Conda software requirements are very
    likely to be "best practice" (e.g. smaller than ad hoc containers). Perhaps the most important
    reasons are less technical however such as reducing the opaqueness of traditional Docker
    containers. 

    Read more about this in our preprint `Practical computational reproducibility in the life sciences
    <https://www.biorxiv.org/content/early/2017/10/10/200683>`__

The Common Workflow Language specification loosely describes
`Software Requirements`_ - a way to map CWL hints to packages, environment
modules, or any other mechanism to describe dependencies for running a tool
outside of a container. The large and active Galaxy tool development community
has built a library and set of best practices for describing dependencies
for Galaxy that should work just as well for CWL. The library has been integrated
with cwltool_ and Toil_ to enable CWL tool authors and users to leverage the
power and flexibility of the Galaxy dependency management and best practices.

While `Software Requirements`_ can be configured to resolve dependencies various ways,
Planemo is configured with opinionated defaults geared at making building CWL tools that
target Conda_ as easy as possible and build tools with requirements compatible with
cwltool_ and Toil_ when running outside containers.

During the tool development `introductory tutorial`_, we called ``planemo tool_init``
with the argument ``--requirement seqtk@1.2`` and the resulting tool contained such
a ``SoftwareRequirement`` in the form the following the YAML fragment::

    SoftwareRequirement:
      packages:
      - package: seqtk
        version:
        - "1.2"

Planemo (and cwltool_ and Toil_) can interpret these ``SoftwareRequirement`` in varoius ways
including as Conda packages and install Conda packages referenced this way (such as ``seqtk``),
and install them as needed for tool testing.

We can check if the requirements on a tool are available in best practice
Conda channels using an extended form of the ``planemo lint`` command (``planemo lint`` was
introduced in the `introductory tutorial`_). Passing ``--conda_requirements`` flag will ensure all
listed requirements are found.

::

    $ planemo lint --conda_requirements seqtk_seq.cwl
    Linting tool /Users/john/workspace/planemo/docs/writing/seqtk_seq.cwl
      ...
    Applying linter requirements_in_conda... CHECK
    .. INFO: Requirement [seqtk@1.2] matches target in best practice Conda channel [https://conda.anaconda.org/bioconda/osx-64].

.. note:: You can download a more complete version of the CWL seqtk_ seq from the Planemo
    tutorial using the command::

        $ planemo project_init --template=seqtk_complete_cwl seqtk_example
        $ cd seqtk_example

We can verify these tool requirements install with the ``conda_install`` command. With
its default parameters ``conda_install`` processes tools and creates isolated environments
for their declared `Software Requirements`_ (mirroring what can be done in production with 
cwltool_ and Toil_).

::

    $ planemo conda_install seqtk_seq.cwl
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
will configure cwltool during testing to reuse this environment. If you wish to interactively explore
the resulting enviornment to explore the installed tool or produce test data the output
of the ``conda_env`` command can be sourced.

::

    $ . <(planemo conda_env seqtk_seq.cwl)
    Deactivate environment with conda_env_deactivate
    (seqtk_seq) $ which seqtk
    /home/planemo/miniconda3/envs/jobdepsiJClEUfecc6d406196737781ff4456ec60975c137e04884e4f4b05dc68192f7cec4656/bin/seqtk
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

Here is a portion of the output from the testing command ``planemo test seqtk_seq.cwl``
demonstrating using this tool.

::

    $ planemo test --no-container --engine cwltool seqtk_seq.cwl
    Enable beta testing mode for testing.
    cwltool INFO: /Users/john/workspace/planemo/.venv/bin/planemo 1.0.20170828135420
    cwltool INFO: Resolved '/Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl' to 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl'
    cwltool INFO: [job seqtk_seq.cwl] /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpaDQ1nK$ seqtk \
        seq \
        -a \
        /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpJtPKCr/stg24cf7e67-5ca6-44a4-a46b-26cbe104e1d4/2.fastq > /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpaDQ1nK/out
    cwltool INFO: [job seqtk_seq.cwl] completed success
    cwltool INFO: Final process status is success
    galaxy.tools.parser.factory INFO: Loading CWL tool - this is experimental - tool likely will not function in future at least in same way.
    All 1 test(s) executed passed.
    seqtk_seq_0: passed

Since ``seqtk`` isn't on the path and we did not use a container, we can see the SoftwareRequirement 
resolution was successful and it found the environment we previously installed with ``conda_install``.

.. include:: _writing_conda_search.rst

----------------------------------------------------------------
Exercise - Leveraging Bioconda
----------------------------------------------------------------

Use the ``project_init`` command to download this exercise.

::

    $ planemo project_init --template conda_exercises_cwl conda_exercises
    $ cd conda_exercises/exercise1
    $ ls 
    pear.cwl              test-data

This project template contains a few exercises. The first uses a CWL tool for
`PEAR - Paired-End reAd mergeR <http://sco.h-its.org/exelixis/web/software/pear/>`__.
This tool however has no ``SoftwareRequirement`` or container annotations and so will
not work properly without modification.

1. Run ``planemo test pear.cwl`` to verify the tool does not function
   without dependencies defined.
2. Use ``--conda_requirements`` flag with ``planemo lint`` to verify it does
   indeed lack requirements.
3. Use ``planemo conda_search`` or the Anaconda_ website to search for the
   correct package and version in a best practice channel.
4. Update ``pear.cwl`` with the correct ``SoftwareRequirement`` hints.
5. Re-run the ``lint`` command from above to verify the tool now has the
   correct dependency definition.
6. Re-run the ``test`` command from above to verify the tool test now
   works properly.


.. include:: _writing_conda_new.rst

::

    $ planemo project_init --template conda_exercises_cwl conda_exercises
    $ cd conda_exercises/exercise2
    $ ls 
    fleeqtk_seq.cwl      fleeqtk_seq_tests.yml         test-data

.. include:: _writing_conda_fleeqtk.rst

1. Clone and branch Bioconda_.
2. Build a recipe for fleeqtk version 1.3. You may wish to use ``conda skeleton``, start from
   scratch, or copy the recipe of seqtk and work from there - any of these strategies should work.  
3. Use ``conda build`` or Bioconda tooling to build the recipe.
4. Run ``planemo conda_install --conda_use_local fleeqtk_seq.cwl`` to verify the resulting package
   can be built into a Galaxy environment.
5. Run ``planemo test fleeqtk_seq.cwl`` to verify the resulting package works as expected.

.. note: The planemo flag ``--conda_use_local`` causes Planemo to use locally built
     packages during dependency resolution and related commands.

.. include:: _writing_conda_recipe_complete.rst

.. _Software Requirements: https://www.commonwl.org/v1.0/CommandLineTool.html#SoftwareRequirement
.. _seqtk: https://github.com/lh3/seqtk
.. _fleeqtk: https://github.com/jmchilton/fleeqtk
.. _Bioconda: https://github.com/bioconda/bioconda-recipes
.. _Conda: https://conda.io/docs/
.. _Anaconda: https://anaconda.org/
.. _cwltool: https://github.com/common-workflow-language/cwltool
.. _Toil: https://github.com/BD2KGenomics/toil
.. _introductory tutorial: http://planemo.readthedocs.io/en/latest/writing_cwl.html
