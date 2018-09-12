
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
    /Users/john/miniconda3/bin/conda search --override-channels --channel iuc --channel conda-forge --channel bioconda --channel defaults '*seqt*'
    Loading channels: done
    # Name                  Version           Build  Channel
    bioconductor-htseqtools          1.26.0        r3.4.1_0  bioconda
    bioconductor-seqtools          1.10.0        r3.3.2_0  bioconda
    bioconductor-seqtools          1.10.0        r3.4.1_0  bioconda
    bioconductor-seqtools          1.12.0        r3.4.1_0  bioconda
    seqtk                       r75               0  bioconda
    seqtk                       r82               0  bioconda
    seqtk                       r82               1  bioconda
    seqtk                       r93               0  bioconda
    seqtk                       1.2               0  bioconda
    seqtk                       1.2               1  bioconda

.. note:: The Planemo command ``conda_search`` is a light wrapper around the underlying
   ``conda search`` command but configured to use the same channels and other options as
   Planemo and Galaxy. The following Conda command would also work to search::

       $ $HOME/miniconda3/bin/conda -c iuc -c conda-forge -c bioconda '*seqt*'

   For Conda versions 4.3.X or less, the search invocation would be something a bit
   different::

       $ $HOME/miniconda3/bin/conda -c iuc -c conda-forge -c bioconda seqt


Alternatively the Anaconda_ website can be used to search for packages. Typing ``seqtk``
into the search form on that page and clicking the top result will bring on to `this page
<https://anaconda.org/bioconda/seqtk>`__ with information about the Bioconda package.

When using the website to search though, you need to aware of what channel you are using. By
default, Planemo and Galaxy will search a few different Conda channels. While it is possible
to configure a local Planemo or Galaxy to target different channels - the current best practice
is to add tools to the existing channels.

The existing channels include:

* Bioconda (`github <https://github.com/bioconda/bioconda-recipes>`__ | `conda <https://anaconda.org/bioconda>`__) - best practice channel for various bioinformatics packages.
* Conda-Forge (`github <https://github.com/conda-forge/staged-recipes>`__ | `conda <https://anaconda.org/conda-forge>`__) - best practice channel for general purpose and widely useful computing packages and libraries.
* iuc (`github <https://github.com/galaxyproject/conda-iuc>`__ | `conda <https://anaconda.org/iuc>`__) - best practice channel for other more Galaxy specific packages.
