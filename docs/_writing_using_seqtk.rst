This guide is going to demonstrate building up tools for commands from Heng
Li's Seqtk_ package - a package for processing sequence data in FASTA_ and
FASTQ_ files.

To get started let's install Seqtk. Here we are going to use ``conda`` to
install Seqtk - but however you obtain it should be fine.

::

    $ conda config --add channels r
    $ conda config --add channels bioconda
    $ conda install seqtk
        ... seqtk installation ...
    $ seqtk seq
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

Next we will download an example FASTQ file and test out the a simple Seqtk
command - ``seq`` which converts FASTQ files into FASTA.

::

    $ wget https://raw.githubusercontent.com/galaxyproject/galaxy-test-data/master/2.fastq
    $ seqtk seq -A 2.fastq > 2.fasta
    $ cat 2.fasta
    >EAS54_6_R1_2_1_413_324
    CCCTTCTTGTCTTCAGCGTTTCTCC
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    >EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG

.. _Seqtk: https://github.com/lh3/seqtk
.. _FASTA: https://en.wikipedia.org/wiki/FASTA_format
.. _FASTQ: https://en.wikipedia.org/wiki/FASTQ_format
