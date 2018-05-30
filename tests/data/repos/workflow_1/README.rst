# Directory is modified BLAST example from Peter Cock
# https://raw.githubusercontent.com/peterjc/galaxy_blast/master/workflows/blast_top_hit_species

Introduction
============

Galaxy is a web-based platform for biological data analysis, supporting
extension with additional tools (often wrappers for existing command line
tools) and datatypes. See http://www.galaxyproject.org/ and the public
server at http://usegalaxy.org for an example.

The NCBI BLAST suite is a widely used set of tools for biological sequence
comparison. It is available as standalone binaries for use at the command
line, and via the NCBI website for smaller searches. For more details see
http://blast.ncbi.nlm.nih.gov/Blast.cgi

This is an example workflow using the Galaxy wrappers for NCBI BLAST+,
see https://github.com/peterjc/galaxy_blast


Galaxy workflow for counting species of top BLAST hits 
======================================================

This Galaxy workflow (file ``blast_top_hit_species.ga``) is intended for an
initial assessment of a transcriptome assembly to give a crude indication of
any major contamination present based on the species of the top BLAST hit
of 1000 representative sequences.

.. image:: https://raw.githubusercontent.com/peterjc/galaxy_blast/master/workflows/blast_top_hit_species/blast_top_hit_species.png

In words, the workflow proceeds as follows:

1. Upload/import your transcriptome assembly or any nucleotide FASTA file.
2. Samples 1000 representative sequences, selected uniformly/evenly though
   the file.
3. Convert the sampled FASTA file into a three column tabular file.
4. Runs NCBI BLASTX of the sampled FASTA file against the latest NCBI ``nr``
   database (assuming this is already available setup on your local Galaxy
   under the alias ``nr``), requesting tabular output including the taxonomy
   fields, and at most one matching target sequence.
5. Remove any duplicate alignments (multiple HSPs for the same match).
6. Combine the filtered BLAST output with the tabular version of the 1000
   sequences to give a new tabular file with exactly 1000 lines, adding
   ``None`` for sequences missing a BLAST hit.
7. Count the BLAST species names in this file.
8. Sort the counts.

Finally we would suggest visualising the sorted tally table as a Pie Chart,
as in the example below.


Sample Data
===========

As an example, you can upload the transcriptome assembly of the nematode
*Nacobbus abberans* from Eves van den Akker *et al.* (2015),
https://doi.org/10.1093/gbe/evu171 using this URL:

http://nematode.net/Data/nacobbus_aberrans_transcript_assembly/N.abberans_reference_no_contam.zip

Running this workflow with a copy of the NCBI non-redundant ``nr`` database
from 16 Oct 2014 (which did **not** contain this *N. abberans* dataset) gave
the following results - note 609 out of the 1000 sequences gave no BLAST hit.

===== ==================
Count Subject Blast Name
----- ------------------
  609 None
  244 nematodes
   30 ascomycetes
   27 eukaryotes
    8 basidiomycetes
    6 aphids
    5 eudicots
    5 flies
  ... ...
===== ==================

As you might guess from	the filename ``N.abberans_reference_no_contam.fasta``,
this transcriptome assembly has already had obvious contamination removed.

At the time of writing, Galaxy's visualizations could not be included in
a workflow. You can generate a pie chart from the final count file using
the counts (c1) and labels (c2), like this:

.. image:: https://raw.githubusercontent.com/peterjc/galaxy_blast/master/workflows/blast_top_hit_species/N_abberans_piechart_mouseover.png

Note the nematode count in this image was shown as a mouse-over effect.


Disclaimer
==========

Species assignment by top BLAST hit is not suitable for any in depth
analysis. It is particularly prone to false positives where contaminants
in public datasets are mislabelled. See for example Ed Yong (2015),
"There's No Plague on the NYC Subway. No Platypuses Either.":

http://phenomena.nationalgeographic.com/2015/02/10/theres-no-plague-on-the-nyc-subway-no-platypuses-either/


Known Issues
============

Counts
------

This workflow uses the Galaxy "Count" tool (tool id ``Count1``) version
1.0.0, as shipped with the current stable release (Galaxy v15.03, i.e.
March 2015).

The updated "Count" tool version 1.0.1 included a fix not to remove spaces
in the fields being counted. In the example above, while the top hits are
not affected, minor entries like "cellular slime molds" are shown as
"cellularslimemolds" instead (look closely at the Pie Chart key).

The updated "Count" tool version 1.0.2 added a new option to sort the
output, which would allow skipping the final sorting step in the current
version of this workflow.

A future update to this workflow will use the revised "Count" tool, once
this is included in the next stable Galaxy release - or migrated to the
Galaxy Tool Shed.

NCBI nr database
----------------

The use of external datasets within Galaxy via the ``*.loc`` configuration
files undermines provenance tracking within Galaxy. This is exacerbated
by the lack of officially versioned BLAST database releases by the NCBI.

This workflow assumes that you have an entry ``nr`` in your ``blastdb_p.loc``
(the configuration file listing locally installed BLAST databases external
to Galaxy - consult the NCBI BLAST+ wrapper documentation for more details),
and that this points to a mirror of the latest NCBI "non-redundant" database
from ftp://ftp.ncbi.nlm.nih.gov/blast/db/

i.e. The workflow is intended to be used against the *latest* nr database,
and thus is not reproducible over the long term as the database changes.

Note that if your ``blastdb_p.loc`` is missing an entry ``nr`` then the
workflow should abort. However as of Galaxy v15.03 (March 2015) there is
a problem with how this is handled: https://trello.com/c/lkYlW14W/


Availability
============

This workflow is available from myExperiment:

http://www.myexperiment.org/workflows/4637

You can also download and/or install it from the main Galaxy Tool Shed:

http://toolshed.g2.bx.psu.edu/view/peterjc/blast_top_hit_species

Test releases (which should not normally be used) are on the Test Tool Shed:

http://testtoolshed.g2.bx.psu.edu/view/peterjc/blast_top_hit_species

Development is being done on github here:

https://github.com/peterjc/galaxy_blast/tree/master/workflows/blast_top_hit_species


Citation
========

Please cite the following paper (currently available as a preprint):

NCBI BLAST+ integrated into Galaxy.
P.J.A. Cock, J.M. Chilton, B. Gruening, J.E. Johnson, N. Soranzo
bioRxiv DOI: https://doi.org/10.1101/014043 (preprint)

You should also cite Galaxy, and the NCBI BLAST+ tools:

BLAST+: architecture and applications.
C. Camacho et al. BMC Bioinformatics 2009, 10:421.
DOI: https://doi.org/10.1186/1471-2105-10-421


Automated Installation
======================

Installation via the Galaxy Tool Shed should take care of the dependencies
on Galaxy tools including the NCBI BLAST+ wrappers and associated binaries.

However, this workflow requires a current version of the NCBI nr protein
BLAST database to be listed in ``blastdb_p.loc`` with the key ``nr`` (lower
case).


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.1.0  - Initial MyExperiment and Tool Shed release.
        - Targetting NCBI BLAST+ 2.2.29
======= ======================================================================


Developers
==========

This workflow is under source code control here:

https://github.com/peterjc/galaxy_blast/tree/master/workflows/blast_top_hit_species

To prepare the tar-ball for uploading to the Tool Shed, I use this::

    $ tar -cf blast_top_hit_species.tar.gz README.rst repository_dependencies.xml blast_top_hit_species.ga blast_top_hit_species.png N_abberans_piechart_mouseover.png

Check this::

    $ tar -tzf blast_top_hit_species.tar.gz
    README.rst
    repository_dependencies.xml
    blast_top_hit_species.ga
    blast_top_hit_species.png
    N_abberans_piechart_mouseover.png


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
