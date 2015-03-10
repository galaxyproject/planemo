Collections
==============================

Galaxy has the concept of dataset collections to group together and operate
over them as single units with tools and in workflows.

Galaxy collections are hierarchical and composed from two simple collection
types - ``list`` and ``paired``.

A ``list`` is a simple a collection of datasets (or other collections) where
each element has an ``identifier``. Unlike Galaxy dataset names which are
transformed throughout complex analyses - the ``identifier`` is generally
perserved and can be used for concepts such ``sample`` name that one wants to
perserve the sample name in earlier mapping steps of a workflow and use it
during reduction steps and reporting later in workflows.

The ``paired`` collection type is much simpler and more specific to sequencing
applications. Each ``paired`` collection consists of a ``forward`` and
``reverse`` dataset.

.. note:: Read more about creating and managing collections on the `Galaxy Wiki <https://wiki.galaxyproject.org/Histories#Dataset_Collections>`__.

Composite types include for instance the ``list:paired`` collection type -
which represents a list of dataset pairs. In this case, instead of each
dataset having a list idenifier, each pair of datasets does.

Many Galaxy tools can in conjuction with collections used without
modification. Galaxy users can take a collection and `map over` any tool that
consumes individual datasets. For instance, early in typical bioinformatics
workflows you may have steps that filter raw data, convert to standard
formats, perform QC on individual files - users can take lists, pairs, or
lists of paired datasets and map over such tools that consume individual
files. Galaxy will then run the tool once for each dataset in the collection
and for each output of that tool Galaxy will rebuild a new collection with the
same `identifier` structure (so sample name or forward/reverse structure is
perserved).

Tools can also consume collections if they must or should process multiple
files at once. We will discuss three cases - consuming pairs of datasets,
consuming lists, and consuming arbitrary collections.

.. warning:: If you find yourself consuming a collection of files and calling the underlying application multiple times within the tool command block, you are likely doing something wrong. Just process and pair or a single dataset and allow the user to map over the collection.

Dataset collections are in their infancy - so for tools which process datasets
the recommended best practice is to allow users to either supply paired
collections or two individual datasets. Furthermore, many tools which process
pairs of datasets can also process single datasets. The following
``conditional`` captures this idiom.

::

    <conditional name="fastq_input">
      <param name="fastq_input_selector" type="select" label="Single or Paired-end reads" help="Select between paired and single end data">
        <option value="paired">Paired</option>
        <option value="single">Single</option>
        <option value="paired_collection">Paired Collection</option>
        <option value="paired_iv">Paired Interleaved</option>
      </param>
      <when value="paired">
        <param name="fastq_input1" type="data" format="fastqsanger" label="Select first set of reads" help="Specify dataset with forward reads"/>
        <param name="fastq_input2" type="data" format="fastqsanger" label="Select second set of reads" help="Specify dataset with reverse reads"/>
      </when>     
      <when value="single">
        <param name="fastq_input1" type="data" format="fastqsanger" label="Select fastq dataset" help="Specify dataset with single reads"/>
      </when>
      <when value="paired_collection">
        <param name="fastq_input" format="fastqsanger" type="data_collection" collection_type="paired" label="Select a paired collection" />
      </when>
    </conditional>

This introduces a new ``param`` type - ``data_collection``. The optional
attribute ``collection_type`` can be specified to specify which kinds of
collections are appropriate for this input. Additional ``data`` attributes
such as ``format`` can be specified to further restrict valid collections.
Here we specified that both items of the paired collection must be of datatype
``fastqsanger``.

In Galaxy's ``command`` block, the individual datasets can be accessed using
``$fastq_input1.forward`` and ``$fastq_input1.reverse``. If processing
arbitrary collection types an array syntax can also be used (e.g.
``$fastq_input['forward']``).

Some example tools which consume paired datasets include:

 - `collection_paired_test <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_paired_test.xml>`__ (minimal test tool in Galaxy test suite)
 - `Bowtie 2 <https://github.com/galaxyproject/tools-devteam/blob/master/tools/bowtie2/bowtie2_wrapper.xml>`__
 - `BWA MEM <https://github.com/galaxyproject/tools-devteam/blob/master/tools/bwa/bwa-mem.xml>`__
 - `Tophat <https://github.com/galaxyproject/tools-devteam/blob/master/tools/tophat2/tophat2_wrapper.xml>`__

-------------------------------
Processing Lists (Reductions)
-------------------------------

The ``data_collection`` parameter type can specify a ``collection_type`` or
``list`` but whenever possible, it is actually better to not explicitly
consume lists as a tool author. Parameters of type ``data`` can include a
``multiple="True"`` attribute to allow many datasets to be selected
simultaneously. While the default UI will then have Galaxy users pick
individual datsets, they can easily substitute a collections the tool can
process both as individual datasets. This has the benefit of allowing tools to
process either individual datasets or collections.

::

    <param type="data" name="inputs" label="Input BAM(s)" format="bam" multiple="true" />

The ``command`` tag can use ``for`` `loops <http://www.cheetahtemplate.org/docs/users_guide_html/users_guide.html#SECTION0001010000000000000000>`__ to build command lines using these parameters.

For instance:

::

    #for $input in $inputs 
    --input "$input"
    #end for

or using the single form of this expression:

::

    #for $input in $inputs# $input #end for#

Will produce command strings with an argument for each input (e.g. ``--input
"/path/to/input1" --input "/path/to/input2"``). Other programs may require all
inputs to be supplied in a single parameter. This can be accomplished using
the idiom:

::

    --input "${",".join(map(str, $inputs))}"

TODO: test that statement!

TODO: mention identifiers

Some example tools which consume collections include:

 - `multi_data_param <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/multi_data_param.xml>`__ (small test tool in Galaxy test suite)
 - `cuffmerge <https://github.com/galaxyproject/tools-devteam/blob/master/tool_collections/cufflinks/cuffmerge/cuffmerge_wrapper.xml>`__
 - `unionBedGraphs <https://github.com/galaxyproject/tools-iuc/blob/master/tools/bedtools/unionBedGraphs.xml>`__

Also see the tools-devteam repository `Pull Request #20 <https://github.com/galaxyproject/tools-devteam/pull/20>`__ modifying the cufflinks suite of tools for collection compatible reductions.

-------------------------------
Processing Collections
-------------------------------

Some example tools which consume collections include:

 - `collection_nested_test <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_nested_test.xml>`_ (small test tool demonstrating consumption of nested collections)

----------------------
Further Reading
----------------------

 - Galaxy Community Conference Talk by John Chilton [`Slides <http://bit.ly/gcc2014workflows>`__][`Video <http://jh.hosted.panopto.com/Panopto/Pages/Viewer.aspx?id=f626696c-e68e-4aa4-870b-f224aa60c47a>`__].
 - `Creating and Managing Collections <https://wiki.galaxyproject.org/Histories#Dataset_Collections>`__
 - `Pull Request #386 <https://bitbucket.org/galaxy/galaxy-central/pull-request/386/dataset-collections-initial-models-api>`__ (the initial implementation)
 - `Pull Request #634 <https://bitbucket.org/galaxy/galaxy-central/pull-request/634/allow-tools-to-explicitly-produce-dataset>`__ (implementing ability for tools to explicitly output collections)
