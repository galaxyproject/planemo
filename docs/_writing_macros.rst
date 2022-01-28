Macros - Reusable Elements
==============================

Frequently, tools may require the same XML fragments be repeated in a file
(for instance similar conditional branches, repeated options, etc...) or
between tools in the same repository (for instance, nearly all of the GATK
tools contain the same standard options). Galaxy tools have a macroing system
to address this problem.

-----------------------------------------
Direct XML Macros
-----------------------------------------

The following examples are taken from `Pull Request 129
<https://bitbucket.org/galaxy/galaxy-central/pull-requests/129>`__ the initial
implementation of macros. Prior to to the inclusion of macros, the tophat2
wrapper defined several outputs each which had the following identical actions
block associated with them:

.. code-block:: xml

   <actions>
     <conditional name="refGenomeSource.genomeSource">
       <when value="indexed">
         <action type="metadata" name="dbkey">
           <option type="from_data_table" name="tophat2_indexes" column="1" offset="0">
             <filter type="param_value" column="0" value="#" compare="startswith" keep="False"/>
             <filter type="param_value" ref="refGenomeSource.index" column="0"/>
           </option>
         </action>
       </when>
       <when value="history">
         <action type="metadata" name="dbkey">
           <option type="from_param" name="refGenomeSource.ownFile" param_attribute="dbkey" />
         </action>
       </when>
     </conditional>
   </actions>

To reuse this action definition, first a macros section has been defined in the tophat2_wrpper.xml file.

.. code-block:: xml

   <tool>
      ...
      <macros>
        <xml name="dbKeyActions">
          <action><!-- Whole big example above. -->
            ....
          </action>
        </xml>
      </macros>

With this in place, each output data element can include this block using the expand XML element as follows:

.. code-block:: xml

   <outputs>
       <data format="bed" name="insertions" label="${tool.name} on ${on_string}: insertions" from_work_dir="tophat_out/insertions.bed">
           <expand macro="dbKeyActions" />
       </data>
       <data format="bed" name="deletions" label="${tool.name} on ${on_string}: deletions" from_work_dir="tophat_out/deletions.bed">
         <expand macro="dbKeyActions" />
       </data>
       <data format="bed" name="junctions" label="${tool.name} on ${on_string}: splice junctions" from_work_dir="tophat_out/junctions.bed">
         <expand macro="dbKeyActions" />
       </data>
       <data format="bam" name="accepted_hits" label="${tool.name} on ${on_string}: accepted_hits" from_work_dir="tophat_out/accepted_hits.bam">
         <expand macro="dbKeyActions" />
       </data>
   </outputs>

This has reduced the size of the XML file by dozens of lines and reduces the long term maintenance associated with copied and pasted code.

-----------------------------------------
Imported Macros
-----------------------------------------

The ``macros`` element described above, can also contain any number of
``import`` elements. This allows a directory/repository of tool XML files to
contain shared macro definitions that can be used by any number of actual tool
files in that directory/repository.

Revisiting the tophat example, all three tophat wrappers (``tophat_wrapper.xml``,
``tophat_color_wrapper.xml``, ``and tophat2_wrapper.xml``) shared some common
functionality. To reuse XML elements between these files, a
``tophat_macros.xml`` file was added to that directory.

The following block is a simplified version of that macros file's contents:

.. code-block:: xml

   <macros>
     <xml name="own_junctionsConditional">
       <conditional name="own_junctions">
         <param name="use_junctions" type="select" label="Use Own Junctions">
           <option value="No">No</option>
           <option value="Yes">Yes</option>
         </param>
         <when value="Yes">
           <conditional name="gene_model_ann">
             <param name="use_annotations" type="select" label="Use Gene Annotation Model">
               <option value="No">No</option>
               <option value="Yes">Yes</option>
             </param>
             <when value="No" />
             <when value="Yes">
               <param format="gtf,gff3" name="gene_annotation_model" type="data" label="Gene Model Annotations" help="TopHat will use the exon records in this file to build a set of known splice junctions for each gene, and will attempt to align reads to these junctions even if they would not normally be covered by the initial mapping."/>
             </when>
           </conditional>
           <expand macro="raw_juncsConditional" />
           <expand macro="no_novel_juncsParam" />
         </when>
         <when value="No" />
       </conditional> <!-- /own_junctions -->
     </xml>
     <xml name="raw_juncsConditional">
       <conditional name="raw_juncs">
         <param name="use_juncs" type="select" label="Use Raw Junctions">
           <option value="No">No</option>
           <option value="Yes">Yes</option>
         </param>
         <when value="No" />
         <when value="Yes">
           <param format="interval" name="raw_juncs" type="data" label="Raw Junctions" help="Supply TopHat with a list of raw junctions. Junctions are specified one per line, in a tab-delimited format. Records look like: [chrom] [left] [right] [+/-] left and right are zero-based coordinates, and specify the last character of the left sequenced to be spliced to the first character of the right sequence, inclusive."/>
         </when>
       </conditional>
     </xml>
     <xml name="no_novel_juncsParam">
       <param name="no_novel_juncs" type="select" label="Only look for supplied junctions">
         <option value="No">No</option>
         <option value="Yes">Yes</option>
       </param>
     </xml>
   </macros>

Any tool definition in that directory can use the macros contained therein once imported as shown below.

.. code-block:: xml

   <tool>
     ...
     <macros>
       <import>tophat_macros.xml</import>
     </macros>
     ...
     <inputs>
       <expand macro="own_junctionsConditional" />
       ...
     </inputs>
     ...
   </tool>

This example also demonstrates that macros may themselves expand macros.

-------------------------------------------
Parameterizing XML Macros (with ``yield``)
-------------------------------------------

In some cases, tools may contain similar though not exact same definitions. Some parameterization can be performed by declaring expand elements with child elements and expanding them in the macro definition with a yield element.

For instance, previously the tophat wrapper contained the following definition:

.. code-block:: xml

   <conditional name="refGenomeSource">
     <param name="genomeSource" type="select" label="Will you select a reference genome from your history or use a built-in index?" help="Built-ins were indexed using default options">
       <option value="indexed">Use a built-in index</option>
       <option value="history">Use one from the history</option>
     </param>
     <when value="indexed">
       <param name="index" type="select" label="Select a reference genome" help="If your genome of interest is not listed, contact the Galaxy team">
         <options from_data_table="tophat_indexes_color">
           <filter type="sort_by" column="2"/>
           <validator type="no_options" message="No indexes are available for the selected input dataset"/>
         </options>
       </param>
     </when>
     <when value="history">
       <param name="ownFile" type="data" format="fasta" metadata_name="dbkey" label="Select the reference genome" />
     </when>  <!-- history -->
   </conditional>  <!-- refGenomeSource -->

and the tophat2 wrapper contained the highly analogous definition:

.. code-block:: xml

   <conditional name="refGenomeSource">
     <param name="genomeSource" type="select" label="Will you select a reference genome from your history or use a built-in index?" help="Built-ins were indexed using default options">
       <option value="indexed">Use a built-in index</option>
       <option value="history">Use one from the history</option>
     </param>
     <when value="indexed">
       <param name="index" type="select" label="Select a reference genome" help="If your genome of interest is not listed, contact the Galaxy team">
         <options from_data_table="tophat2_indexes_color">
           <filter type="sort_by" column="2"/>
           <validator type="no_options" message="No indexes are available for the selected input dataset"/>
         </options>
       </param>
     </when>
     <when value="history">
       <param name="ownFile" type="data" format="fasta" metadata_name="dbkey" label="Select the reference genome" />
     </when>  <!-- history -->
   </conditional>  <!-- refGenomeSource -->


These blocks differ only in the from_data_table attribute on the options element. To capture this pattern, tophat_macros.xml contains the following macro definition:


.. code-block:: xml

   <xml name="refGenomeSourceConditional">
     <conditional name="refGenomeSource">
       <param name="genomeSource" type="select" label="Use a built in reference genome or own from your history" help="Built-ins genomes were created using default options">
         <option value="indexed" selected="True">Use a built-in genome</option>
         <option value="history">Use a genome from history</option>
       </param>
       <when value="indexed">
         <param name="index" type="select" label="Select a reference genome" help="If your genome of interest is not listed, contact the Galaxy team">
           <yield />
         </param>
       </when>
       <when value="history">
         <param name="ownFile" type="data" format="fasta" metadata_name="dbkey" label="Select the reference genome" />
       </when>  <!-- history -->
     </conditional>  <!-- refGenomeSource -->
   </xml>

Notice the yield statement in lieu of an options declaration. This allows the nested options element to be declared when expanding the macro:

The following expand declarations have replaced the original conditional elements.

.. code-block:: xml

   <expand macro="refGenomeSourceConditional">
     <options from_data_table="tophat_indexes">
       <filter type="sort_by" column="2"/>
       <validator type="no_options" message="No genomes are available for the selected input dataset"/>
     </options>
   </expand>

.. code-block:: xml

   <expand macro="refGenomeSourceConditional">
     <options from_data_table="tophat2_indexes">
       <filter type="sort_by" column="2"/>
       <validator type="no_options" message="No genomes are available for the selected input dataset"/>
     </options>
   </expand>

From Galaxy release 22.01 named yields are supported. That is, if the macro contains
``<yield name="xyz"/>`` it is replaced by the content of the ``token`` tag with the same name.
Token tags need to be direct children of the ``expand`` element. This is useful if different
parts of the macro should be parametrized. 

In the following example two named yield and one unnamed yield are used to
parametrize the options of the select of a conditional, the options of another
select, and additional when block(s). Given the following macro:

.. code-block:: xml
   :emphasize-lines: 6,10,14

   <xml name="named_yields_example">
     <conditional>
       <param type="select">
         <option value="a">A</option>
         <option value="b">B</option>
         <yield name="more_options"/>
       </param>
       <when value="a">
         <param name="aselect" type="select">
           <yield />
         </param>
       </when>
       <when value="b"/>
       <yield name="more_whens">
     </conditional>
   </xml>

and expanding the macro in the following way:

.. code-block:: xml
   :emphasize-lines: 2,5,8,12

   <expand macro="named_yields_example">
     <token name="more_options">
       <option value="c">C</option>
     </token>
     <token name="more_whens">
       <when value="c">
         <param type="select">
           <yield />
         </param>
       </when>
     </token>
     <options from_data_table="tophat2_indexes">
       <filter type="sort_by" column="2"/>
       <validator type="no_options" message="No genomes are available for the selected input dataset"/>
     </options>
   </expand>

we get the following expanded definition:

.. code-block:: xml
   :emphasize-lines: 6,9-14,17-24

   <xml name="named_yields_example">
     <conditional>
       <param type="select">
         <option value="a">A</option>
         <option value="b">B</option>
         <option value="c">C</option>
       </param>
       <when value="a">
         <param name="aselect" type="select">
           <options from_data_table="tophat2_indexes">
             <filter type="sort_by" column="2"/>
             <validator type="no_options" message="No genomes are available for the selected input dataset"/>
           </options>
         </param>
       </when>
       <when value="b"/>
       <when value="c">
         <param type="select">
           <options from_data_table="tophat2_indexes">
             <filter type="sort_by" column="2"/>
             <validator type="no_options" message="No genomes are available for the selected input dataset"/>
           </options>
         </param>
       </when>
     </conditional>
   </xml>

Named yields are replaced in the order of the tokens defined in the ``expand``
tag. Unamed yields are replaced after all named tokens have been replaced (by
the non-token child elements of the expand tag). If there are named
yields that have no corresponding token, then they are treated like unnamed
yields. Note that, unnamed and named tokens can be used multiple times in a
macro, then each occurence is preplaced by the corresponding content definedin
the ``expand``.

Further, note that the order of the replacements offers some possibilities
to achieve recursion-like replacements, since a token may contain further named
or unnamed ``yield`` tags (see for instance the ``yield`` tag contained in the
named token ``more_whens``).

-----------------------------------------
Parameterizing XML Macros (with tokens)
-----------------------------------------

In addition to using ``yield`` blocks, there is another way to parameterize
macros by specifying:

- required parameters as comma-separated list of parameter names using the
  ``tokens`` attribute (e.g. ``tokens="foo,bar"``) of the ``xml`` element and then using
  ``@FOO@`` and ``@BAR@`` in the macro definition;
- optional parameters as ``token_xyz="default_value"`` attributes of the ``xml`` 
  element, and then using ``@XYZ@`` in the macro definition.

.. code-block:: xml

   <macros>
     <xml name="color" tokens="varname" token_default_color="#00ff00" token_label="Pick a color">
         <param name="@VARNAME@" type="color" label="@LABEL@" value="@DEFAULT_COLOR@" />
     </xml>
   </macros>

This defines a macro with a required parameter ``varname`` and two optional
parameters ``default_color`` and ``label``.  When invoking this macro, you can
pass values for those parameters and produce varying results.

.. code-block:: xml

   <inputs>
       <expand macro="color" varname="myvar" default_color="#ff0000" />
       <expand macro="color" varname="c2" default_color="#0000ff" label="Choose a different color" />
   </inputs>

The attributes passed to the macro definition will be filled in (or defaults used if not provided).
Effectively this yields:

.. code-block:: xml

   <inputs>
       <param name="myvar" type="color" label="Pick a color" value="#ff0000" />
       <param name="c2" type="color" label="Choose a different color" value="#0000ff" />
   </inputs>

Macro tokens can be used in the text content of tags, attribute values, and (with a
little trick also in attribute names). The problem is that the default delimiting character of
macro tokens is ``@`` and the XML must still be valid before processing the macros (and ``@``
is invalid in attribute names). Luckily the delimiting character(s) can be changed by adding
``token_quote`` to the macro definition:

.. code-block:: xml

   <macros>
     <xml name="color" tokens="attr,attr_value" token_quote="__" token_label="label">
         <param __VAR__="__ATTR_VALUE__" label="__LABEL__"/>
     </xml>
   </macros>

Note that, this forbids to use tokens with the name ``quote``.

-----------------------------------------
Macro Tokens
-----------------------------------------

You can use

.. code-block:: xml

   <token name="@IS_PART_OF_VCFLIB@">is a part of VCFlib toolkit developed by Erik Garrison (https://github.com/ekg/vcflib).</token>

and then call the token within any element of the file like this:

.. code-block:: shell

   Vcfallelicprimitives @IS_PART_OF_VCFLIB@
