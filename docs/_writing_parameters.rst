Simple Parameters
====================================================

We have built a tool wrapper for the ``seqtk seq`` command - but this tool
actually has additional options that we may wish to expose the Galaxy user.

Lets take a few of the parameters from the help command and build Galaxy
``param`` blocks to stick in the tool's ``inputs`` block.

::

    -V        shift quality by '(-Q) - 33'

In the previous section we saw ``param`` block of type ``data`` for input
files, but there are many different kinds of parameters one can use.
Flag parameters such as the above ``-V`` parameter are frequently
represented by ``boolean`` parameters in Galaxy tool XML.

.. code-block:: xml

    <param name="shift_quality" type="boolean" label="Shift quality"
           truevalue="-V" falsevalue=""
           help="shift quality by '(-Q) - 33' (-V)" />

We can then stick ``$shift_quality`` in our ``command`` block and if the user
has selected this option it will be expanded as ``-V`` (since we have defined
this as the ``truevalue``). If the user hasn't selected this option
``$shift_quality`` will just expand as an empty string and not affect the
generated command line.

Now consider the following ``seqtk seq`` parameters:

::

    -q INT    mask bases with quality lower than INT [0]
    -X INT    mask bases with quality higher than INT [255]

These can be translated into Galaxy parameters as:

.. code-block:: xml

    <param name="quality_min" type="integer" label="Mask bases with quality lower than"
           value="0" min="0" max="255" help="(-q)" />
    <param name="quality_max" type="integer" label="Mask bases with quality higher than"
           value="255" min="0" max="255" help="(-X)" />

These can be add to the command tag as ``-q $quality_min -X $quality_max``.

At this point the tool would look like:

.. literalinclude:: writing/seqtk_seq_v4.xml
   :language: xml
   :emphasize-lines: 7-10,14-20


Conditional Parameters
====================================================

The previous parameters were simple because they always appeared, now consider.

::

    -M FILE   mask regions in BED or name list FILE [null]

We can mark this ``data`` type ``param`` as optional by adding the attribute
``optional="true"``.

.. code-block:: xml

    <param name="mask_regions" type="data" label="Mask regions in BED"
           format="bed" help="(-M)" optional="true" />

Then instead of just using ``$mask_regions`` directly in the ``command``
block, one can wrap it in an ``if`` statement (because tool XML files
support `Cheetah <https://cheetahtemplate.org/users_guide/index.html>`_).

::

    #if $mask_regions
    -M '$mask_regions'
    #end if

Next consider the parameters:

::

    -s INT    random seed (effective with -f) [11]
    -f FLOAT  sample FLOAT fraction of sequences [1]

In this case, the ``-s`` random seed parameter should only be seen or used if
the sample parameter is set. We can express this using a ``conditional``
block.

.. code-block:: xml

    <conditional name="sample">
        <param name="sample_selector" type="boolean" label="Sample fraction of sequences" />
        <when value="true">
            <param name="fraction" label="Fraction" type="float" value="1.0"
                   help="(-f)" />
            <param name="seed" label="Random seed" type="integer" value="11"
                   help="(-s)" />
        </when>
        <when value="false">
        </when>
    </conditional>

In our command block, we can again use an ``if`` statement to include these
parameters.

::

    #if $sample.sample_selector
    -f $sample.fraction -s $sample.seed
    #end if

Notice we must reference the parameters using the ``sample.`` prefix since
they are defined within the ``sample`` conditional block.

The newest version of this tool is now

.. literalinclude:: writing/seqtk_seq_v5.xml
   :language: xml
   :emphasize-lines: 10-16,28-40

For tools like this where there are many options but in the most uses the defaults
are preferred - a common idiom is to break the parameters into simple and
advanced sections using a ``conditional``.

Updating this tool to use that idiom might look as follows.

.. literalinclude:: writing/seqtk_seq_v6.xml
   :language: xml
   :emphasize-lines: 7-18,23-30,51-52
