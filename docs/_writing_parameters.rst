Adding Parameters
====================================================

We have built a tool wrapper for the ``seqtk seq`` command - but this tool
actually has additional options that we may wish to expose the Galaxy user.

Lets take a few of the parameters from the help command and build Galaxy ``param`` blocks to stick in the tool's ``inputs`` block.

::

    -V        shift quality by '(-Q) - 33'

In the previous section we saw ``param`` block of type ``data`` for input
files, but there are many different kinds of parameters one can use. Flag parameters such as the above ``-V`` parameter are frequently represented by ``boolean`` parameters in Galaxy tool XML.

::

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

::

    <param name="quality_min" type="integer" label="Mask bases with quality lower than" 
           value="0" min="0" max="255" help="(-q)" />
    <param name="quality_max" type="integer" label="Mask bases with quality higher than" 
           value="255" min="0" max="255" help="(-X)" />

These can be add to the command tag as ``-q $quality_min -X $quality_max``.

At this point the tool would look like:

.. literalinclude:: writing/seqtk_seq_v4.xml
   :language: xml
