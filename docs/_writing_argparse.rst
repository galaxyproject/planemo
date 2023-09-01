Creating Galaxy tools from python scripts using argparse
========================================================

If a python script should be wrapped that creates its command line interface
using the argparse library, then planemo can generate a good starting point.

All you need to do is to pass the path to the python source file containing
the argparse definition. Lets use ``tests/data/autopygen/autopygen_end_to_end_sub.py``
from the planemo tests as an example:

::

    $ planemo tool_init -i e2e -n e2e
                        --autopygen tests/data/autopygen/autopygen_end_to_end_sub.py
                        -- tool tool.xml

Here we only provided the ``id`` and ``name`` in addition to the python sources
and planemo creates a pretty tool xml files that can serve as a good starting point
to create a functional Galaxy tool in ``tool.xml``.
Additional information, e.g. requirements, help, etc, can and should be given with
additional parameters to planemo.

There are a few points that need to be edited in any case.  Most importantly
since in most cases ``argparse`` does not distinguish between input and output
files all file parameters will be rendered as input parameters of the tool.

Please open an issue if you have ideas on how to improve the generated tools:
https://github.com/galaxyproject/planemo/issues


