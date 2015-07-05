
::

    $ planemo t --galaxy_root=/path/to/galaxy
    ...
    All 1 test(s) executed passed.
    seqtk_seq[0]: passed

Now we can open Galaxy with the ``serve`` (or just ``s``).

::

    $ planemo s --galaxy_root=/path/to/galaxy
    ...
    serving on http://127.0.0.1:9090

Open up http://127.0.0.1:9090 in a web browser to view your new tool.
