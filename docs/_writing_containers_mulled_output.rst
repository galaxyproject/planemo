
As the output indicates, this command built the container named
``quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0``.
This is the same namespace / URL that would be used if or when published by
the BioContainers_ project.

.. note:: The first part of this ``mulled-v2`` hash is a hash of the package names
    that went into it, the second the packages used and build number. Check out
    the `Multi-package Containers <http://biocontainers.pro/multi-package-containers/>`__
    web application to explore best practice channels and build such hashes.

We can see this new container when running the Docker command ``images`` and
explore the new container interactively with ``docker run``.

::

    $ docker images
    REPOSITORY                                                                 TAG                                          IMAGE ID            CREATED              SIZE
    quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40   03dc1d2818d9de56938078b8b78b82d967c1f820-0   a740fe1e6a9e        16 hours ago         104 MB
    quay.io/biocontainers/seqtk                                                1.2--0                                       10bc359ebd30        2 days ago           7.34 MB
    continuumio/miniconda                                                      latest                                       6965a4889098        3 weeks ago          437 MB
    bgruening/busybox-bash                                                     0.1                                          3d974f51245c        9 months ago         6.73 MB
    $ docker run -i -t quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0 /bin/bash
    bash-4.2# which samtools
    /usr/local/bin/samtools
    bash-4.2# which bwa
    /usr/local/bin/bwa
