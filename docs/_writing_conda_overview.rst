.. note:: *Why Conda?*

    Many different package managers could potentially be targeted here, but we focus on Conda_
    for a few key reasons.

    * No compilation at install time - binaries with their dependencies and libraries
    * Support for all operating systems
    * Easy to manage multiple versions of the same recipe
    * HPC-ready: no root privileges needed
    * Easy-to-write YAML recipes
    * Viberant communities

.. note:: **Conda Terminology**

    .. figure:: http://galaxyproject.github.io/training-material/topics/dev/images/miniconda_vs_anaconda.png
       :alt: Diagram describing the relationship between Conda, Miniconda, and Anaconda.

    Conda *recipes* build *packages* that are published to *channels*.

Planemo is setup to target a few channels by default, these include ``iuc``, ``bioconda``,
``conda_forge``, ``defaults`` - the whole dependency management scheme outlined here works a lot
better if packages can be found in one of these "best practice" channels.
