====================================================
Building Galaxy Tools (Using the Planemo Appliance)
====================================================

This tutorial is a gentle introduction to writing Galaxy tools using the
Planemo virtual appliance  (available as OVA_, Docker_ and Vagrant_). Check out
`these instructions <https://planemo.readthedocs.org/en/latest/appliance.html>`__
for obtaining the virtual appliance if you have not done so already.

**Use the Clipboard sharing.** Please note that you can leverage the clipboard for sharing text between the
virtual image environment and your host system. To copy in the VM **terminal** use
``ctrl`` + ``shift`` + ``C`` and to paste use ``ctrl`` + ``shift`` + ``V``. To copy
in the VM **Firefox browser** use ``ctrl`` + ``C``.
Use the corresponding commands on your host system (e.g. ``Command`` + ``C`` on MacOS).

.. include:: _writing_intro.rst
.. include:: _writing_test_and_serve_appliance.rst
.. include:: _writing_parameters.rst
.. include:: _writing_publish_intro.rst
.. include:: _writing_scripts.rst
.. include:: _writing_suites.rst
.. include:: _writing_conclusion.rst

.. _Galaxy: http://galaxyproject.org/
.. _GitHub: https://github.com/
.. _Docker: https://www.docker.com/
.. _Homebrew: http://brew.sh/
.. _linuxbrew: https://github.com/Homebrew/linuxbrew
.. _Vagrant: https://www.vagrantup.com/
.. _OVA: https://en.wikipedia.org/wiki/Open_Virtualization_Format
