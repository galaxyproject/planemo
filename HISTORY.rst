.. :changelog:

History
-------

.. to_doc

---------------------
0.9.0 (2015-05-03)
---------------------

* Add new logo to the README thanks to @petrkadlec from `puradesign.cz
  <http://puradesign.cz/en>`__ and @carlfeberhard from the Galaxy Project.
  `Issue 108`_
* Implement smarter ``shed_diff`` command - it now produces a meaningful
  exit codes and doesn't report differences if these correspond to attributes
  that will be automatically populated by the Tool Shed. `Issue 167`_
* Use new smarter ``shed_diff`` code to implement a new ``--check_diff``
  option for ``shed_upload`` - to check for meaningful differences before
  updating repositories. `Issue 168`_
* Record git commit hash during ``shed_upload`` if the ``.shed.yml`` is
  located in a git repository. `Issue 170`_
* Allow ``shed_`` operations to operate on git URLs directly. `Issue 169`_
* Fail if missing file inclusion statements encountered during ``.shed.yml``
  repository resolution - bug reported by @peterjc. `Issue 158`_
* Improved exception handling for tool shed operations including new 
  ``--fail_fast`` command-line option. * `Issue 114`_, `Pull Request 173`_
* Implement more validation when using the ``shed_init`` command. 1cd0e2d_
* Add ``-r/--recursive`` option to ``shed_download`` and ``shed_diff`` 
  commands and allow these commands to work with ``.shed.yml`` files defining
  multipe repositories. 40a1f57_
* Add ``--port`` option to the ``serve`` and ``tool_factory`` commands.
  15804be_
* Fix problem introduced with ``setup.py`` during the 0.9.0 development cycle
  - thanks to @peterjc. `Pull Request 171`_
* Fix clone bug introduced during 0.9.0 development cycle - thanks to
  @bgruening. `Pull Request 175`_

---------------------
0.8.4 (2015-04-30)
---------------------

* Fix for Travis CI testing picking up invalid tests (reported by @takadonet). `Issue 161`_
* Fix tar ordering for consistency (always sort by name) - thanks to @peterjc.  `Pull Request 164`_, `Issue 159`_
* Fix exception handling related to tool shed operations - thanks to @peterjc. `Pull Request 155`_, b86fe1f_

---------------------
0.8.3 (2015-04-29)
---------------------

* Fix bug where ``shed_lint`` was not respecting the ``-r/--recursive`` flag.
  9ff0d2d_
* Fix bug where planemo was producing tar files incompatible with the Tool
  Shed for package and suite repositories. a2ee135_

---------------------
0.8.2 (2015-04-29)
---------------------

* Fix bug with ``config_init`` command thanks to @bgruening. `Pull Request 151`_
* Fix unnessecary ``lint`` warning about ``parallelism`` tag reported by
  @peterjc. 9bf1eab_

---------------------
0.8.1 (2015-04-28)
---------------------

* Fixes for the source distribution to allow installation of 0.8.0 via Homebrew.

---------------------
0.8.0 (2015-04-27)
---------------------

* Implement the new ``shed_lint`` command that verifies various aspects of tool
  shed repositories - including XSD_ validation of ``repository_dependencies.xml``
  and ``tool_dependencies.xml`` files, best practices for README files, and the
  contents of ``.shed.yml`` files. This requires the lxml_ library to be available
  to Planemo or the application xmllint_ to be on its ``PATH``. `Pull Request 130`_
  `Issue 89`_ `Issue 91`_ 912df02_ d26929e_ 36ac6d8_
* Option to enable experimental XSD_ based validation of tools when ``lint``
  is executed with the new ``--xsd`` flag. This validation occurs against the
  unofficial `Galaxy Tool XSD project <https://github.com/JeanFred/Galaxy-XSD>`__
  maintained by @JeanFred. This requires the lxml_ library to be
  available to Planemo or the application xmllint_ to be on its ``PATH``.
  `Pull Request 130`_ 912df02_
* Allow skipping specific linters when using the ``lint`` command using the new
  ``--skip`` option. 26e3cdb_
* Implement sophisticated options in ``.shed.yml`` to map a directory to many,
  custom Tool Shed repositories during shed operaitons such ``shed_upload``
  including automatically mapping tools to their own directories and automatically
  building suites repositories. `Pull Request 143`_
* Make ``shed_upload`` more intelligent when building tar files so that package
  and suite repositories may have README files in source control and they will
  just be filtered out during upload. 53edd99_
* Implement a new ``shed_init`` command that will help bootstrap ``.shed.yml``
  files in the specified directory. cc1a447_
* Extend ``shed_init`` to automatically build a ``repository_rependencies.xml``
  file corresponding to a Galaxy workflow (``.ga`` file). `Issue 118`_ 988de1d_
* In addition to a single file or directory, allow ``lint`` to be passed multiple
  files. 343902d_ `Issue 139`_
* Add ``-r/--recursive`` option to ``shed_create`` and ``lint`` commands. 63cd431_
  01f2af9_
* Improved output formatting and option to write diffs to a file for the
  ``shed_diff`` command. 965511d_
* Fix lint problem when using new Galaxy testing features such as expecting
  job failures and verifing job output. `Issue 138`_
* Fix typo in ``test`` help thanks to first time contributor @pvanheus.
  `Pull Request 129`_ 1982076_
* Fix NPE on empty ``help`` element when linting tools. `Issue 124`_
* Fix ``lint`` warnings when ``configfiles`` are defined in a tool. 1a85493_
* Fix for empty ``.shed.yml`` files. b7d9e96_
* Fix the ``test`` command for newer versions of nose_. 33294d2_
* Update help content and documentation to be clear ``normalize`` should not
  be used to update the contents of tool files at this time. 08de8de_
* Warn on unknown ``command`` attributes when linting tools (anything but
  ``interpreter``). 4f61025_
* Various design, documentation (including new documentation on Tool Shed
  `publishing <http://planemo.readthedocs.org/en/latest/publishing.html>`__),
  and testing related improvements (test coverage has risen from 65% to over
  80% during this release cycle).

---------------------
0.7.0 (2015-04-13)
---------------------

* Implement `shed_create` command to create Tool Shed repositories from
  ``.shed.yml`` files (thanks to Eric Rasche). `Pull Request 101`_
* Allow automatic creation of missing repositories  during ``shed_upload``
  with the new ``--force_repository_creation`` flag (thanks to Eric Rasche).
  `Pull Request 102`_
* Allow specifying files to exclude in ``.shed.yml`` when creating tar files
  for ``shed_upload`` (thanks to Björn Grüning). `Pull Request 99`_
* Resolve symbolic links when building Tool Shed tar files with
  ``shed_upload`` (thanks to Dave Bouvier). `Pull Request 104`_
* Add a `Contributor Code of Conduct
  <https://planemo.readthedocs.org/en/latest/conduct.html>`__.
  `Pull Request 113`_
* Omit ``tool_test_output.json`` from Tool Shed tar file created with
  ``shed_upload`` (thanks to Dave Bouvier). `Pull Request 111`_
* Update required version of bioblend_ to ``0.5.3``. Fixed `Issue 88`_.
* Initial work on implementing tests cases for Tool Shed functionality.
  182fe57_
* Fix incorrect link in HTML test report (thanks to Martin Čech). 4c71299_
* Download Galaxy from the new, official Github repository. 7c69bf6_
* Update travis_test to install stable planemo from PyPI. 39fedd2_
* Enable caching on ``--install_galaxy`` by default (disable with
  ``--no_cache_galaxy``). d755fe7_

---------------------
0.6.0 (2015-03-16)
---------------------

* Many enhancements to the tool building documentation - descriptions of macros, collections, simple and conditional parameters, etc...
* Fix ``tool_init`` to quote file names (thanks to Peter Cock).  `Pull Request 98`_.
* Allow ignoring file patterns in ``.shed.yml`` (thanks to Björn Grüning). `Pull Request 99`_
* Add ``--macros`` flag to ``tool_init`` command to generate a macro file as part of tool generation. ec6e30f_
* Add linting of tag order for tool XML files. 4823c5e_
* Add linting of ``stdio`` tags in tool XML files. 8207026_
* More tests, much higher test coverage. 0bd4ff0_

---------------------
0.5.0 (2015-02-22)
---------------------

* Implement ``--version`` option. `Issue 78`_
* Implement ``--no_cleanup`` option for ``test`` and ``serve`` commands to
  persist temp files. 2e41e0a_
* Fix bug that left temp files undeleted. `Issue 80`_
* More improvements to release process. fba3874_

---------------------
0.4.2 (2015-02-21)
---------------------

* Fix setup.py for installing non-Python data from PyPI (required newer
  for ``tool_factory`` command and reStructuredText linting). Thanks to
  Damion Dooley for the bug report. `Issue 83`_

---------------------
0.4.1 (2015-02-16)
---------------------

* Fix README.rst so it renders properly on PyPI.

---------------------
0.4.0 (2015-02-16)
---------------------

* Implement ``tool_init`` command for bootstrapping creation of new
  tools (with `tutorial <http://planemo.readthedocs.org/en/latest/writing.html>`_.) 78f8274_
* Implement ``normalize`` command for reorganizing tool XML and macro
  debugging. e8c1d45_
* Implement ``tool_factory`` command to spin up Galaxy pre-configured the
  `Tool Factory
  <http://bioinformatics.oxfordjournals.org/content/early/2012/09/27/bioinformatics.bts573.full.pdf>`_. 9e746b4_
* Added basic linting of ``command`` blocks. b8d90ab_
* Improved linting of ``help`` blocks, including verifying valid
  `reStructuredText`. 411a8da_
* Fix bug related to ``serve`` command not killing Galaxy properly when complete. 53a6766_
* Have ``serve`` command display tools at the top level instead of in shallow sections. badc25f_
* Add additional dependencies to ``setup.py`` more functionality works out
  of the box. 85b9614_
* Fix terrible error message related to ``bioblend`` being unavailable.
  `Issue 70`_
* Various smaller documentation and project structure improvements.

---------------------
0.3.1 (2015-02-15)
---------------------

* Fixes to get PyPI workflow working properly.

---------------------
0.3.0 (2015-02-13)
---------------------

* Add option (``-r``) to the ``shed_upload`` command to recursively upload
  subdirectories (thanks to Eric Rasche). `Pull Request 68`_
* Fix diff formatting in test reports (thanks to Eric Rasche).
  `Pull Request 63`_
* Grab updated test database to speed up testing (thanks to approach from
  Eric Rasche and Dannon Baker). `Issue 61`_, dff4f33_
* Fix test data command-line argument name (was ``test-data`` now it is
  ``test_data``). 834bfb2_
* Use ``tool_data_table_conf.xml.sample`` file if
  ``tool_data_table_conf.xml.test`` is unavailable. Should allow some
  new tools to be tested without modifying Galaxy's global
  ``tool_data_table_conf.xml`` file. ac4f828_

---------------------
0.2.0 (2015-01-13)
---------------------

* Improvements to way Planemo loads its own copy of Galaxy modules to prevent
  various conflicts when launching Galaxy from Planemo. `Pull Request 56`_
* Allow setting various test output options in ``~/.planemo.yml`` and disabling
  JSON output. 21bb463_
* More experimental Brew and Tool Shed options that should not be considered
  part of Planemo's stable API. See bit.ly/gxbrew1 for more details.
* Fix ``project_init`` for BSD tar (thanks to Nitesh Turaga for the bug
  report.) a4110a8_
* Documentation fixes for tool linting command (thanks to Nicola Soranzo).
  `Pull Request 51`_

---------------------
0.1.0 (2014-12-16)
---------------------

* Moved repository URL to https://github.com/galaxyproject/planemo.
* Support for publishing to the Tool Shed. `Pull Request 6`_
* Support for producing diffs (``shed_diff``) between local repositories and
  the Tool Shed (based on scripts by Peter Cock). `Pull Request 33`_
* Use tool's local test data when available - add option for configuring
  ``test-data`` target. `Pull Request 1`_
* Support for testing tool features dependent on cached data. 44de95c_
* Support for generating XUnit tool test reports. 82e8b1f_
* Prettier HTML reports for tool tests. 05cc9f4_
* Implement ``share_test`` command for embedding test result links in pull
  requests. `Pull Request 40`_
* Fix for properly resolving links during Tool Shed publishing (thanks to Dave
  Bouvier). `Pull Request 29`_
* Fix for citation linter (thanks to Michael Crusoe for the bug report). af39061_
* Fix tool scanning for tool files with fewer than 10 lines (thanks to Dan
  Blankenberg). a2c13e4_
* Automate more of Travis CI testing so the scripts added to tool repository
  can be smaller. 20a8680_
* Documentation fixes for Travis CI (thanks to Peter Cock). `Pull Request 22`_,
  `Pull Request 23`_
* Various documentation fixes (thanks to Martin Čech). 36f7cb11_, b9232e55_
* Various smaller fixes for Docker support, tool linting, and documentation.

---------------------
0.0.1 (2014-10-04)
---------------------

* Initial work on the project - commands for testing, linting, serving Galaxy
  tools - and more experimental features involving Docker and Homebrew. 7d07782_

.. github_links
.. _Issue 114: https://github.com/galaxyproject/planemo/issues/114
.. _Pull Request 173: https://github.com/galaxyproject/planemo/pull/173
.. _Issue 108: https://github.com/galaxyproject/planemo/issues/108
.. _15804be: https://github.com/galaxyproject/planemo/commit/15804be
.. _Issue 158: https://github.com/galaxyproject/planemo/issues/158
.. _Pull Request 171: https://github.com/galaxyproject/planemo/pull/171
.. _1cd0e2d: https://github.com/galaxyproject/planemo/commit/1cd0e2d
.. _40a1f57: https://github.com/galaxyproject/planemo/commit/40a1f57
.. _Pull Request 175: https://github.com/galaxyproject/planemo/pull/175
.. _Issue 167: https://github.com/galaxyproject/planemo/issues/167
.. _Issue 170: https://github.com/galaxyproject/planemo/issues/170
.. _Issue 169: https://github.com/galaxyproject/planemo/issues/169
.. _Issue 168: https://github.com/galaxyproject/planemo/issues/168
.. _b86fe1f: https://github.com/galaxyproject/planemo/commit/b86fe1f
.. _Pull Request 155: https://github.com/galaxyproject/planemo/pull/155
.. _Pull Request 164: https://github.com/galaxyproject/planemo/pull/164
.. _Issue 159: https://github.com/galaxyproject/planemo/issues/159
.. _Issue 161: https://github.com/galaxyproject/planemo/issues/161
.. _a2ee135: https://github.com/galaxyproject/planemo/commit/a2ee135
.. _9ff0d2d: https://github.com/galaxyproject/planemo/commit/9ff0d2d
.. _Pull Request 151: https://github.com/galaxyproject/planemo/pull/151
.. _9bf1eab: https://github.com/galaxyproject/planemo/commit/9bf1eab
.. _Pull Request 143: https://github.com/galaxyproject/planemo/pull/143
.. _Issue 139: https://github.com/galaxyproject/planemo/issues/139
.. _Issue 89: https://github.com/galaxyproject/planemo/issues/#89
.. _Issue 91: https://github.com/galaxyproject/planemo/issues/#91
.. _d26929e: https://github.com/galaxyproject/planemo/commit/d26929e
.. _36ac6d8: https://github.com/galaxyproject/planemo/commit/36ac6d8
.. _08de8de: https://github.com/galaxyproject/planemo/commit/08de8de
.. _4f61025: https://github.com/galaxyproject/planemo/commit/4f61025
.. _1982076: https://github.com/galaxyproject/planemo/commit/1982076
.. _Pull Request 129: https://github.com/galaxyproject/planemo/pull/129
.. _912df02: https://github.com/galaxyproject/planemo/commit/912df02
.. _Pull Request 130: https://github.com/galaxyproject/planemo/pull/130
.. _1a85493: https://github.com/galaxyproject/planemo/commit/1a85493
.. _53edd99: https://github.com/galaxyproject/planemo/commit/53edd99
.. _988de1d: https://github.com/galaxyproject/planemo/commit/988de1d
.. _Issue 118: https://github.com/galaxyproject/planemo/issues/118
.. _cc1a447: https://github.com/galaxyproject/planemo/commit/cc1a447
.. _b7d9e96: https://github.com/galaxyproject/planemo/commit/b7d9e96
.. _Issue 138: https://github.com/galaxyproject/planemo/issues/#138
.. _Issue 124: https://github.com/galaxyproject/planemo/issues/#124
.. _26e3cdb: https://github.com/galaxyproject/planemo/commit/26e3cdb
.. _63cd431: https://github.com/galaxyproject/planemo/commit/63cd431
.. _965511d: https://github.com/galaxyproject/planemo/commit/965511d
.. _01f2af9: https://github.com/galaxyproject/planemo/commit/01f2af9
.. _343902d: https://github.com/galaxyproject/planemo/commit/343902d
.. _33294d2: https://github.com/galaxyproject/planemo/commit/33294d2
.. _4c71299: https://github.com/galaxyproject/planemo/commit/4c71299
.. _Pull Request 111: https://github.com/galaxyproject/planemo/pull/111
.. _Pull Request 99: https://github.com/galaxyproject/planemo/pull/99
.. _Pull Request 101: https://github.com/galaxyproject/planemo/pull/101
.. _Pull Request 102: https://github.com/galaxyproject/planemo/pull/102
.. _Issue 88: https://github.com/galaxyproject/planemo/issues/88
.. _182fe57: https://github.com/galaxyproject/planemo/commit/182fe57
.. _Pull Request 104: https://github.com/galaxyproject/planemo/pull/104
.. _7c69bf6: https://github.com/galaxyproject/planemo/commit/7c69bf6
.. _39fedd2: https://github.com/galaxyproject/planemo/commit/39fedd2
.. _d755fe7: https://github.com/galaxyproject/planemo/commit/d755fe7
.. _Pull Request 113: https://github.com/galaxyproject/planemo/pull/113
.. _Pull Request 98: https://github.com/galaxyproject/planemo/pull/98
.. _0bd4ff0: https://github.com/galaxyproject/planemo/commit/0bd4ff0
.. _Pull Request 99: https://github.com/galaxyproject/planemo/pull/99
.. _ec6e30f: https://github.com/galaxyproject/planemo/commit/ec6e30f
.. _8207026: https://github.com/galaxyproject/planemo/commit/8207026
.. _4823c5e: https://github.com/galaxyproject/planemo/commit/4823c5e
.. _2e41e0a: https://github.com/galaxyproject/planemo/commit/2e41e0a
.. _fba3874: https://github.com/galaxyproject/planemo/commit/fba3874
.. _Issue 78: https://github.com/galaxyproject/planemo/issues/78
.. _Issue 80: https://github.com/galaxyproject/planemo/issues/80


.. _Issue 83: https://github.com/galaxyproject/planemo/issues/83
.. _Issue 70: https://github.com/galaxyproject/planemo/issues/70
.. _Pull Request 68: https://github.com/galaxyproject/planemo/pull/68
.. _Issue 61: https://github.com/galaxyproject/planemo/issues/61
.. _Pull Request 63: https://github.com/galaxyproject/planemo/pull/63
.. _Pull Request 56: https://github.com/galaxyproject/planemo/pull/56
.. _Pull Request 51: https://github.com/galaxyproject/planemo/pull/51
.. _Pull Request 40: https://github.com/galaxyproject/planemo/pull/40
.. _Pull Request 29: https://github.com/galaxyproject/planemo/pull/29
.. _Pull Request 22: https://github.com/galaxyproject/planemo/pull/22
.. _Pull Request 23: https://github.com/galaxyproject/planemo/pull/23
.. _Pull Request 33: https://github.com/galaxyproject/planemo/pull/33
.. _Pull Request 6: https://github.com/galaxyproject/planemo/pull/6
.. _Pull Request 1: https://github.com/galaxyproject/planemo/pull/1

.. _3499ca0: https://github.com/galaxyproject/planemo/commit/3499ca0a15affcaf8ac9efc55880da40b0626679
.. _85b9614: https://github.com/galaxyproject/planemo/commit/85b961465f46351507f80ddc3758349535060502
.. _53a6766: https://github.com/galaxyproject/planemo/commit/53a6766cdebdddc976189f6dc6a264bb4105c4bf
.. _badc25f: https://github.com/galaxyproject/planemo/commit/badc25fca495b61457ffb2e027f3fe9cf17c798f
.. _411a8da: https://github.com/galaxyproject/planemo/commit/411a8da21c92ba37c7ad95bfce9928d9b8fd998e
.. _b8d90ab: https://github.com/galaxyproject/planemo/commit/b8d90abab8bf53ae2e7cca4317223c01af9ab68c
.. _e8c1d45: https://github.com/galaxyproject/planemo/commit/e8c1d45f0c9a11bcf69ec2967836c3b8f432dd97
.. _78f8274: https://github.com/galaxyproject/planemo/commit/78f82747996e4a28f96c85ad72efe5e54c8c74bd
.. _9e746b4: https://github.com/galaxyproject/planemo/commit/9e746b455e3b15219878cddcdeda722979639401
.. _ac4f828: https://github.com/galaxyproject/planemo/commit/ac4f82898f7006799142503a33c3978428660ce7
.. _834bfb2: https://github.com/galaxyproject/planemo/commit/834bfb2929d367892a3abe9c0b88d5a0277d7905
.. _dff4f33: https://github.com/galaxyproject/planemo/commit/dff4f33c750a8dbe651c38e149a26dd42e706a82
.. _a4110a8: https://github.com/galaxyproject/planemo/commit/a4110a85a770988e5cd3c31ccc9475717897d59c
.. _21bb463: https://github.com/galaxyproject/planemo/commit/21bb463ad6c321bcb669603049a5e89a69766ad9
.. _af39061: https://github.com/galaxyproject/planemo/commit/af390612004dab636d8696839bb723d39f97c85d
.. _20a8680: https://github.com/galaxyproject/planemo/commit/20a86807cb7ea87db2dbc0197ae08a40df3ab2bc
.. _44de95c: https://github.com/galaxyproject/planemo/commit/44de95c0d7087a5822941959f9a062f6382e329b
.. _82e8b1f: https://github.com/galaxyproject/planemo/commit/82e8b1f17eae526aeb341cb4fffb8d09d73bb419
.. _05cc9f4: https://github.com/galaxyproject/planemo/commit/05cc9f485ee87bc344e3f43bb1cfd025a16a6247
.. _32c6e7f: https://github.com/galaxyproject/planemo/commit/32c6e7f78bb8f04d27615cfd8948b0b89f27b4e6
.. _7d07782: https://github.com/galaxyproject/planemo/commit/7d077828559c9c9c352ac814f9e3b86b1b3a2a9f
.. _a2c13e4: https://github.com/galaxyproject/planemo/commit/a2c13e46259e3be35de1ecaae858ba818bb94734
.. _36f7cb11: https://github.com/galaxyproject/planemo/commit/36f7cb114f77731f90860d513a930e10ce5c1ba5
.. _b9232e55: https://github.com/galaxyproject/planemo/commit/b9232e55e713abbd1d9ce8b0b34cbec6c701dc17

.. _bioblend: https://github.com/galaxyproject/bioblend/
.. _XSD: http://www.w3schools.com/schema/
.. _lxml: http://lxml.de/
.. _xmllint: http://xmlsoft.org/xmllint.html
.. _nose: https://nose.readthedocs.org/en/latest/
