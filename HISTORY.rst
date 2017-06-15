.. :changelog:

History
-------

.. to_doc

---------------------
0.42.0 (2017-06-15)
---------------------

* Conda/Container documentation and option naming improvements. `Pull Request
  684`_
* Sync `galaxy.xsd` with latest upstream Galaxy updates (thanks to `@nsoranzo`_).
  `Pull Request 687`_
* Fix `ci_find_repos` command to not filter repos whose only modifications where
  in subdirs (thanks to `@nsoranzo`_).
  `Pull Request 688`_
* Update `container_register` for mulled version 2 and repository name changes.
  `Pull Request 689`_
* Better pull request messages for the `container_register` command.
  `Pull Request 690`_

---------------------
0.41.0 (2017-06-05)
---------------------

* Fix ``shed_update`` not fail if there is nothing to update
  (thanks to `@nsoranzo`_). `Issue 494`_, `Pull Request 680`_
* Conda documentation and option naming improvements.
  `Pull Request 683`_
* Implement ``container_register`` for tool repositories.
  `Pull Request 675`_
* Fix ``hub`` binary installation for Mac OS X.
  `Pull Request 682`_

---------------------
0.40.1 (2017-05-03)
---------------------

* Fix data manager configuration to not conflict with original Galaxy at
  ``galaxy_root`` (thanks to `@nsoranzo`_). `Pull Request 662`_
* Fix ``filter_paths()`` to not partial match paths when filtering shed repositories
  (thanks to `@nsoranzo`_). `Pull Request 665`_
* Fix description when creating ``.shed.yml`` files (thanks to `@RJMW`_).
  `Pull Request 664`_

---------------------
0.40.0 (2017-03-16)
---------------------

* Implement instructions and project template for GA4GH Tool Execution
  Challenge Phase 1. 84c4a73_
* Eliminate Conda hack forcing ``/tmp`` as temp directory. b4ae44d_
* Run dependency script tests in isolated directories. 32f41c9_
* Fix OS X bug in ``planemo run`` by reworking it to wait using urllib instead of sockets.
  3129216_

---------------------
0.39.0 (2017-03-15)
---------------------

* Implement documentation and examples for Conda-based dependency development (under
  "Advanced" topics).
  `Pull Request 642`_, `Pull Request 643`_
* Implement documentation and examples for container-based dependency development (under
  "Advanced" topics).
  0a1abfe_
* Implement a ``planemo conda_search`` command for searching best practice channels
  from the command line.
  `Pull Request 642`_
* Allow Planemo to work with locally built Conda packages using the ``--conda_use_local``
  command.
  `Pull Request 643`_, `Issue 620`_
* Implement an ``open`` (or just ``o``) command to quickly open the last test results
  (or any file if supplied). `Pull Request 641`_
* Linting improvements and fixes due to `galaxy-lib`_ update.
  * WARN on test output names not found or not matching.
  * INFO correct information about stdio if profile version is found.
  * WARN if profile version is incorrect.
  * INFO profile version
  * Fix ``assert_command`` not detected as a valid test (fixes  `Issue 260`_).
* Have ``lint --conda_requirements`` check that at least one actual requirement is found.
  6638caa_
* Allow ``conda_install`` to work with packages as well as just tools.
  8faf661_
* Add ``--global`` option to conda_install to install requirements into global Conda setup
  instead of using an environment.
  8faf661_
* Implement ``planemo lint --biocontainer`` that checks that a tool has an available BioContainer
  registered.
  0a1abfe_
* Add more options and more documentation to the ``planemo mull`` command.
  0a1abfe_
* Hack around a bug in Conda 4.2 that makes it so ``planemo mull`` doesn't work out of the box on
  Mac OS X.
  0a1abfe_
* Allow URIs to be used instead of paths for a couple operations. ce0dc4e_
* Implement non-strict CWL parsing option. 4c0f100_
* Fixes for changes to cwltool_ and general CWL-relate functionality.
  3c95b7b_, 06bcf19_, 525de8f_, 9867e56_, 9ab4a0d_
* Eliminate deprecated XML-based abstraction from ``planemo.tools``. 04238d3_
* Fix ``MANIFEST.in`` entry that was migrated to galaxy-lib. ced5ce2_
* Various fixes for the command ``conda_env``. `Pull Request 640`_
* Improved command help - both formatting and content. `Pull Request 639`_
* Implement a ``--no_dependency_resolution`` option disabling conda dependency
  resolver.
  `Pull Request 635`_, `Issue 633`_
* Tests for new linting logic. `Pull Request 638`_
* Fix bug where tool IDs needs to be lowercase for the shed (thanks to
  `@bgruening`_).
  `Pull Request 649`_
* Update seqtk version targetted by intro docs. e343b67_
* Various other Conda usability improvements. `Pull Request 634`_

---------------------
0.38.1 (2017-02-06)
---------------------

* Fix bug with ``shed_lint --urls`` introduced in 0.38.0.
  84ebc1f_

---------------------
0.38.0 (2017-02-06)
---------------------

* Trim down the default amount of logging during testing.
  `Pull Request 629`_, `Issue 515`_
* Improved log messages during shed operations. 08c067c_
* Update tool XSD against latest Galaxy.
  fca4183_, 03c9658_
* Fix bug where ``shed_lint --tools`` for a suite lints the same tools multiple
  times.
  `Issue 564`_, `Pull Request 628`_

---------------------
0.37.0 (2017-01-25)
---------------------

* Update to the latest `galaxy-lib`_ release. This means new installs start with
  Miniconda 3 instead of Minicoda 2 and at a newer version. This fixes many
  Conda_ related bugs.
* Change defaults so that Conda automatically initializes and performs tool installs
  by default from within the spawned Galaxy server. The trio of flags
  ``--conda_dependency_resolution``, ``--conda_auto_install``, and ``--conda_auto_init``
  are effectively enabled by default now. 4595953_
* Use the Galaxy cached dependency manager by default (thanks to `@abretaud`_).
  `Pull Request 612`_
* Test Conda dependency resolution for more versions of Galaxy including the forthcoming
  release of 17.01.
* Update to the latest Galaxy tool XSD for various tool linting fixes. 32acd68_
* Fix pip ignores for ``bioconda_scripts`` (thanks to `@nturaga`_)
  `Pull Request 614`_

---------------------
0.36.1 (2016-12-12)
---------------------

* Fix move error when using ``project_init``.
  `Issue 388`_, `Pull Request 610`_
* Improved integration testing for ``test`` command. `Pull Request 609`_
* Update CWL links to v1.0 (thanks to `@mr-c`_).
  `Pull Request 608`_

---------------------
0.36.0 (2016-12-11)
---------------------

* Bring in latest tool XSD file from Galaxy (thanks to `@peterjc`_).
  `Pull Request 605`_
* PEP8 fixes for various linting problems 
  (thanks to `@peterjc`_).
  `Pull Request 606`_
* Update tool syntax URL to new URL (thanks to `@mvdbeek`_).
  `Pull Request 602`_

---------------------
0.35.0 (2016-11-14)
---------------------

* Native support for building bioconductor tools and recipes
  (thanks to `@nturaga`_). `Pull Request 570`_
* Fixes for running Galaxy via docker-galaxy-stable (thanks to
  `@bgruening`_). 50d3c4a_
* Import order linting fixes (thanks to `@bgruening`_).

---------------------
0.34.1 (2016-10-12)
---------------------

* Mimic web browser to validate user help URLs fixing `Issue 578`_
  (thanks to `@peterjc`_). `Pull Request 591`_
* Fix for Bioconda recipes depending on ``conda-forge`` (thanks to `@nsoranzo`_).
  `Pull Request 590`_


---------------------
0.34.0 (2016-10-05)
---------------------

* Implement ``mull`` command to build containers for tools based on Conda_
  recipes matching requirement definitions. 08cef54_
* Implement ``--mulled_containers`` flag on ``test``, ``serve``, and ``run``
  commands to run tools in "mulled" containers. Galaxy will first search
  locally cache containers (such as ones built with ``mull``), then search
  the mulled namespace of `quay.io`_, and finally build one on-demand if
  needed using `galaxy-lib`_ and Involucro_ developed by `@thriqon`_.
* Implement ``--conda_requirements`` flag on ``lint`` command to ensure requirements
  can be resolved in best practice channels. 9da8387_
* Allow ``conda_install`` command over multiple tool paths. 2e4e5fc_
* Update pip_ as part of setting virtual environment in ``Makefile`` target.
  19b2ee9_
* Add script to auto-update Bioconda_ recipe for Planemo and open a pull request.
  f0da66f_

---------------------
0.33.2 (2016-09-28)
---------------------

* Fix HISTORY.rst link problem that prevented correct display of content on PyPI.

---------------------
0.33.1 (2016-09-28)
---------------------

* Fix ``lint --urls`` false positives by being more restrictive with what is considered a URL
  (fixed by `@erasche`_ after detailed report from `@peterjc`_).
  `Issue 573`_, `Pull Request 579`_

---------------------
0.33.0 (2016-09-23)
---------------------

* Enable XSD validation of tools by default (restore old behavior with
  ``planemo lint --no_xsd``). 1ef05d2_
* Implement a ``conda_lint`` command to lint Conda_ recipes based
  on `anaconda-verify`_. 6a6f164_
* Implement ``clone`` and ``pull_request`` commands to ease PRs
  (with documentation fixes from `@martenson`_).
  e925ba1_, ea5324f_
* Update `galaxy.xsd`_ to allow version_command's to have an interpreter
  attribute. 7cca2e4_
* Apply improvement from `@nsoranzo`_ for Planemo's use of git_ 
  `diff <https://git-scm.com/docs/git-diff>`__.
  6f91719_
* Pull in downstream refactoring of ``tool_init`` code from `@nturaga`_'s 
  Bioconductor_ work. ccdd2d5_
* Update to latest `Tool Factory`_ code from `tools-iuc`_. ca88b0c_
* Small code cleanups. b6d8294_, d6da3a8_
* Fixup docs in ``planemo.xml.validation``.
* Allow skipping newly required lxml_ dependency in `setup.py`_. 34538de_
    

---------------------
0.32.0 (2016-09-16)
---------------------

* Enhance ``planemo lint --xsd`` to use a fairly complete and newly official XSD
  definition. `Pull Request 566`_
* Migrate and update documentation related to tool XML macros and handling 
  multiple outputs from the Galaxy wiki (with help from `@bgruening`_, `@mvdbeek`_,
  and `@nsoranzo`_). `Pull Request 559`_
* Documentation fixes (thanks to `@ramezrawas`_). `Pull Request 561`_
* Do not fail URL linting in case of too many requests (thanks to `@nsoranzo`_).
  `Pull Request 565`_

---------------------
0.31.0 (2016-09-06)
---------------------

* Implement new commands to ``ci_find_repos`` and ``ci_find_tools`` to ease
  CI scripting.
  `Pull Request 555`_
    

---------------------
0.30.2 (2016-09-01)
---------------------

* Fix another problem with Conda_ prefix handling when using
  ``--conda_dependency_resolution``. f7b6c7e_

---------------------
0.30.1 (2016-09-01)
---------------------

* Fix a problem with Conda_ prefix handling when using
  ``--conda_dependency_resolution``. f7b6c7e_
* Fix for quote problem in ``update_planemo_recipe.bash``. 6c03de8_
* Fix to restore linting of ``tests/`` directory and fix import order 
  throughout module. ef4b9f4_

---------------------
0.30.0 (2016-09-01)
---------------------

* Update to the latest `galaxy-lib`_ release and change Conda_ semantics to match
  recent updates to Galaxy. For the most robust Conda_ usage - use planemo 0.30+
  with Galaxy 16.07 or master.
  07d94bd_
* Implement the ``--conda_auto_init`` flag for ``conda_install``. ca19910_
* Allow the environment variable ``PLANEMO_CONDA_PREFIX`` to set a default
  for ``--conda_prefix``.
  24008ab_
* Fixup documentation regarding installs and Conda_. ce44e87_
* Fix and lint Python module import order throughout project.
  `Pull Request 550`_
* Use ``cp`` rather than symlink to ``$DOWNLOAD_CACHE`` in the
  ``dependency_script`` command (thanks to `@peterjc`_).  c2204b3_
* Fixes for the Homebrew recipe updater. c262b6d_

---------------------
0.29.1 (2016-08-19)
---------------------

* Improved handling of Python 2.7 specific dependencies.

---------------------
0.29.0 (2016-08-19)
---------------------

* Look for sha256sum checksums during shed_lint (thanks to `@peterjc`_).
  `Pull Request 539`_
* An assortment fixes and enhancements to the ``dependency_script`` command
  (thanks to `@peterjc`_). `Pull Request 541`_, `Pull Request 545`_
* Fix shed_build to respect exclude: in .shed.yml (thanks to `@nsoranzo`_).
  `Pull Request 540`_
* Fix linting of tool URLs (thanks to `@nsoranzo`_). `Pull Request 546`_ 

---------------------
0.28.0 (2016-08-17)
---------------------

* Fixes for bioblend_ v0.8.0 (thanks to `@nsoranzo`_). 9fdf490_ 
* Enable shed repo type update (thanks to `@nsoranzo`_). 3ceaa40_
* Create suite repositories with repository_suite_definition type by default
  (thanks to `@nsoranzo`_).
  057f4f0_
* Include ``shed_lint`` in script run by ``travis_init`` (thanks to `@peterjc`_).
  `Pull Request 528`_
* Minor polish to the ``travis_init`` command (thanks to `@peterjc`_).
  `Pull Request 512`_
* Update pip_ and setuptools on TravisCI; fix travis_init (thanks to `@peterjc`_).
  `Pull Request 521`_
* Shorten command one line descriptions for main help (thanks to `@peterjc`_).
  `Pull Request 510`_
* Use ``planemo test --no_cache_galaxy`` under TravisCI (thanks to `@peterjc`_).
  `Pull Request 513`_
* Improve and fix docs ahead of GCC 2016 (thanks to `@martenson`_).
  `Pull Request 498`_, 725b232_
* Add description of ``expect_num_outputs`` to planemo FAQ. a066afb_
* Revise planemo tools docs to be more explicit about collection identifiers.
  a811e65_
* Add more docs on existing dynamic tool output features. `Pull Request 526`_
* Fix serve command doc (thanks to `@nsoranzo`_). 8c088c6_
* Fix `make lint-readme` (RST link errors) (thanks to `@peterjc`_).
  `Pull Request 525`_
* Add union bedgraph example to project templates (for GCC demo example). 
  d53bcd6_
* Add Flow Cytometry Analysis, Data Export, and Constructive Solid Geometry as
  shed categories (thanks to `@bgruening`_, `@gregvonkuster`_, and `@nsoranzo`_).
  e890ab5_, 08bb354_, e2398fb_
* Remove duplicated attribute in docs/writing/bwa-mem_v5.xml (thanks to
  Paul Stewart `@pstew`_).
  `Pull Request 507`_

---------------------
0.27.0 (2016-06-22)
---------------------

* Use ephemeris to handle syncing shed tools for workflow actions.
  1c6cfbb_
* More planemo testing enhancements for testing artifacts that aren't
  Galaxy tools. `Pull Request 491`_
* Implement ``docker_galaxy`` engine type. eb039c0_, `Issue 15`_
* Enhance profiles to be Dockerized Galaxy-aware. `Pull Request 488`_
* Add linter for DOI type citation - thanks to `@mvdbeek`_.
  `Pull Request 484`_

---------------------
0.26.0 (2016-05-20)
---------------------

* Implement ``Engine`` and ``Runnable`` abstractions - Planemo now has
  beta support for testing Galaxy workflows and CWL_ tools with Galaxy and
  any CWL artifact with cwltool_.
  `Pull Request 454`_, 7be1bf5_
* Fix missing command_line in test output json. e38c436_
* More explicit Galaxy ``job_conf.xml`` handling, fixes bugs caused by
  ``galaxy_root`` having existing and incompatible ``job_conf.xml`` files
  and makes it possible to specify defaults with fixed server name. c4dfd55_
* Introduce profile commands (``profile_create``, ``profile_delete``, and
  ``profile_list``) and profile improvements (automatic postgres database
  creation support). `Pull Request 480`_, a87899b_
* Rework Galaxy test reporting to use structured data instead of XUnit
  data. 4d29bf1_
* Refactor Galaxy configuration toward support for running Galaxy in
  docker-galaxy-stable. `Pull Request 479`_    

---------------------
0.25.1 (2016-05-11)
---------------------

* Tweak dependencies to try to fix cwltool_ related issues - such
  as `Issue 475`_.

---------------------
0.25.0 (2016-05-11)
---------------------

* Implement Galaxy "profiles" - the ability to configure 
  perisistent, named environments for ``serve`` and ``test``.
  5d08b67_
* Greatly improved ``serve`` command - make ``test-data``
  available as an FTP folder, (on 16.07) automatically log
  in an admin user, and many more options (such as those 
  required for "profiles" and a ``--daemon`` mode).
* Two fixes to ensure more consistent, dependable ``test`` output.
  `Pull Request 472`_, f3c6917_
* Add code and documentation for linting (``lint``) and
  building (``tool_init``) CWL_ tools. a4e6958_, b0b867e_,
  4cd571c_
* If needed for Conda_ workaround, shorten ``config_directory`` 
  path (thanks to `@mvdbeek`_). efc5f30_
* Fix ``--no_cache_galaxy`` option (thanks to Gildas Le 
  Corguillé). d8f2038_
* Target draft 3 of CWL_ instead of draft 2. 775bf49_
* Fix ``cwltool`` dependency version - upstream changes broke
  compatibility. `65b999d`_
* Add documentation section and slides about recent Galaxy
  tool framework changes (with fix from `@remimarenco`_). 069e7ba_
* Add IUC standards to Planemo docs. 2ae2b49_
* Improve collection-related contents in documentation
  (thanks in part to `@martenson`_).
  fea51fc_, 13a5ae7_
* Add documentation on ``GALAXY_SLOTS`` and running planemo
  on a cluster. 45135ff_, e0acf91_
* Revise command-line handling framework for consistency and
  extension - allow extra options to be configured as 
  defaults ``~/.planemo.yml`` including ``--job_config_file``
  and Conda_ configuration options. e769118_, 26e378e_
* Fix ``tool_init`` commans options typos (thanks to
  Nitesh Turaga). 826d371_
* Refactor galaxy-related modules into submodules of a new
  ``planemo.galaxy`` package. 8e96864_
* Fix error message typo (thanks to `@blankenberg`_). b1c8f1d_
* Update documentation for recent command additions. 3f4ab44_
* Rename option ``--galaxy_sqlite_database`` option to
  ``--galaxy_database_seed`` and fix it so it actually works. 
  f7554d1_
* Add ``--extra_tools`` option to ``serve`` command. 02a08a0_
* Update project testing to include linting documentation
  (``docs/``), Python import order, and docstrings.
  a13a120_, 6e1e726_, 95d5cba_


---------------------
0.24.2 (2016-04-25)
---------------------

* Revert "check ``.shed.yml`` owner against credentials during shed
  creation", test was incorrect and preventing uploads.
  `Pull Request 425`_, `Issue 246`_

---------------------
0.24.1 (2016-04-08)
---------------------

* Fix test summary report. `Pull Request 429`_
* Improve error reporting when running ``shed_test``. ce8e1be_
* Improved code comments and tests for shed related functionality.
  89674cb_
* Rev `galaxy-lib`_ dependency to 16.4.1 to fix wget usage in
  newer versions of wget. d76b489_

---------------------
0.24.0 (2016-03-29)
---------------------
    
* Drop support for Python 2.6. 93b7bda_
* A variety of fixes for ``shed_update``.
  `Pull Request 428`_, `Issue 416`_
* Fix reporting of metadata updates for invalid shed updates.
  `Pull Request 426`_, `Issue 420`_
* Check ``.shed.yml`` owner against credentials during shed creation.
  `Pull Request 425`_, `Issue 246`_
* Fix logic error if there is a problem with ``shed_create``. 358a42c_
* Tool documentation improvements. 0298510_, a58a3b8_

---------------------
0.23.0 (2016-02-15)
---------------------

* Fix duplicated attributes with Conda_ resolver (thanks
  to Björn Grüning). `Pull Request 403`_
* Upgrade to latest version of `galaxy-lib`_ for more linting.
* Attempt to better handle conditional dependency on cwltool.

---------------------
0.22.2 (2016-01-14)
---------------------

* Fixed bug targetting forthcoming release of Galaxy 16.01.

---------------------
0.22.1 (2016-01-14)
---------------------

* Fixed problem with PyPI_ build artifacts due to submodule's not
  being initialized during previous release.

---------------------
0.22.0 (2016-01-13)
---------------------

* Add ``--skip_venv`` to support running Galaxy 16.01 inside of
  conda environments. 9f3957d_
* Implement conda support. f99f6c1_, ad3b2f0_, 5e0b6d1_
* Update LICENSE for Planemo to match Galaxy. 15d33c7_
* Depend on new `galaxy-lib`_ on PyPI_ instead of previous hacks....
  `Pull Request 394`_
* Fix egg caching against master/15.10. 6d0f502_
* Fix bug causing shed publishing of ``.svn`` directories.
  `Issue 391`_
* Bug fixes for Conda_ support thanks to `@bgruening`_. 63e456c_
* Fix document issues thanks to `@einon`_.
  `Pull Request 390`_
* Improve client for shed publishing to support newer shed backend
  being developed by `@erasche`_. `Pull Request 394`_
* Tool Shed ``repo_id`` change, `@erasche`_. `Pull Request 398`_
* Various other small changes to testing, project structure, and
  Python 3 support.

---------------------
0.21.1 (2015-11-29)
---------------------

* Fix serious regression to ``test`` command. 94097c7_
* Small fixes to release process. 4e1377c_, 94645ed_

---------------------
0.21.0 (2015-11-29)
---------------------

* If ``virtualenv`` not on ``PATH``, have Planemo create one for Galaxy.
  5b97f2e_
* Add documentation section on testing tools installed in an existing
  Galaxy instance. 1927168_
* When creating a virtualenv for Galaxy, prefer Python 2.7.
  e0577e7_
* Documentation fixes and improvements thanks to `@martenson`_.
  0f8cb10_, 01584c5_, b757791_
* Specify a minimum ``six`` version requirement. 1c7ee5b_
* Add script to test a planemo as a wheel. 6514ff5_, `Issue 184`_
* Fix empty macro loading. `Issue 362`_
* Fix an issue when you run ``shed_diff --shed_target local`` thanks
  to Gwendoline Andres and Gildas Le Corguillé at ABiMS Roscoff.
  `Pull Request 375`_
* Fix ``shed_diff`` printing to stdout if ``-o`` isn't specified.
  f3394e7_
* Small ``shed_diff`` improvements to XML diffing and XUnit reporting.
  af7448c_, 83e227a_
* More logging of ``shed_diff`` results if ``--verbose`` flagged.
  9427b47_
* Add ``test_report`` command for rebuilding reports from structured JSON.
  99ee51a_
* Fix option bug with Click 6.0 thanks to `@bgruening`_. 2a7c792_
* Improved error messages for test commands. fdce74c_
* Various fixes for Python 3. 2f66fc3_, 7572e99_, 8eda729_, 764ce01_
* Use newer travis container infrastructure for testing. 6d81a94_
* Test case fixes. 98fdc8c_, 0e4f70a_
    


---------------------
0.20.0 (2015-11-11)
---------------------

* More complete I/O capturing for XUnit. 6409449_
* Check for select parameter without options when linting tools.
  `Issue 373`_
* Add ``--cwl_engine`` argument to ``cwl_run`` command. dd94ddc_
* Fixes for select parameter linting. 8b31850_
* Fix to demultiplexing repositories after tool uploads. `Issue 361`_
* Fix to update planemo for Galaxy wheels. 25ef0d5_
* Various fixes for Python 2.6 and Python 3.
  c1713d2_, 916f610_, c444855_
    

---------------------
0.19.0 (2015-11-03)
---------------------

* Initial implementation of ``cwl_run`` command that runs a
  CWL tool and job file through Galaxy. 49c5c1e_
* Add ``--cwl`` flag to ``serve`` to experimentally serve CWL tools
  in Galaxy.
  `Pull Request 339`_
* Implement highly experimental ``cwl_script`` command to convert
  a CWL job to a bash script. 508dce7_
* Add name to all XUnit reports (thanks to `@erasche`_).
  `Pull Request 343`_
* Capture stdout and stderr for ``shed_diff`` and ``shed_update`` 
  XUnit reports. `Pull Request 344`_
* More tool linting (conditionals) thanks to `@erasche`_.
  `Pull Request 350`_
* UTF-8 fixes when handling XUnit reports. `Pull Request 345`_
* Add `Epigenetics` as Tool Shed category. `Pull Request 351`_
* Merge changes to common modules shared between Galaxy, Planemo, and Pulsar (thanks to `@natefoo`_).
  `Pull Request 356`_
* Add ``--cite_url`` to ``tool_init``. fdb1b51_
* ``tool_init`` bug fix. f854138_
* Fix `setup.py`_ for cwltool and bioblend_ changes. 1a157d4_
* Add option to specify template sqlite database locally. c23569f_
* Add example IPython notebooks to docs. c8640b6_

---------------------
0.18.1 (2015-10-22)
---------------------

* Fix issue with test reporting not being populated. 19900a6_

---------------------
0.18.0 (2015-10-20)
---------------------

* Improvements to ``docker_shell`` usability (thanks to `@kellrott`_).
  `Pull Request 334`_
* Add docker pull attempt when missing Dockerfile (thanks to `@kellrott`_).
  `Pull Request 333`_
* Fix bug inferring which files are tool files (thanks to `@erasche`_).
  `Pull Request 335`_, `Issue 313`_
* Initial work toward automating brew recipe update. 4d6f7d9_, `Issue 329`_

---------------------
0.17.0 (2015-10-19)
---------------------

* Implement basic XUnit report option for ``shed_update`` (thanks to `@martenson`_).
  `Pull Request 322`_
* Fix issues with producing test outputs. 572e754_
* Xunit reporting improvements - refactoring, times, diff output (thanks to `@erasche`_).
  `Pull Request 330`_
* Implement project governance policy and update developer code of conduct to
  match that of the Galaxy project. `Pull Request 316`_
* Update filters for account for new ``.txt`` and ``.md`` test outputs
  (thanks to `@erasche`_). `Pull Request 327`_
* Add verbose logging to galaxy test output handling problems. 5d7db92_
* Flake8 fixes (thanks to `@martenson`_). 949a36d_
* Remove uses of deprecated ``mktemp`` Python standard library function
  (thanks to `@erasche`_). `Pull Request 330`_
    

---------------------
0.16.0 (2015-10-07)
---------------------

* Adding new command ``dependency_script`` to convert Tool Shed dependencies
  into shell scripts - thanks to `@peterjc`_.
  `Pull Request 310`_, f798c7e_, `Issue 303`_
* Implement profiles in sheds section of the ``~/.planemo.yml``.
  `Pull Request 314`_

---------------------
0.15.0 (2015-10-01)
---------------------

* Template framework for reporting including new markdown and plain
  text reporting options for testing - thanks to `@erasche`_.
  `Pull Request 304`_
* XUnit style reporting for ``shed_diff`` command - thanks to
  `@erasche`_. `Pull Request 305`_
* Add new ``shed_build`` command for building repository tarballs -
  thanks to `@kellrott`_. `Pull Request 297`_
* Fix exit code handling for ``lint`` commands - thanks to `@mvdbeek`_.
  `Pull Request 292`_    
* Improved documentation for ``serve`` command - thanks to `@lparsons`_.
  `Pull Request 312`_
* Tiny backward compatible Python 3 tweaks for `Tool Factory`_ - thanks
  to `@peterjc`_. dad2d9d_
* Fixed detection of virtual environment in ``Makefile`` - thanks to
  `@lparsons`_. `Pull Request 311`_
* Updates to Galaxy XSD - thanks to `@mr-c`_. `Pull Request 309`_
* Allow reading shed key option from an environment variable.
  `Pull Request 307`_
* Allow specifying host to serve Galaxy using ``-host`` - thanks in
  part to `@chambm`_. `Pull Request 301`_
* Allow specifying defaults for ``-host`` and ``--port`` in
  ``~/.planemo.yml``. `Pull Request 301`_
* Improve ``~/.planemo.yml`` sample comments - thanks to `@martenson`_.
  `Pull Request 287`_
* Update tool shed categories - thanks to `@bgruening`_. `Pull Request 285`_
* Improved output readibility for ``diff`` command - thanks to `@martenson`_. `Pull Request 284`_

---------------------
0.14.0 (2015-08-06)
---------------------

* Allow ``-t`` as shorthand for ``--shed_target`` (thanks to Peter Cock).
  `Pull Request 278`_
* Fix ``tool_init`` command to use ``from_work_dir`` only if file in command
  (thanks to bug report and initial fix outline by Gildas Le Corguillé).
  `Pull Request 277`_
* Various documentation fixes (thanks in part to Peter Cock and Daniel
  Blankenberg). `Pull Request 256`_, `Pull Request 253`_, `Pull Request 254`_, 
  `Pull Request 255`_, `Pull Request 251`_, `Issue 272`_

---------------------
0.13.2 (2015-07-06)
---------------------

* Fix project_init for missing files. cb5b906_
* Various documentation improvements.    

---------------------
0.13.1 (2015-07-01)
---------------------

* Fix for ``shed_init`` producing non-standard type hints. `Issue 243`_,
  f0610d7_
* Fix tool linting for parameters that define an ``argument`` but not a
  ``name``. `Issue 245`_, aad1eed_
* Many doc updates including a tutorial for developing tools in a test-driven
  fashion and instructions for using the planemo appliance through Kitematic
  (with Kitematic screenshots from Eric Rasche).

---------------------
0.13.0 (2015-06-28)
---------------------

* If planemo cannot find a Galaxy root, it will now automatically fetch
  one (specifing ``--galaxy_install`` will still force a fetch).
  `Pull Request 235`_
* `Docuementation <http://planemo.readthedocs.org/en/latest/appliance.html>`__
  has been updated to reflect new and vastly improved Docker and Vagrant
  virtual appliances are now available, as well as a new VirtualBox OVA
  variant.
* Update linting for new tool XML features (including ``detect_errors``
  and output collections). `Issue 233`_, 334f2d4_
* Fix ``shed_test`` help text. `Issue 223`_
* Fix code typo (thanks to Nicola Soranzo). `Pull Request 230`_
* Improvements to algorithm used to guess if an XML file is a tool XML file.
  `Issue 231`_
* Fix configuration file handling bug. `Issue 240`_

---------------------
0.12.2 (2015-05-23)
---------------------

* Fix ``shed_test`` and ``shed_serve`` for test and local tool sheds.
  f3cafaa_

---------------------
0.12.1 (2015-05-21)
---------------------

* Fix to ensure the tab completion script is in the Python source tarball
  (required for setting up tab-completion for Homebrew). 6b4e7a6_

---------------------
0.12.0 (2015-05-21)
---------------------

* Implement a ``--failed`` flag for the ``test`` command to rerun
  previously faied tests. `Pull Request 210`_
* Implement ``shed_update`` to upload contents and update repository
  metadata. `Pull Request 216`_
* Implement ``shed_test`` and ``shed_serve`` commands to test and view
  published artifacts in the Tool Shed. `Pull Request 213`_, `Issue 176`_
* Add shell tab-completion script. 37dcc07_
* Many more commands allow specifing multiple tool and/or repository targets.
  `Issue 150`_
* Add -m as alias for --message in planemo shed_upload (thanks to
  Peter Cock). `Pull Request 200`_
* Add ``--ensure_metadata`` option to ``shed_lint`` to ensure ``.shed.yml``
  files contain many repository. `Pull Request 215`_
* More developer documentation, additional ``make`` targets including ones
  for setting up git pre-commit hooks. cc8abb6_, `Issue 209`_
* Small README improvement (thanks to Martin Čech) b53006d_
* Fixes for shed operation error handling (thanks to Martin Čech).
  `Pull Request 203`_,  `Pull Request 206`_
* Fix for "smart" ``shed_diff`` not in the repository root directory
  (thanks to Peter Cock). `Pull Request 207`_, `Issue 205`_
* Recursive ``shed_diff`` with directories not yet in Tool Shed.
  `Pull Request 208`_
* Improve error handling and reporting for problematic ``--shed_target``
  values. `Issue 217`_
* Fix typos in lint messages. `Issue 211`_


---------------------
0.11.1 (2015-05-12)
---------------------

* Fix default behavior for ``planemo lint`` to use current directory if
  explicit paths are not supplied. 1e3668a_

---------------------
0.11.0 (2015-05-12)
---------------------

* More compact syntax for defining multiple custom inclusions in ``.shed.yml``
  files - thanks to Peter Cock. `Issue 180`_, `Pull Request 185`_,
  `Pull Request 196`_
* Prevent ambigous destinations when defining custom inclusions in
  ``.shed.yml``- thanks to Peter Cock. `Pull Request 186`_
* ``lint`` now warns if tool ids contain whitespace. `Pull Request 190`_
* Handle empty tar-balls gracefully on older Python versions - thanks
  to Peter Cock. `Pull Request 187`_
* Tweak quoting in ``cp`` command - thanks to Peter Cock. 6bcf699_
* Fix regression causing testing to no longer produce "pretty" test
  results under certain circumstances. `Issue 188`_
* Fix for recursive ``shed_diff`` folder naming. `Issue 192`_
* Fix output definitions to ``tool_init`` command. `Issue 189`_

---------------------
0.10.0 (2015-05-06)
---------------------

* Extend ``shed_lint`` to check for valid actions in tool_dependencies.xml
  files. 8117e03_
* Extend ``shed_lint`` to check for required files based on repository type.
  `Issue 156`_
* Ignore common editor backup files during ``shed_upload``. `Issue 179`_
* Fix missing file when installing from source via PyPI_. `Issue 181`_
* Fix ``lint`` to verify ``data`` inputs specify a ``format`` attribute.
  8117e03_
* Docstring fix thanks to `@peterjc`_. fe7ad46_


---------------------
0.9.0 (2015-05-03)
---------------------

* Add new logo to the README thanks to `@petrkadlec`_ from `puradesign.cz
  <http://puradesign.cz/en>`__ and `@carlfeberhard`_ from the Galaxy Project.
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
  repository resolution - bug reported by `@peterjc`_. `Issue 158`_
* Improved exception handling for tool shed operations including new 
  ``--fail_fast`` command-line option. * `Issue 114`_, `Pull Request 173`_
* Implement more validation when using the ``shed_init`` command. 1cd0e2d_
* Add ``-r/--recursive`` option to ``shed_download`` and ``shed_diff`` 
  commands and allow these commands to work with ``.shed.yml`` files defining
  multipe repositories. 40a1f57_
* Add ``--port`` option to the ``serve`` and ``tool_factory`` commands.
  15804be_
* Fix problem introduced with `setup.py`_ during the 0.9.0 development cycle
  - thanks to `@peterjc`_. `Pull Request 171`_
* Fix clone bug introduced during 0.9.0 development cycle - thanks to
  `@bgruening`_. `Pull Request 175`_

---------------------
0.8.4 (2015-04-30)
---------------------

* Fix for Travis CI testing picking up invalid tests (reported by `@takadonet`_). `Issue 161`_
* Fix tar ordering for consistency (always sort by name) - thanks to `@peterjc`_.  `Pull Request 164`_, `Issue 159`_
* Fix exception handling related to tool shed operations - thanks to `@peterjc`_. `Pull Request 155`_, b86fe1f_

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

* Fix bug with ``config_init`` command thanks to `@bgruening`_. `Pull Request 151`_
* Fix unnessecary ``lint`` warning about ``parallelism`` tag reported by
  `@peterjc`_. 9bf1eab_

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
  maintained by `@JeanFred`_. This requires the lxml_ library to be
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
* Fix typo in ``test`` help thanks to first time contributor `@pvanheus`_.
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
* Update travis_test to install stable planemo from PyPI_. 39fedd2_
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

* Fix `setup.py`_ for installing non-Python data from PyPI_ (required newer
  for ``tool_factory`` command and reStructuredText linting). Thanks to
  Damion Dooley for the bug report. `Issue 83`_

---------------------
0.4.1 (2015-02-16)
---------------------

* Fix README.rst so it renders properly on PyPI_.

---------------------
0.4.0 (2015-02-16)
---------------------

* Implement ``tool_init`` command for bootstrapping creation of new
  tools (with `tutorial <http://planemo.readthedocs.org/en/latest/writing.html>`_.) 78f8274_
* Implement ``normalize`` command for reorganizing tool XML and macro
  debugging. e8c1d45_
* Implement ``tool_factory`` command to spin up Galaxy pre-configured the
  `Tool Factory`_. 9e746b4_
* Added basic linting of ``command`` blocks. b8d90ab_
* Improved linting of ``help`` blocks, including verifying valid
  `reStructuredText`. 411a8da_
* Fix bug related to ``serve`` command not killing Galaxy properly when complete. 53a6766_
* Have ``serve`` command display tools at the top level instead of in shallow sections. badc25f_
* Add additional dependencies to ``setup.py`` more functionality works out
  of the box. 85b9614_
* Fix terrible error message related to bioblend_ being unavailable.
  `Issue 70`_
* Various smaller documentation and project structure improvements.

---------------------
0.3.1 (2015-02-15)
---------------------

* Fixes to get PyPI_ workflow working properly.

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
* Various documentation fixes (thanks to Martin Čech). 36f7cb1_, b9232e5_
* Various smaller fixes for Docker support, tool linting, and documentation.

---------------------
0.0.1 (2014-10-04)
---------------------

* Initial work on the project - commands for testing, linting, serving Galaxy
  tools - and more experimental features involving Docker and Homebrew. 7d07782_

.. github_links
.. _Pull Request 684: https://github.com/galaxyproject/planemo/pull/684
.. _Pull Request 687: https://github.com/galaxyproject/planemo/pull/687
.. _Pull Request 688: https://github.com/galaxyproject/planemo/pull/688
.. _Pull Request 689: https://github.com/galaxyproject/planemo/pull/689
.. _Pull Request 690: https://github.com/galaxyproject/planemo/pull/690
.. _Issue 494: https://github.com/galaxyproject/planemo/issues/494
.. _Pull Request 665: https://github.com/galaxyproject/planemo/pull/665
.. _Pull Request 675: https://github.com/galaxyproject/planemo/pull/675
.. _Pull Request 680: https://github.com/galaxyproject/planemo/pull/680
.. _Pull Request 682: https://github.com/galaxyproject/planemo/pull/682
.. _Pull Request 683: https://github.com/galaxyproject/planemo/pull/683
.. _Pull Request 662: https://github.com/galaxyproject/planemo/pull/662
.. _Pull Request 665: https://github.com/galaxyproject/planemo/pull/665
.. _Pull Request 664: https://github.com/galaxyproject/planemo/pull/664
.. _84c4a73: https://github.com/galaxyproject/planemo/commit/84c4a73
.. _32f41c9: https://github.com/galaxyproject/planemo/commit/32f41c9
.. _3129216: https://github.com/galaxyproject/planemo/commit/3129216
.. _b4ae44d: https://github.com/galaxyproject/planemo/commit/b4ae44d
.. _3c95b7b: https://github.com/galaxyproject/planemo/commit/3c95b7b
.. _06bcf19: https://github.com/galaxyproject/planemo/commit/06bcf19
.. _525de8f: https://github.com/galaxyproject/planemo/commit/525de8f
.. _9867e56: https://github.com/galaxyproject/planemo/commit/9867e56
.. _ce0dc4e: https://github.com/galaxyproject/planemo/commit/ce0dc4e
.. _4c0f100: https://github.com/galaxyproject/planemo/commit/4c0f100
.. _04238d3: https://github.com/galaxyproject/planemo/commit/04238d3
.. _ced5ce2: https://github.com/galaxyproject/planemo/commit/ced5ce2
.. _9ab4a0d: https://github.com/galaxyproject/planemo/commit/9ab4a0d
.. _Pull Request 640: https://github.com/galaxyproject/planemo/pull/640
.. _0a1abfe: https://github.com/galaxyproject/planemo/commit/0a1abfe
.. _Pull Request 649: https://github.com/galaxyproject/planemo/pull/649
.. _Issue 620: https://github.com/galaxyproject/planemo/issues/620
.. _Pull Request 643: https://github.com/galaxyproject/planemo/pull/643
.. _Pull Request 642: https://github.com/galaxyproject/planemo/pull/642
.. _Pull Request 641: https://github.com/galaxyproject/planemo/pull/641
.. _Pull Request 639: https://github.com/galaxyproject/planemo/pull/639
.. _Pull Request 635: https://github.com/galaxyproject/planemo/pull/635
.. _Issue 633: https://github.com/galaxyproject/planemo/issues/633
.. _Issue 260: https://github.com/galaxyproject/planemo/issues/260
.. _Pull Request 638: https://github.com/galaxyproject/planemo/pull/638
.. _6638caa: https://github.com/galaxyproject/planemo/commit/6638caa
.. _8faf661: https://github.com/galaxyproject/planemo/commit/8faf661
.. _e343b67: https://github.com/galaxyproject/planemo/commit/e343b67
.. _Pull Request 634: https://github.com/galaxyproject/planemo/pull/634
.. _84ebc1f: https://github.com/galaxyproject/planemo/commit/84ebc1f
.. _03c9658: https://github.com/galaxyproject/planemo/commit/03c9658
.. _08c067c: https://github.com/galaxyproject/planemo/commit/08c067c
.. _fca4183: https://github.com/galaxyproject/planemo/commit/fca4183
.. _Issue 564: https://github.com/galaxyproject/planemo/issues/564
.. _Pull Request 628: https://github.com/galaxyproject/planemo/pull/628
.. _Issue 515: https://github.com/galaxyproject/planemo/issues/515
.. _Pull Request 629: https://github.com/galaxyproject/planemo/pull/629
.. _Pull Request 614: https://github.com/galaxyproject/planemo/pull/614
.. _32acd68: https://github.com/galaxyproject/planemo/commit/32acd68
.. _4595953: https://github.com/galaxyproject/planemo/commit/4595953
.. _Pull Request 612: https://github.com/galaxyproject/planemo/pull/612
.. _Issue 388: https://github.com/galaxyproject/planemo/issues/388
.. _Pull Request 610: https://github.com/galaxyproject/planemo/pull/610
.. _Pull Request 609: https://github.com/galaxyproject/planemo/pull/609
.. _Pull Request 608: https://github.com/galaxyproject/planemo/pull/608
.. _Pull Request 605: https://github.com/galaxyproject/planemo/pull/605
.. _Pull Request 606: https://github.com/galaxyproject/planemo/pull/606
.. _Pull Request 602: https://github.com/galaxyproject/planemo/pull/602
.. _Pull Request 570: https://github.com/galaxyproject/planemo/pull/570
.. _9228416: https://github.com/galaxyproject/planemo/commit/9228416
.. _50d3c4a: https://github.com/galaxyproject/planemo/commit/50d3c4a
.. _Issue 578: https://github.com/galaxyproject/planemo/issues/578
.. _Pull Request 591: https://github.com/galaxyproject/planemo/pull/591
.. _Pull Request 590: https://github.com/galaxyproject/planemo/pull/590
.. _f0da66f: https://github.com/galaxyproject/planemo/commit/f0da66f
.. _19b2ee9: https://github.com/galaxyproject/planemo/commit/19b2ee9
.. _9da8387: https://github.com/galaxyproject/planemo/commit/9da8387
.. _08cef54: https://github.com/galaxyproject/planemo/commit/08cef54
.. _2e4e5fc: https://github.com/galaxyproject/planemo/commit/2e4e5fc
.. _2e4e5fc: https://github.com/galaxyproject/planemo/commit/2e4e5fc
.. _Issue 573: https://github.com/galaxyproject/planemo/issues/573
.. _Pull Request 579: https://github.com/galaxyproject/planemo/pull/579
.. _ccdd2d5: https://github.com/galaxyproject/planemo/commit/ccdd2d5
.. _e925ba1: https://github.com/galaxyproject/planemo/commit/e925ba1
.. _ea5324f: https://github.com/galaxyproject/planemo/commit/ea5324f
.. _ca88b0c: https://github.com/galaxyproject/planemo/commit/ca88b0c
.. _b6d8294: https://github.com/galaxyproject/planemo/commit/b6d8294
.. _6a6f164: https://github.com/galaxyproject/planemo/commit/6a6f164
.. _d6da3a8: https://github.com/galaxyproject/planemo/commit/d6da3a8
.. _1ef05d2: https://github.com/galaxyproject/planemo/commit/1ef05d2
.. _7cca2e4: https://github.com/galaxyproject/planemo/commit/7cca2e4
.. _34538de: https://github.com/galaxyproject/planemo/commit/34538de
.. _6f91719: https://github.com/galaxyproject/planemo/commit/6f91719
.. _Pull Request 566: https://github.com/galaxyproject/planemo/pull/566
.. _Pull Request 559: https://github.com/galaxyproject/planemo/pull/559
.. _Pull Request 561: https://github.com/galaxyproject/planemo/pull/561
.. _Pull Request 565: https://github.com/galaxyproject/planemo/pull/565
.. _Pull Request 555: https://github.com/galaxyproject/planemo/pull/555
.. _a8e797b: https://github.com/galaxyproject/planemo/commit/a8e797b
.. _6c03de8: https://github.com/galaxyproject/planemo/commit/6c03de8
.. _ef4b9f4: https://github.com/galaxyproject/planemo/commit/ef4b9f4
.. _f7b6c7e: https://github.com/galaxyproject/planemo/commit/f7b6c7e
.. _07d94bd: https://github.com/galaxyproject/planemo/commit/07d94bd
.. _ca19910: https://github.com/galaxyproject/planemo/commit/ca19910
.. _24008ab: https://github.com/galaxyproject/planemo/commit/24008ab
.. _ce44e87: https://github.com/galaxyproject/planemo/commit/ce44e87
.. _Pull Request 550: https://github.com/galaxyproject/planemo/pull/550
.. _c2204b3: https://github.com/galaxyproject/planemo/commit/c2204b3
.. _c262b6d: https://github.com/galaxyproject/planemo/commit/c262b6d
.. _Pull Request 539: https://github.com/galaxyproject/planemo/pull/539
.. _Pull Request 541: https://github.com/galaxyproject/planemo/pull/541
.. _Pull Request 540: https://github.com/galaxyproject/planemo/pull/540
.. _Pull Request 545: https://github.com/galaxyproject/planemo/pull/545
.. _Pull Request 546: https://github.com/galaxyproject/planemo/pull/546
.. _3ceaa40: https://github.com/galaxyproject/planemo/commit/3ceaa40
.. _057f4f0: https://github.com/galaxyproject/planemo/commit/057f4f0
.. _9fdf490: https://github.com/galaxyproject/planemo/commit/9fdf490
.. _8c088c6: https://github.com/galaxyproject/planemo/commit/8c088c6
.. _e2398fb: https://github.com/galaxyproject/planemo/commit/e2398fb
.. _Pull Request 526: https://github.com/galaxyproject/planemo/pull/526
.. _Pull Request 528: https://github.com/galaxyproject/planemo/pull/528
.. _Pull Request 525: https://github.com/galaxyproject/planemo/pull/525
.. _a811e65: https://github.com/galaxyproject/planemo/commit/a811e65
.. _Pull Request 521: https://github.com/galaxyproject/planemo/pull/521
.. _a066afb: https://github.com/galaxyproject/planemo/commit/a066afb
.. _Pull Request 512: https://github.com/galaxyproject/planemo/pull/512
.. _08bb354: https://github.com/galaxyproject/planemo/commit/08bb354
.. _Pull Request 513: https://github.com/galaxyproject/planemo/pull/513
.. _Pull Request 510: https://github.com/galaxyproject/planemo/pull/510
.. _e890ab5: https://github.com/galaxyproject/planemo/commit/e890ab5
.. _Pull Request 507: https://github.com/galaxyproject/planemo/pull/507
.. _d53bcd6: https://github.com/galaxyproject/planemo/commit/d53bcd6
.. _725b232: https://github.com/galaxyproject/planemo/commit/725b232
.. _Pull Request 498: https://github.com/galaxyproject/planemo/pull/498
.. _01584c5: https://github.com/galaxyproject/planemo/commit/01584c5
.. _01f2af9: https://github.com/galaxyproject/planemo/commit/01f2af9
.. _0298510: https://github.com/galaxyproject/planemo/commit/0298510
.. _02a08a0: https://github.com/galaxyproject/planemo/commit/02a08a0
.. _05cc9f4: https://github.com/galaxyproject/planemo/commit/05cc9f485ee87bc344e3f43bb1cfd025a16a6247
.. _069e7ba: https://github.com/galaxyproject/planemo/commit/069e7ba
.. _08de8de: https://github.com/galaxyproject/planemo/commit/08de8de
.. _0bd4ff0: https://github.com/galaxyproject/planemo/commit/0bd4ff0
.. _0e4f70a: https://github.com/galaxyproject/planemo/commit/0e4f70a
.. _0f8cb10: https://github.com/galaxyproject/planemo/commit/0f8cb10
.. _13a5ae7: https://github.com/galaxyproject/planemo/commit/13a5ae7
.. _15804be: https://github.com/galaxyproject/planemo/commit/15804be
.. _15d33c7: https://github.com/galaxyproject/planemo/commit/15d33c7
.. _182fe57: https://github.com/galaxyproject/planemo/commit/182fe57
.. _1927168: https://github.com/galaxyproject/planemo/commit/1927168
.. _1982076: https://github.com/galaxyproject/planemo/commit/1982076
.. _19900a6: https://github.com/galaxyproject/planemo/commit/19900a6
.. _1a157d4: https://github.com/galaxyproject/planemo/commit/1a157d4
.. _1a85493: https://github.com/galaxyproject/planemo/commit/1a85493
.. _1c6cfbb: https://github.com/galaxyproject/planemo/commit/1c6cfbb
.. _1c7ee5b: https://github.com/galaxyproject/planemo/commit/1c7ee5b
.. _1cd0e2d: https://github.com/galaxyproject/planemo/commit/1cd0e2d
.. _1e3668a: https://github.com/galaxyproject/planemo/commit/1e3668a
.. _2052db0: https://github.com/galaxyproject/planemo/commit/2052db0
.. _20a8680: https://github.com/galaxyproject/planemo/commit/20a86807cb7ea87db2dbc0197ae08a40df3ab2bc
.. _21bb463: https://github.com/galaxyproject/planemo/commit/21bb463ad6c321bcb669603049a5e89a69766ad9
.. _25ef0d5: https://github.com/galaxyproject/planemo/commit/25ef0d5
.. _26e378e: https://github.com/galaxyproject/planemo/commit/26e378e
.. _26e3cdb: https://github.com/galaxyproject/planemo/commit/26e3cdb
.. _2a7c792: https://github.com/galaxyproject/planemo/commit/2a7c792
.. _2ae2b49: https://github.com/galaxyproject/planemo/commit/2ae2b49
.. _2e41e0a: https://github.com/galaxyproject/planemo/commit/2e41e0a
.. _2f66fc3: https://github.com/galaxyproject/planemo/commit/2f66fc3
.. _30a9c3f: https://github.com/galaxyproject/planemo/commit/30a9c3f
.. _32c6e7f: https://github.com/galaxyproject/planemo/commit/32c6e7f78bb8f04d27615cfd8948b0b89f27b4e6
.. _33294d2: https://github.com/galaxyproject/planemo/commit/33294d2
.. _334f2d4: https://github.com/galaxyproject/planemo/commit/334f2d4
.. _343902d: https://github.com/galaxyproject/planemo/commit/343902d
.. _3499ca0: https://github.com/galaxyproject/planemo/commit/3499ca0a15affcaf8ac9efc55880da40b0626679
.. _358a42c: https://github.com/galaxyproject/planemo/commit/358a42c
.. _36ac6d8: https://github.com/galaxyproject/planemo/commit/36ac6d8
.. _36f7cb1: https://github.com/galaxyproject/planemo/commit/36f7cb114f77731f90860d513a930e10ce5c1ba5
.. _37dcc07: https://github.com/galaxyproject/planemo/commit/37dcc07
.. _39fedd2: https://github.com/galaxyproject/planemo/commit/39fedd2
.. _3f4ab44: https://github.com/galaxyproject/planemo/commit/3f4ab44
.. _40a1f57: https://github.com/galaxyproject/planemo/commit/40a1f57
.. _411a8da: https://github.com/galaxyproject/planemo/commit/411a8da21c92ba37c7ad95bfce9928d9b8fd998e
.. _44de95c: https://github.com/galaxyproject/planemo/commit/44de95c0d7087a5822941959f9a062f6382e329b
.. _45135ff: https://github.com/galaxyproject/planemo/commit/45135ff
.. _4823c5e: https://github.com/galaxyproject/planemo/commit/4823c5e
.. _49c5c1e: https://github.com/galaxyproject/planemo/commit/49c5c1e
.. _4c71299: https://github.com/galaxyproject/planemo/commit/4c71299
.. _4cd571c: https://github.com/galaxyproject/planemo/commit/4cd571c
.. _4d29bf1: https://github.com/galaxyproject/planemo/commit/4d29bf1
.. _4d6f7d9: https://github.com/galaxyproject/planemo/commit/4d6f7d9
.. _4e1377c: https://github.com/galaxyproject/planemo/commit/4e1377c
.. _4f61025: https://github.com/galaxyproject/planemo/commit/4f61025
.. _508dce7: https://github.com/galaxyproject/planemo/commit/508dce7
.. _53a6766: https://github.com/galaxyproject/planemo/commit/53a6766cdebdddc976189f6dc6a264bb4105c4bf
.. _53edd99: https://github.com/galaxyproject/planemo/commit/53edd99
.. _552059f: https://github.com/galaxyproject/planemo/commit/552059f
.. _572e754: https://github.com/galaxyproject/planemo/commit/572e754
.. _5b97f2e: https://github.com/galaxyproject/planemo/commit/5b97f2e
.. _5d08b67: https://github.com/galaxyproject/planemo/commit/5d08b67
.. _5d7db92: https://github.com/galaxyproject/planemo/commit/5d7db92
.. _5e0b6d1: https://github.com/galaxyproject/planemo/commit/5e0b6d1
.. _63cd431: https://github.com/galaxyproject/planemo/commit/63cd431
.. _63e456c: https://github.com/galaxyproject/planemo/commit/63e456c
.. _6409449: https://github.com/galaxyproject/planemo/commit/6409449
.. _6514ff5: https://github.com/galaxyproject/planemo/commit/6514ff5
.. _65b999d: https://github.com/galaxyproject/planemo/commit/65b999d21bacc133a80ecf5f61e0728afec66ccc
.. _6b4e7a6: https://github.com/galaxyproject/planemo/commit/6b4e7a6
.. _6bcf699: https://github.com/galaxyproject/planemo/commit/6bcf699
.. _6d0f502: https://github.com/galaxyproject/planemo/commit/6d0f502
.. _6d81a94: https://github.com/galaxyproject/planemo/commit/6d81a94
.. _6e1e726: https://github.com/galaxyproject/planemo/commit/6e1e726
.. _7572e99: https://github.com/galaxyproject/planemo/commit/7572e99
.. _764ce01: https://github.com/galaxyproject/planemo/commit/764ce01
.. _775bf49: https://github.com/galaxyproject/planemo/commit/775bf49
.. _776773c: https://github.com/galaxyproject/planemo/commit/776773c
.. _78f8274: https://github.com/galaxyproject/planemo/commit/78f82747996e4a28f96c85ad72efe5e54c8c74bd
.. _7be1bf5: https://github.com/galaxyproject/planemo/commit/7be1bf5
.. _7c69bf6: https://github.com/galaxyproject/planemo/commit/7c69bf6
.. _7d07782: https://github.com/galaxyproject/planemo/commit/7d077828559c9c9c352ac814f9e3b86b1b3a2a9f
.. _8117e03: https://github.com/galaxyproject/planemo/commit/8117e03
.. _8207026: https://github.com/galaxyproject/planemo/commit/8207026
.. _826d371: https://github.com/galaxyproject/planemo/commit/826d371
.. _82e8b1f: https://github.com/galaxyproject/planemo/commit/82e8b1f17eae526aeb341cb4fffb8d09d73bb419
.. _834bfb2: https://github.com/galaxyproject/planemo/commit/834bfb2929d367892a3abe9c0b88d5a0277d7905
.. _83e227a: https://github.com/galaxyproject/planemo/commit/83e227a
.. _85b9614: https://github.com/galaxyproject/planemo/commit/85b961465f46351507f80ddc3758349535060502
.. _89674cb: https://github.com/galaxyproject/planemo/commit/89674cb
.. _8b31850: https://github.com/galaxyproject/planemo/commit/8b31850
.. _8e96864: https://github.com/galaxyproject/planemo/commit/8e96864
.. _8eda729: https://github.com/galaxyproject/planemo/commit/8eda729
.. _912df02: https://github.com/galaxyproject/planemo/commit/912df02
.. _916f610: https://github.com/galaxyproject/planemo/commit/916f610
.. _93b7bda: https://github.com/galaxyproject/planemo/commit/93b7bda
.. _94097c7: https://github.com/galaxyproject/planemo/commit/94097c7
.. _9427b47: https://github.com/galaxyproject/planemo/commit/9427b47
.. _94645ed: https://github.com/galaxyproject/planemo/commit/94645ed
.. _949a36d: https://github.com/galaxyproject/planemo/commit/949a36d
.. _95d5cba: https://github.com/galaxyproject/planemo/commit/95d5cba
.. _965511d: https://github.com/galaxyproject/planemo/commit/965511d
.. _988de1d: https://github.com/galaxyproject/planemo/commit/988de1d
.. _98fdc8c: https://github.com/galaxyproject/planemo/commit/98fdc8c
.. _99ee51a: https://github.com/galaxyproject/planemo/commit/99ee51a
.. _9bf1eab: https://github.com/galaxyproject/planemo/commit/9bf1eab
.. _9e746b4: https://github.com/galaxyproject/planemo/commit/9e746b455e3b15219878cddcdeda722979639401
.. _9f3957d: https://github.com/galaxyproject/planemo/commit/9f3957d
.. _9ff0d2d: https://github.com/galaxyproject/planemo/commit/9ff0d2d
.. _CWL: http://www.commonwl.org/
.. _Issue 108: https://github.com/galaxyproject/planemo/issues/108
.. _Issue 114: https://github.com/galaxyproject/planemo/issues/114
.. _Issue 118: https://github.com/galaxyproject/planemo/issues/118
.. _Issue 124: https://github.com/galaxyproject/planemo/issues/#124
.. _Issue 138: https://github.com/galaxyproject/planemo/issues/#138
.. _Issue 139: https://github.com/galaxyproject/planemo/issues/139
.. _Issue 150: https://github.com/galaxyproject/planemo/issues/150
.. _Issue 156: https://github.com/galaxyproject/planemo/issues/156
.. _Issue 158: https://github.com/galaxyproject/planemo/issues/158
.. _Issue 159: https://github.com/galaxyproject/planemo/issues/159
.. _Issue 15: https://github.com/galaxyproject/planemo/issues/15
.. _Issue 161: https://github.com/galaxyproject/planemo/issues/161
.. _Issue 167: https://github.com/galaxyproject/planemo/issues/167
.. _Issue 168: https://github.com/galaxyproject/planemo/issues/168
.. _Issue 169: https://github.com/galaxyproject/planemo/issues/169
.. _Issue 170: https://github.com/galaxyproject/planemo/issues/170
.. _Issue 176: https://github.com/galaxyproject/planemo/issues/176
.. _Issue 179: https://github.com/galaxyproject/planemo/issues/179
.. _Issue 180: https://github.com/galaxyproject/planemo/issues/180
.. _Issue 181: https://github.com/galaxyproject/planemo/issues/181
.. _Issue 184: https://github.com/galaxyproject/planemo/issues/184
.. _Issue 186: https://github.com/galaxyproject/planemo/issues/186
.. _Issue 188: https://github.com/galaxyproject/planemo/issues/188
.. _Issue 189: https://github.com/galaxyproject/planemo/issues/189
.. _Issue 192: https://github.com/galaxyproject/planemo/issues/192
.. _Issue 205: https://github.com/galaxyproject/planemo/issues/205
.. _Issue 206: https://github.com/galaxyproject/planemo/issues/206
.. _Issue 209: https://github.com/galaxyproject/planemo/issues/209
.. _Issue 211: https://github.com/galaxyproject/planemo/issues/211
.. _Issue 217: https://github.com/galaxyproject/planemo/issues/217
.. _Issue 223: https://github.com/galaxyproject/planemo/issues/223
.. _Issue 231: https://github.com/galaxyproject/planemo/issues/231
.. _Issue 233: https://github.com/galaxyproject/planemo/issues/233
.. _Issue 240: https://github.com/galaxyproject/planemo/issues/240
.. _Issue 243: https://github.com/galaxyproject/planemo/issues/243
.. _Issue 245: https://github.com/galaxyproject/planemo/issues/245
.. _Issue 246: https://github.com/galaxyproject/planemo/issues/246
.. _Issue 272: https://github.com/galaxyproject/planemo/issues/272
.. _Issue 303: https://github.com/galaxyproject/planemo/issues/303
.. _Issue 313: https://github.com/galaxyproject/planemo/issues/313
.. _Issue 329: https://github.com/galaxyproject/planemo/issues/329
.. _Issue 333: https://github.com/galaxyproject/planemo/issues/333
.. _Issue 361: https://github.com/galaxyproject/planemo/issues/361
.. _Issue 362: https://github.com/galaxyproject/planemo/issues/362
.. _Issue 373: https://github.com/galaxyproject/planemo/issues/373
.. _Issue 391: https://github.com/galaxyproject/planemo/issues/391
.. _Issue 416: https://github.com/galaxyproject/planemo/issues/416
.. _Issue 420: https://github.com/galaxyproject/planemo/issues/420
.. _Issue 475: https://github.com/galaxyproject/planemo/issues/475
.. _Issue 61: https://github.com/galaxyproject/planemo/issues/61
.. _Issue 70: https://github.com/galaxyproject/planemo/issues/70
.. _Issue 78: https://github.com/galaxyproject/planemo/issues/78
.. _Issue 80: https://github.com/galaxyproject/planemo/issues/80
.. _Issue 83: https://github.com/galaxyproject/planemo/issues/83
.. _Issue 88: https://github.com/galaxyproject/planemo/issues/88
.. _Issue 89: https://github.com/galaxyproject/planemo/issues/#89
.. _Issue 91: https://github.com/galaxyproject/planemo/issues/#91
.. _Pull Request 101: https://github.com/galaxyproject/planemo/pull/101
.. _Pull Request 102: https://github.com/galaxyproject/planemo/pull/102
.. _Pull Request 104: https://github.com/galaxyproject/planemo/pull/104
.. _Pull Request 111: https://github.com/galaxyproject/planemo/pull/111
.. _Pull Request 113: https://github.com/galaxyproject/planemo/pull/113
.. _Pull Request 129: https://github.com/galaxyproject/planemo/pull/129
.. _Pull Request 130: https://github.com/galaxyproject/planemo/pull/130
.. _Pull Request 143: https://github.com/galaxyproject/planemo/pull/143
.. _Pull Request 151: https://github.com/galaxyproject/planemo/pull/151
.. _Pull Request 155: https://github.com/galaxyproject/planemo/pull/155
.. _Pull Request 164: https://github.com/galaxyproject/planemo/pull/164
.. _Pull Request 171: https://github.com/galaxyproject/planemo/pull/171
.. _Pull Request 173: https://github.com/galaxyproject/planemo/pull/173
.. _Pull Request 175: https://github.com/galaxyproject/planemo/pull/175
.. _Pull Request 185: https://github.com/galaxyproject/planemo/pull/185
.. _Pull Request 186: https://github.com/galaxyproject/planemo/pull/186
.. _Pull Request 187: https://github.com/galaxyproject/planemo/pull/187
.. _Pull Request 190: https://github.com/galaxyproject/planemo/pull/190
.. _Pull Request 196: https://github.com/galaxyproject/planemo/pull/196
.. _Pull Request 1: https://github.com/galaxyproject/planemo/pull/1
.. _Pull Request 200: https://github.com/galaxyproject/planemo/pull/200
.. _Pull Request 203: https://github.com/galaxyproject/planemo/pull/203
.. _Pull Request 206: https://github.com/galaxyproject/planemo/pull/206
.. _Pull Request 207: https://github.com/galaxyproject/planemo/pull/207
.. _Pull Request 208: https://github.com/galaxyproject/planemo/pull/208
.. _Pull Request 210: https://github.com/galaxyproject/planemo/pull/210
.. _Pull Request 213: https://github.com/galaxyproject/planemo/pull/213
.. _Pull Request 215: https://github.com/galaxyproject/planemo/pull/215
.. _Pull Request 216: https://github.com/galaxyproject/planemo/pull/216
.. _Pull Request 22: https://github.com/galaxyproject/planemo/pull/22
.. _Pull Request 230: https://github.com/galaxyproject/planemo/pull/230
.. _Pull Request 235: https://github.com/galaxyproject/planemo/pull/235
.. _Pull Request 23: https://github.com/galaxyproject/planemo/pull/23
.. _Pull Request 251: https://github.com/galaxyproject/planemo/pull/251
.. _Pull Request 253: https://github.com/galaxyproject/planemo/pull/253
.. _Pull Request 254: https://github.com/galaxyproject/planemo/pull/254
.. _Pull Request 255: https://github.com/galaxyproject/planemo/pull/255
.. _Pull Request 256: https://github.com/galaxyproject/planemo/pull/256
.. _Pull Request 277: https://github.com/galaxyproject/planemo/pull/277
.. _Pull Request 278: https://github.com/galaxyproject/planemo/pull/278
.. _Pull Request 284: https://github.com/galaxyproject/planemo/pull/284
.. _Pull Request 285: https://github.com/galaxyproject/planemo/pull/285
.. _Pull Request 287: https://github.com/galaxyproject/planemo/pull/287
.. _Pull Request 292: https://github.com/galaxyproject/planemo/pull/292
.. _Pull Request 297: https://github.com/galaxyproject/planemo/pull/297
.. _Pull Request 29: https://github.com/galaxyproject/planemo/pull/29
.. _Pull Request 301: https://github.com/galaxyproject/planemo/pull/301
.. _Pull Request 304: https://github.com/galaxyproject/planemo/pull/304
.. _Pull Request 305: https://github.com/galaxyproject/planemo/pull/305
.. _Pull Request 307: https://github.com/galaxyproject/planemo/pull/307
.. _Pull Request 309: https://github.com/galaxyproject/planemo/pull/309
.. _Pull Request 310: https://github.com/galaxyproject/planemo/pull/310
.. _Pull Request 311: https://github.com/galaxyproject/planemo/pull/311
.. _Pull Request 312: https://github.com/galaxyproject/planemo/pull/312
.. _Pull Request 314: https://github.com/galaxyproject/planemo/pull/314
.. _Pull Request 316: https://github.com/galaxyproject/planemo/pull/316
.. _Pull Request 322: https://github.com/galaxyproject/planemo/pull/322
.. _Pull Request 327: https://github.com/galaxyproject/planemo/pull/327
.. _Pull Request 330: https://github.com/galaxyproject/planemo/pull/330
.. _Pull Request 333: https://github.com/galaxyproject/planemo/pull/333
.. _Pull Request 334: https://github.com/galaxyproject/planemo/pull/334
.. _Pull Request 335: https://github.com/galaxyproject/planemo/pull/335
.. _Pull Request 339: https://github.com/galaxyproject/planemo/pull/339
.. _Pull Request 33: https://github.com/galaxyproject/planemo/pull/33
.. _Pull Request 343: https://github.com/galaxyproject/planemo/pull/343
.. _Pull Request 344: https://github.com/galaxyproject/planemo/pull/344
.. _Pull Request 345: https://github.com/galaxyproject/planemo/pull/345
.. _Pull Request 350: https://github.com/galaxyproject/planemo/pull/350
.. _Pull Request 351: https://github.com/galaxyproject/planemo/pull/351
.. _Pull Request 356: https://github.com/galaxyproject/planemo/pull/356
.. _Pull Request 375: https://github.com/galaxyproject/planemo/pull/375
.. _Pull Request 390: https://github.com/galaxyproject/planemo/pull/390
.. _Pull Request 394: https://github.com/galaxyproject/planemo/pull/394
.. _Pull Request 398: https://github.com/galaxyproject/planemo/pull/398
.. _Pull Request 403: https://github.com/galaxyproject/planemo/pull/403
.. _Pull Request 40: https://github.com/galaxyproject/planemo/pull/40
.. _Pull Request 425: https://github.com/galaxyproject/planemo/pull/425
.. _Pull Request 426: https://github.com/galaxyproject/planemo/pull/426
.. _Pull Request 428: https://github.com/galaxyproject/planemo/pull/428
.. _Pull Request 429: https://github.com/galaxyproject/planemo/pull/429
.. _Pull Request 454: https://github.com/galaxyproject/planemo/pull/454
.. _Pull Request 472: https://github.com/galaxyproject/planemo/pull/472
.. _Pull Request 479: https://github.com/galaxyproject/planemo/pull/479
.. _Pull Request 480: https://github.com/galaxyproject/planemo/pull/480
.. _Pull Request 484: https://github.com/galaxyproject/planemo/pull/484
.. _Pull Request 488: https://github.com/galaxyproject/planemo/pull/488
.. _Pull Request 491: https://github.com/galaxyproject/planemo/pull/491
.. _Pull Request 51: https://github.com/galaxyproject/planemo/pull/51
.. _Pull Request 56: https://github.com/galaxyproject/planemo/pull/56
.. _Pull Request 63: https://github.com/galaxyproject/planemo/pull/63
.. _Pull Request 68: https://github.com/galaxyproject/planemo/pull/68
.. _Pull Request 6: https://github.com/galaxyproject/planemo/pull/6
.. _Pull Request 98: https://github.com/galaxyproject/planemo/pull/98
.. _Pull Request 99: https://github.com/galaxyproject/planemo/pull/99
.. _XSD: http://www.w3schools.com/schema/
.. _a13a120: https://github.com/galaxyproject/planemo/commit/a13a120
.. _a2c13e4: https://github.com/galaxyproject/planemo/commit/a2c13e46259e3be35de1ecaae858ba818bb94734
.. _a2ee135: https://github.com/galaxyproject/planemo/commit/a2ee135
.. _a4110a8: https://github.com/galaxyproject/planemo/commit/a4110a85a770988e5cd3c31ccc9475717897d59c
.. _a4e6958: https://github.com/galaxyproject/planemo/commit/a4e6958
.. _a58a3b8: https://github.com/galaxyproject/planemo/commit/a58a3b8
.. _a87899b: https://github.com/galaxyproject/planemo/commit/a87899b
.. _aad1eed: https://github.com/galaxyproject/planemo/commit/aad1eed
.. _ac4f828: https://github.com/galaxyproject/planemo/commit/ac4f82898f7006799142503a33c3978428660ce7
.. _ad3b2f0: https://github.com/galaxyproject/planemo/commit/ad3b2f0
.. _af39061: https://github.com/galaxyproject/planemo/commit/af390612004dab636d8696839bb723d39f97c85d
.. _af7448c: https://github.com/galaxyproject/planemo/commit/af7448c
.. _b0b867e: https://github.com/galaxyproject/planemo/commit/b0b867e
.. _b1c8f1d: https://github.com/galaxyproject/planemo/commit/b1c8f1d
.. _b53006d: https://github.com/galaxyproject/planemo/commit/b53006d
.. _b757791: https://github.com/galaxyproject/planemo/commit/b757791
.. _b7d9e96: https://github.com/galaxyproject/planemo/commit/b7d9e96
.. _b86fe1f: https://github.com/galaxyproject/planemo/commit/b86fe1f
.. _b8d90ab: https://github.com/galaxyproject/planemo/commit/b8d90abab8bf53ae2e7cca4317223c01af9ab68c
.. _b9232e5: https://github.com/galaxyproject/planemo/commit/b9232e55e713abbd1d9ce8b0b34cbec6c701dc17
.. _badc25f: https://github.com/galaxyproject/planemo/commit/badc25fca495b61457ffb2e027f3fe9cf17c798f
.. _bioblend: https://github.com/galaxyproject/bioblend/
.. _c1713d2: https://github.com/galaxyproject/planemo/commit/c1713d2
.. _c23569f: https://github.com/galaxyproject/planemo/commit/c23569f
.. _c444855: https://github.com/galaxyproject/planemo/commit/c444855
.. _c4dfd55: https://github.com/galaxyproject/planemo/commit/c4dfd55
.. _c8640b6: https://github.com/galaxyproject/planemo/commit/c8640b6
.. _cb5b906: https://github.com/galaxyproject/planemo/commit/cb5b906
.. _cc1a447: https://github.com/galaxyproject/planemo/commit/cc1a447
.. _cc8abb6: https://github.com/galaxyproject/planemo/commit/cc8abb6
.. _ce8e1be: https://github.com/galaxyproject/planemo/commit/ce8e1be
.. _cwltool: https://github.com/common-workflow-language/cwltool/.. _d26929e: https://github.com/galaxyproject/planemo/commit/d26929e
.. _d26929e: https://github.com/galaxyproject/planemo/commit/d26929e
.. _d755fe7: https://github.com/galaxyproject/planemo/commit/d755fe7
.. _d76b489: https://github.com/galaxyproject/planemo/commit/d76b489
.. _d8f2038: https://github.com/galaxyproject/planemo/commit/d8f2038
.. _dad2d9d: https://github.com/galaxyproject/planemo/commit/dad2d9d
.. _dd94ddc: https://github.com/galaxyproject/planemo/commit/dd94ddc
.. _dff4f33: https://github.com/galaxyproject/planemo/commit/dff4f33c750a8dbe651c38e149a26dd42e706a82
.. _e0577e7: https://github.com/galaxyproject/planemo/commit/e0577e7
.. _e0acf91: https://github.com/galaxyproject/planemo/commit/e0acf91
.. _e38c436: https://github.com/galaxyproject/planemo/commit/e38c436
.. _e769118: https://github.com/galaxyproject/planemo/commit/e769118
.. _e8c1d45: https://github.com/galaxyproject/planemo/commit/e8c1d45f0c9a11bcf69ec2967836c3b8f432dd97
.. _eb039c0: https://github.com/galaxyproject/planemo/commit/eb039c0
.. _ec6e30f: https://github.com/galaxyproject/planemo/commit/ec6e30f
.. _efc5f30: https://github.com/galaxyproject/planemo/commit/efc5f30
.. _f0610d7: https://github.com/galaxyproject/planemo/commit/f0610d7
.. _f3394e7: https://github.com/galaxyproject/planemo/commit/f3394e7
.. _f3c6917: https://github.com/galaxyproject/planemo/commit/f3c6917
.. _f3cafaa: https://github.com/galaxyproject/planemo/commit/f3cafaa
.. _f7554d1: https://github.com/galaxyproject/planemo/commit/f7554d1
.. _f798c7e: https://github.com/galaxyproject/planemo/commit/f798c7e
.. _f854138: https://github.com/galaxyproject/planemo/commit/f854138
.. _f99f6c1: https://github.com/galaxyproject/planemo/commit/f99f6c1
.. _fba3874: https://github.com/galaxyproject/planemo/commit/fba3874
.. _fdb1b51: https://github.com/galaxyproject/planemo/commit/fdb1b51
.. _fdce74c: https://github.com/galaxyproject/planemo/commit/fdce74c
.. _fe7ad46: https://github.com/galaxyproject/planemo/commit/fe7ad46
.. _fea51fc: https://github.com/galaxyproject/planemo/commit/fea51fc
.. _lxml: http://lxml.de/
.. _nose: https://nose.readthedocs.org/en/latest/
.. _xmllint: http://xmlsoft.org/xmllint.html
.. _Conda: http://conda.pydata.org/
.. _Tool Factory: http://bioinformatics.oxfordjournals.org/content/early/2012/09/27/bioinformatics.bts573.full.pdf
.. _git: https://git-scm.com/
.. _anaconda-verify: https://github.com/ContinuumIO/anaconda-verify
.. _galaxy.xsd: https://github.com/galaxyproject/planemo/blob/master/planemo/xml/xsd/tool/galaxy.xsd
.. _setup.py: https://github.com/galaxyproject/planemo/blob/master/setup.py
.. _Bioconductor: https://www.bioconductor.org/
.. _tools-iuc: https://github.com/galaxyproject/tools-iuc
.. _PyPI: https://pypi.python.org/pypi
.. _Involucro: https://github.com/involucro/involucro
.. _Bioconda: https://bioconda.github.io/
.. _pip: https://pip.pypa.io/en/stable/
.. _quay.io: https://quay.io/
.. _galaxy-lib: https://github.com/galaxyproject/galaxy-lib
.. _@abretaud: https://github.com/abretaud
.. _@erasche: https://github.com/erasche
.. _@peterjc: https://github.com/peterjc
.. _@mr-c: https://github.com/mr-c
.. _@martenson: https://github.com/martenson
.. _@nsoranzo: https://github.com/nsoranzo
.. _@nturaga: https://github.com/nturaga
.. _@bgruening: https://github.com/bgruening
.. _@carlfeberhard: https://github.com/carlfeberhard
.. _@lparsons: https://github.com/lparsons
.. _@kellrott: https://github.com/kellrott
.. _@mvdbeek: https://github.com/mvdbeek
.. _@natefoo: https://github.com/natefoo
.. _@pstew: https://github.com/pstew
.. _@ramezrawas: https://github.com/ramezrawas
.. _@chambm: https://github.com/chambm
.. _@takadonet: https://github.com/takadonet
.. _@petrkadlec: https://github.com/petrkadlec
.. _@pvanheus: https://github.com/pvanheus
.. _@einon: https://github.com/einon
.. _@blankenberg: https://github.com/blankenberg
.. _@JeanFred: https://github.com/JeanFred
.. _@gregvonkuster: https://github.com/gregvonkuster
.. _@remimarenco: https://github.com/remimarenco
.. _@thriqon: https://github.com/thriqon
.. _@RJMW: https://github.com/RJMW

