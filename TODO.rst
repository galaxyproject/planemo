============
TODO
============

Critical Bugs
-------------

- Use repositories test-data instead of Galaxy's.
- planemo serve doesn't kill Galaxy server properly.

linting
-------

- Check test data files exist.
- Ensure supplied test parameters reference inputs.
- Compile RST in Help ensure compiles cleanly.

Galaxy Interactions
-------------------

- Handle tool data tables.
- Add --dependency_resolution.
  wih options as at least `brew`, `<path to resolution file>`
- serve should pre-create an admin user.
- serve should prepolate a history with all test data or at least
  an FTP directory.
- Option to run galaxy in Bjoern's Docker images instead of local.

Travis CI
---------

- Handle custom datatypes. 
- Handle tool data tables.

New Commands
-------------

- watch - retest and lint tools whenever they change.
- tool_init - Prompt for a bunch stuff and fill out form.


Other stuff
-----------

- project_init should also do a travis_init.
- travis_init should report more and have -f option.
