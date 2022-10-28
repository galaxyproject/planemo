
``workflow_job_init`` command
======================================

This section is auto-generated from the help text for the planemo command
``workflow_job_init``. This help message can be generated with ``planemo workflow_job_init
--help``.

**Usage**::

    planemo workflow_job_init [OPTIONS] WORKFLOW_PATH_OR_ID

**Help**

Initialize a Galaxy workflow job description for supplied workflow.

Be sure to your lint your workflow with ``workflow_lint`` before calling this
to ensure inputs and outputs comply with best practices that make workflow
testing easier.

Jobs can be run with the planemo run command (``planemo run workflow.ga job.yml``).
Planemo run works with Galaxy tools and CWL artifacts (both tools and workflows)
as well so this command may be renamed to to job_init at something along those
lines at some point.

**Options**::


      -f, --force                     Overwrite existing files if present.
      -o, --output FILE
      --galaxy_url TEXT               Remote Galaxy URL to use with external Galaxy
                                      engine.
      --galaxy_user_key TEXT          User key to use with external Galaxy engine.
      --from_invocation / --from_uri  Build a workflow test or job description from
                                      an invocation ID run on an external Galaxy.A
                                      Galaxy URL and API key must also be specified.
                                      This allows test data to be downloadedand
                                      inputs and parameters defined automatically.
                                      Alternatively, the default is to build
                                      thedescriptions from a provided workflow URI.
      --profile TEXT                  Name of profile (created with the
                                      profile_create command) to use with this
                                      command.
      --help                          Show this message and exit.
    
