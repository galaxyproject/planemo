Best Practices for Maintaining Galaxy Workflows
===============================================

There are a number of things the user interface of Galaxy will allow that are not
considered best practices because they make the workflow harder to test, use within
subworkflows and invocation reports, and consume via the API. These best practices
can be checking with:

::

    $ planemo workflow_lint path/to/workflow.ga

Workflow Structure
------------------

Outputs
~~~~~~~

Workflows should define explicit, labelled outputs. Galaxy doesn't require you to
declare outputs explicitly or label them - but doing so provides a lot of advantages.
A workflow with declared, labelled outputs specifies and explicit interface that is
much easier to consume when building a report for the workflow, testing the workflow,
using the workflow via the API, and using the workflow as a subworkflow in Galaxy.

TODO: picture of this in the UI.

Inputs
~~~~~~

Similarly to outputs and for similar reasons, all inputs should be explicit (with
labelled input nodes) and tool steps should not have disconnected data inputs (even
though the GUI can handle this) or consume workflow parameters. Runtime parameters
should only be used for post job actions and newer type workflow parameter inputs
should be used to manipulate tool logic.

TODO: picture of new input parameters

In addition to making the interface easier to use in the context of subworkflows,
the API, testing, etc.. future enhancements to Galaxy will allow a much simpler
UI for workflows that only use explicit input parameters in this fashion.

https://github.com/galaxyproject/galaxy/pull/9151

Tools
~~~~~

The tools used within a workflow should be packaged with Galaxy by default or published
to the main Galaxy ToolShed. Using private tool sheds or the test tool shed limits the
ability of other Galaxy's to use the workflow.

Tests
-----

Writing workflow tests allows consumers of your workflow to know it works in their
Galaxy environment and can allow for richer continuous integration (CI). Check out
the `Planemo Test Format <http://planemo.readthedocs.io/en/latest/test_format.html>`__
documentation for more information on the format and how to test workflows with Planemo.

Planemo can help you stud out tests for a workflow developed within the UI quickly
with the ``workflow_test_init`` command.

::

    $ planemo workflow_test_init path/to/workflow.ga

Publishing 
----------

Unlike with Galaxy tools - the Galaxy team doesn't endorse a specific registry for
Galaxy workflows. But also unlike Galaxy tools, any user can just paste a URL for
a workflow right into the user interface so sharing a workflow can be as easy as
passing around a GitHub link.

Github
~~~~~~

Even if you're publishing your workflows to other registries or website, we always
recommend publishing workflows to Github (or a publicly available Gitlab server).

Dockstore
~~~~~~~~~

A repository containing Galaxy workflows and published to GitHub can be registered
with `Dockstore <https://dockstore.org/>`__. This allows others to search for the
workflow and access it using standard GA4GH APIs. In the future, deep bi-directional
integration between Galaxy and Dockstore will be available that will make these
workflows even more useful.

A ``.dockstore.yml`` file should be placed in the root of your workflow repository before
registering the repository with Dockstore. This will allow Dockstore to find your workflows
and their tests automatically.

Planemo can create this file for you by executing the ``dockstore_init`` command from
the root of your workflow repository

::

    $ planemo dockstore_init

Planemo's ``workflow_lint`` will check the contents of your ``.dockstore.yml`` file during
execution if this file is present.

WorkFlow Hub
~~~~~~~~~~~~

Coming soon.
