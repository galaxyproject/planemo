Workflow execution against an external Galaxy
===============================================

The first requirement for executing workflows on an external Galaxy server is
a user account for that server. If you don't already have one, `<https://usegalaxy.org>`_,
`<https://usegalaxy.eu>`_ and `<https://usegalaxy.org.au>`_ all provide free accounts which
can be used for this tutorial.

Assuming you have selected a server for this tutorial and have an account, you
need to retrieve the API key associated with that account. This can be found at
``{server_url}/user/api_key``, or by going to the 'User' dropdown menu, selecting
'Preferences' and then clicking on 'Manage API key'.

Now you can run the workflow with:

::

    $ planemo run tutorial.ga tutorial-job.yml --engine external_galaxy --galaxy_url SERVER_URL --galaxy_user_key YOUR_API_KEY


If you want to set the name of the history in which the workflow executes, add
``--history_name NAME`` to the command. You should be able to see the workflow
executing in the web browser, if you navigate to the 'List all histories' view. 
If you prefer to download data without interacting with the web interface at all,
you can add ``--download_outputs --output_directory . --output_json output.json``
to the command as before.

Typing ``--engine external_galaxy --galaxy_url SERVER_URL --galaxy_user_key YOUR_API_KEY``
each time you want to execute a workflow is a bit annoying. Fortunately, Planemo
provides the option to create 'profiles' which save this information for you.
To create a new profile called ``tutorial_profile``, you can run a command like
the following:

::

    $ planemo profile_create tutorial_profile --galaxy_url SERVER_URL --galaxy_user_key YOUR_API_KEY --engine external_galaxy
    Profile [tutorial_profile] created.


This allows creation of multiple profiles (e.g. for different Galaxy servers).
A list of all created profiles is provided by the ``profile_list`` subcommand:

::

    $ planemo profile_list
    Looking for profiles...
    tutorial_profile
    usegalaxy-eu
    usegalaxy-org
    3 configured profiles are available.


Once the new ``tutorial_profile`` is created, a workflow can be executed with:

::

    $ planemo run tutorial.ga tutorial-job.yml --profile tutorial_profile


Generating the job file
===============================================

The example workflow used so far provides not only the workflow, but also the
job file which specifies the inputs to be used. If you have created and downloaded
your own workflow, you need to create this job file yourself. As a first step,
ensure that your workflow is linted correctly:

::

    $ planemo workflow_lint tutorial.ga
    Applying linter tests... CHECK
    .. CHECK: Tests appear structurally correct


In this case, linting completes successfully, but you might see a message such as
``WARNING: Workflow contained output without a label`` or ``WARNING: Test referenced
File path not found``.

To generate the job file, you can now run:

::

    $ planemo workflow_job_init tutorial.ga -o tutorial-init-job.yml


This generates a template for the job file which you can modify yourself. Opening
``tutorial-init-job.yml`` should show the following:

::

    $ cat tutorial-init-job.yml
    Dataset 1:
      class: File
      path: todo_test_data_path.ext
    Dataset 2:
      class: File
      path: todo_test_data_path.ext
    Number of lines: todo_param_value


For each of the specified inputs in the workflow, an entry is created in the
output YAML file. The two dataset inputs are classified as ``class: File``,
with a placeholder path included, which you should change to the paths of your
chosen input files. You can also specify the URL of a file available online,
by replacing the ``path`` attribute with ``location`` (e.g. ``location: https://website.org/file.txt``).
The placeholder value for the ``Number of lines`` parameter should also be replaced,
ensuring it is of the correct type, i.e. in this case an integer.

Another more complex example, also including a collection as input, might look
like the following:

::

    input_dataset:
      class: File
      path: todo_test_data_path.ext
    input_collection:
      class: Collection
      collection_type: list
      elements:
      - class: File
        identifier: todo_element_name
        path: todo_test_data_path.ext
      - class: File
        identifier: todo_element_name
        path: todo_test_data_path.ext
    input_parameter: todo_param_value


For the collection, each dataset is listed, with a path and identifier specified.

If you are creating a workflow for the first time, you should include tests to
ensure it works in the way intended. These tests can be run using the ``planemo test``,
command, just like Galaxy tool testing (for more information, see `here <best_practices_workflows.html#tests>`_).
These tests require a test file, similar to the job file used so far, which also
specifies expected outputs which can be used to validate the workflow. An
equivalent planemo command for creating a template for these test files is also
available:

::

    $ planemo workflow_test_init tutorial.ga -o tutorial-init-test.yml
    $ cat tutorial-init-test.yml
    - doc: Test outline for tutorial.ga
      job:
        Dataset 1:
          class: File
          path: todo_test_data_path.ext
        Dataset 2:
          class: File
          path: todo_test_data_path.ext
        Number of lines: todo_param_value
      outputs:
        output:
          class: ''


Using workflow and dataset IDs
===============================================

If you ran all the commands above then you probably noticed that both the
workflow and the input datasets get newly uploaded at each execution. If you
want to run the same workflow multiple times, you may prefer to avoid this.
In the examples given so far, all workflows and datasets are specified by means
of a local path, but Planemo also allows you to use the IDs created by Galaxy
as well. These IDs are unique to each Galaxy server, so this approach isn't
transferrable if you want to run your workflows on multiple servers.

The first step is to ensure all the datasets which are required for the
workflow are already uploaded. You can either do this by running the workflow
once in the normal way, as described above, or just manually uploading through
the web interface.

To get dataset IDs, you can click on the dataset's 'View details' button (a
small letter 'i' in a circle). This provides various information about the
dataset and the job which created it. Under the 'Job information' section,
there is a row named 'History Content API ID'. For each input dataset, copy
this string (it will probably look something like ``457d46215431cc37baf96108ad87f351``)
and paste it into the workflow job file so it looks something like the following:

::

    Dataset 1:
      class: File
      galaxy_id: "457d46215431cc37baf96108ad87f351"
    Dataset 2:
      class: File
      galaxy_id: "55f30adf41ae36455431abeaa185ed89"
    Number of lines: 3


i.e. just replace the ``path`` line with ``galaxy_id``.

You can do exactly the same with a collection; either of the following will
work:

::

    input_collection1:
      class: Collection
      galaxy_id: "9d362c51f575db89"
    input_collection2:
      class: Collection
      collection_type: list
      elements:
      - class: File
        identifier: element 1
        galaxy_id: "457d46215431cc37baf96108ad87f351"


For ``input_collection1``, an existing collection will be used (by specifying its
collection ID), whereas for ``input_collection2``, a new collection will be created
from a list of existing datasets.

Once the job file has been modified, run ``planemo run`` as before. The result
should be the same, though it should be a bit faster, since the upload step was
skipped. Instead, the selected datasets get copied to a new history, which
unlike a new upload, doesn't result in any additional storage being used.

To run the workflow using a workflow ID, replace the workflow file path with
the workflow ID from the Galaxy server:

::

    $ planemo run 501da2f0ba775fd0 tutorial-job.yml --profile tutorial_profile


Using aliases
===============================================

Once you are dealing with a large number of workflows and datasets, you may
find that it becomes difficult to keep track of the file paths or IDs
which you are using for execution, particularly if you are executing workflows
based on their ID. Planemo offers the option to create aliases, or easily
memorable mnemonics, for Galaxy workflows, with the following command:

::

    $ planemo create_alias 501da2f0ba775fd0 --alias my_favorite_workflow --profile tutorial_profile


You can then execute the workflow with:

::

    $ planemo run my_favorite_workflow tutorial-job.yml --profile tutorial_profile


Note that aliases are associated with a particular profile, so if you want to
execute the same workflow with multiple profiles, you should recreate the alias
for each one. Aliases can be created either for workflow IDs (as above) or for
workflow file paths. You can list all aliases associated with a profile with:

::

    $ planemo list_alias --profile tutorial_profile


Checking invocations
===============================================

Assuming you know the workflow ID (or an alias for it), you can get a list of
all created invocations with:

::

    $ planemo list_invocations my_favorite_workflow --profile tutorial_profile


This indicates the number of datasets created, as well as the state they are in
(running, errored, paused, etc.)


Profile configuration files
===============================================

Information about each of the files is located in a configuration file, located
at ``~/.planemo/profiles/{profile_name}/planemo_profile_options.json``.

If you ran all the commands in this tutorial, the contents should be similar to
the following:

::

    $ cat ~/.planemo/profiles/tutorial_profile/planemo_profile_options.json
    {
      "galaxy_url": "SERVER_URL",
      "galaxy_user_key": "YOUR_API_KEY",
      "galaxy_admin_key": null,
      "engine": "external_galaxy",
      "aliases": {
        "my_favorite_workflow": "501da2f0ba775fd0"
      }
    }


You can also delete unwanted profiles or aliases with these commands:

::

    $ planemo delete_alias --alias my_favorite_workflow --profile tutorial_profile
    $ planemo profile_delete tutorial_profile


Rerunning failed jobs
===============================================

A frequent issue that arises when running a complex workflow is that component
jobs can fail, resulting in failure of the entire workflow. These jobs can be
rerun in the graphical interface, selecting the ``Resume dependencies from this
job ?`` option, which restarts the paused workflow (so-called 'remapping' of the
failed job over the previously created output dataset(s)). However, if there are
a large number of failures, and you believe the errors are transitory, e.g. due
to some temporary server issue, you can rerun several failed jobs simultaneously
using the ``planemo rerun`` command:

::

    $ planemo rerun --history 68008488b4fb94de
    $ planemo rerun --invocation 27267240b7d1f22a a9b086729787c907c
    $ planemo rerun --job a2b39deaa34509bb 3318707f2f0ff1fd

In the first two cases, all failed, remappable jobs which are associated with
the specified history(s) or invocation(s) will be rerun. In the third case,
the specified jobs will simply be rerun.
