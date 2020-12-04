Workflow execution against an external Galaxy
===============================================

The first requirement for executing workflows on an external Galaxy server is
a user account for that server. If you don't already have one, usegalaxy.org,
usegalaxy.eu and usegalaxy.org.au all provide free accounts which can be used
for this tutorial.

Assuming you have selected a server for this tutorial and have an account, you
need to retrieve the API key associated with that account. This can be found at
``{server_url}/user/api_key``, or by going to the User dropdown menu, selecting
'Preferences' and then clicking on 'Manage API key'.

Now you can run the workflow with:

::

    $ planemo run wf3-shed-tools.ga wf3-shed-tools-job.yml --engine external --url SERVER_URL --api_key YOUR_API_KEY


If you want to set the name of the history in which the workflow executes, add
``--history_name NAME`` to the command. You should be able to see the workflow
executing in the web browser, if you navigate to the 'List all histories' view. 

Typing ``--engine external --url SERVER_URL --api_key YOUR_API_KEY`` each time
you want to execute a workflow is a bit annoying. Fortunately, Planemo provides
the option to create 'profiles' which save this information for you. To create
a new profile called ``tutorial_profile``, you can run a command like the
following:

::

    $ planemo profile_create tutorial_profile --galaxy_url SERVER_URL --galaxy_user_key YOUR_API_KEY --engine external_galaxy



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

    $ planemo run wf3-shed-tools.ga wf3-shed-tools-job.yml --profile tutorial_profile


Generating the job file
===============================================

The example workflow used so far provides not only the workflow, but also the
job file which specifies the inputs to be used. If you have created and downloaded
your own workflow, you need to create this job file yourself. Planemo provides a
handy command which provides a template for this:

::

    $ planemo workflow_job_init wf3-shed-tools.ga -o test-job-file.yml


Opening ``test-job-file.yml`` should show the following:

::

    GSE37268_mof3:
      class: File
      path: todo_test_data_path.ext
    Genes:
      class: File
      path: todo_test_data_path.ext


For each of the specified inputs in the workflow, an entry is created in the
output YAML file. In this case, both inputs are datasets rather than collections
or parameters, so they are both classified as ``class: File``. A placeholder
path is also included, which you should change to the paths of your chosen
input files.

Another more complex example, also including a collection and parameter as
input, might look like the following:

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
    input_parameter: param_todo


For the collection, each dataset is listed, with a path and identifier specified.
For the parameter, only a single value has to be specified, but ensure it is of
the correct type, e.g. float, if that is what the workflow requires. 

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

    GSE37268_mof3:
      class: File
      galaxy_id: "457d46215431cc37baf96108ad87f351"
    Genes:
      class: File
      galaxy_id: "55f30adf41ae36455431abeaa185ed89"


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
        identifier: todo_element_name
        galaxy_id: "457d46215431cc37baf96108ad87f351"


For ``input_collection1``, an existing collection will be used (by specifying its
collection ID), whereas for ``input_collection2``, a new collection will be created from an
existing dataset.

Once the job file has been modified, run ``planemo run`` as before. The result
should be the same, though it should be a bit faster, since the upload step was
skipped. Instead, the selected datasets get copied to a new history, which
unlike a new upload, doesn't result in any additional storage being used.

To run the workflow using a workflow ID, replace the workflow file path with
the workflow ID from the Galaxy server:

::

    $ planemo run 501da2f0ba775fd0 wf3-shed-tools-job.yml --profile tutorial_profile


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

    $ planemo run my_favorite_workflow wf3-shed-tools-job.yml --profile tutorial_profile


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