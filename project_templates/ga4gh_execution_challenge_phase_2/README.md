# GA4GH Tool Execution Challenge - Phase 2.

This document describes how to use Planemo ([docs](http://planemo.readthedocs.io/en/latest/?badge=latest)
| [github](https://github.com/galaxyproject/planemo)) to enter Phase 2 of the
[GA4GH Tool Execution Challenge](https://www.synapse.org/#!Synapse:syn8507133/wiki/415976) using Galaxy.

1. Register for [Synapse account](https://www.synapse.org/#!RegisterAccount:0).
1. Become a [certified Synapse user](http://docs.synapse.org/articles/accounts_certified_users_and_profile_validation.html#certified-users) - requires passing a short quiz on data usage.
1. [Register for the challenge](https://www.synapse.org/#!Synapse:syn8507133/wiki/415976).
1. Setup a ``~/.synapseConfig`` file in your home directory.

   ```
   [authentication]
   username: me@example.com
   password: secret
   ```

1. Install Planemo and Synapse
   
   ```
   $ virtualenv planemo-venv; . planemo-venv/bin/activate
   $ pip install -U pip
   $ pip install planemo synapseclient
   ```

1. Pick a workflow, fetch it, run it, check it, and submit it (per-workflow instructions below).
1. Go to the resulting workflow entry page and update it with your specific environment settings and ensure it and the relevant files are shared. 

## The Workflows

### hello_world

https://www.synapse.org/#!Synapse:syn9630940

A CWL script is distributed to fetch the workflow from Synapse. We will run this
script using cwltool (installed with Planemo) - since it is a bit quicker than 
starting up a Galaxy server. We will use Galaxy only to run the actual workflow.

```
mkdir hello_world
cd hello_world
cp ~/.synapseConfig .
synapse get syn9733811
synapse get syn9732885
synapse get syn9770802
cwltool --non-strict dockstore-tool-synapse-get.cwl hello_world_get.cwl.json
```

Use Planemo to run the ``hello_world.cwl`` workflow that results.

```
planemo --verbose run --docker --non_strict_cwl --no_dependency_resolution --output_directory `pwd` hello_world.cwl hello_world.cwl.json
```

Check the results using the checker workflow:

```
cwltool --non-strict hello_world_checker.cwl hello_world_checker.cwl.json
```

Update ``hello_world_submit.cwl.json`` with a ``team_name`` field and a ``parent_id`` and then
submit using the submit worfklow. Create a new Project for the challenge (if you haven't 
already) and a folder under that project to create a ``parent_id`` for this submission. If you
wish to use the "GalaxyProject" team - contact John Chilton (jmchilton@gmail.com) to get added
to the project.

```
cwltool --non-strict dockstore-tool-synapse-submit.cwl hello_world_submit.cwl.json
```

Wait for an e-mail with additional instructions and then go to the Synapse entry and update
the fields in the wiki page generated for this new entry.

### md5sum

https://www.synapse.org/#!Synapse:syn10167920

A CWL script is distributed to fetch the workflow from Synapse. We will run this
script using cwltool (installed with Planemo) - since it is a bit quicker than 
starting up a Galaxy server. We will use Galaxy only to run the actual workflow.


```
mkdir md5sum
cd md5sum
cp ~/.synapseConfig .
synapse get syn9771740
synapse get syn9732885
synapse get syn9770802
cwltool --non-strict dockstore-tool-synapse-get.cwl md5sum_get.cwl.json
```

Use Planemo to run the ``hello_world.cwl`` workflow that results.

```
planemo --verbose run --docker --non_strict_cwl --no_dependency_resolution --output_directory `pwd` md5sum.cwl md5sum.cwl.json
```

Check the results with the checker workflow:

```
cwltool md5sum_checker.cwl md5sum_checker.cwl.json
cat results.json | grep -i overall
```

Update ``md5sum_submit.cwl.json`` with a ``team_name`` field and a ``parent_id`` and then
submit using the submit worfklow. Create a new Project for the challenge (if you haven't 
already) and a folder under that project to create a ``parent_id`` for this submission. If you
wish to use the "GalaxyProject" team - contact John Chilton (jmchilton@gmail.com) to get added
to the project.

```
cwltool --non-strict dockstore-tool-synapse-submit.cwl md5sum_submit.cwl.json
```

Wait for an e-mail with additional instructions and then go to the Synapse entry and update
the fields in the wiki page generated for this new entry.

### encode_mapping_workflow

```
mkdir encode_mapping
cd encode_mapping
cp ~/.synapseConfig .
synapse get syn9732885
synapse get -r syn10163025
```

### gdc_dnaseq_transform (WIP)

This worklfow doesn't work yet - see open issues at https://github.com/common-workflow-language/galaxy.

```
mkdir gdc_dnaseq_transform
cd gdc_dnaseq_transform
cp ~/.synapseConfig .
synapse get syn9962438
synapse get syn9732885
synapse get syn9770802 
cwltool --non-strict dockstore-tool-synapse-get.cwl gdc_dnaseq_transform_get.cwl.json
```

This last step downloads around 12 GB of data so this may take some time.

```
planemo --verbose run --docker --cwl_galaxy_root ~/workspace/galaxy --output_directory `pwd` cwl/workflows/dnaseq/transform.cwl transform.cwl.json
```

### biowardrobe_chipseq_se (WIP)

This worklfow doesn't work yet - see open issues at https://github.com/common-workflow-language/galaxy.

mkdir biowardrobe_chipseq_se
cd biowardrobe_chipseq_se
cp ~/.synapseConfig .
synapse get syn10204916
synapse get syn9770802
