layout: true
class: inverse, middle

---
class: special
# Planemo: A Scientific Workflow SDK
John Chilton, Aysam Guerler, and the Galaxy Team

---

### A Galaxy Philosophy

* The most important Galaxy user is the biologist using the GUI, they come first!
  * No one wants to inconvenience a bioinformatician or developer, but we will if absolutely required for the biologist.
* Galaxy workflows will never require an SDK.
  * Planemo, workflow formats, etc... are conveniences for people who prefer developer processes.

???

I'm going to get in trouble from a couple different sides with this talk, and I will state some strong, personal opinions. So I'm going to start with something that amounts to either a justification or a reassurance.

I'm pretty confident this amounts to a Galaxy team perspective - regardless
of how your define that team.

READ SLIDES

The first part of this talk will cover workflow enhancements from an end
user perspective in some ways though.

---

### GUI Enhancements - Workflow Editor Form

![Workflow Editor](images/gx_new_workflow_editor.png)

???

The editor form now uses the same backbone driven MVC components as the
new tool form presented last year.

---

### GUI Enhancements - Workflow Run Form

TODO: UPDATE IMAGE

![Workflow Run Form](images/gx_new_workflow_run.png)

???

The run workflow form has likewise been overhauled and will be merged soon. This
should allow more dynamic option control when running workflows.

---

### GUI Enhancements - Labels

TODO: UPDATE IMAGE

![Workflow Run Form](images/gx_new_workflow_run.png)

???

When reasoning about workflows and connections between steps, persistent and unique 
labels for steps and outputs are important. These are useful in the API and too a 
lesser extent in the GUI today.

A major theme of this presentation is going to be that workflows are programs, they
are a coding artifact. I'm not sure anyone would disagree with me on that - but I 
think the implications may be counter-intuitive at times.

---

### GUI Enhancements - Nested Workflows

TODO: UPDATE IMAGE

![Nested Workflows](images/gx_new_workflow_run.png)

???

Workflows are programs, languages describing programs should provide abstractions for
composition. Nesting workflows was one of the most requested feature requests of 
Galaxy and it now supports this.

---

### What About Planemo?

???

Enough screenshots right - when I present people expect to see long command-lines!

---

### Planemo's Success

It is *the way* to develop Galaxy tools in 2016! Why?

- Artifact-centric - not Galaxy-centric or registry-centric.
- Fits into existing developer workflows - CLI, Git(hub), CI (Travis).
- Flexible
- Well documented w/focus on usage examples.

It is about **developer workflow**.

???

In 2014, Greg von Kuster presented a tool development workflow that involved 
publishing things to a local tool shed and running tests from there and viewing
the results through the web interface. I call this registry or shed centric tool
development (development activities "boot strap the tool shed, upload to the tool 
shed, run tests against the tool shed, view results in the tool shed, export capsule
from the tool shed". Prior to that tool development was Galaxy-centric - "download 
Galaxy, update the Galaxy tool conf, update the Galaxy test data, run the Galaxy 
tests." Planemo is tool centric - lint the tool, test the tool, serve the tool.

---

class: bottom
background-image: url(images/organic_mower_wat.jpg)
background-position: center;
background-repeat: no-repeat;
background-size: contain;

### So Planemo?

Photo Credit: Peter Smith (@skwiot)

---

### Workflows are Programs

When I write programs...

* ... I write tests (and write them first)!
* ... I store them on Github!
* ... use a text editor - my text editor!

???

The hippest, artisinal, free range text editor of my choice. It is a go 
lang rewrite of a node rewrite of an Erlang editor from 1987 - it is super hot right
now but I'm sure you've never heard of.

---

### Workflow Operations with Planemo

Serve them:

```
$ planemo serve <workflow.ga>
```

Brings up Galaxy interface with the workflow loaded.

---

### Load up local tools!

Serve them with tools...

```
$ planemo serve --extra_tools <tool_dir> <workflow.ga>
```

Supply as many tools as you want.

---

### Workflows are harder to serve than tools

- Longer to setup input data
- Takes time to install shed tools
- ``sqlite`` database locks
- Lack of cluster access is more likely to be a problem.

---

### Enter Planemo Profiles

Profile - a persistent, named Galaxy configuration available for serving and testing
across invocations.

```
$ planemo profile_create <name>
$ planemo serve --profile <name> <workflow.ga>
```

---

### Planemo Profiles and Postgres

```
$ planemo profile_create --postgres <name>
$ planemo serve --profile <name> <workflow.ga>
```

Automatically provision a postgres database and configure this profile to use it. All
future ``serve``, ``run``, and ``test`` invocations using this profile will use this 
database.

Makes some reasonable guesses, but Postgres connection settings can be configured with
additional CLI options or in ``~/.planemo.yml``.

---

### Planemo Profiles and Clusters

```
$ planemo profile_create --job_conf <job_conf.xml> <name>
$ planemo serve --profile <name> <workflow.ga>
```

Setup a job configuration object for this profile and use it with all future
``serve``, ``run``, and ``test`` invocations using this profile.

---

### Planemo Docker Profiles

Or just...

```
$ planemo profile_create --engine_type docker_galaxy <name>
$ planemo serve --profile <name> <workflow.ga>
```

Leverage the omnipresent ``docker-galaxy-stable`` community project spearheaded by the 
omnipresent Björn Grüning to ``serve`` and ``run`` using Galaxy in a Docker container.
This container comes pre-configured with slurm and Postgres as available locally -
but offers lots of additionals goodies such as Condor, ProFTP, supervisor, and uwsgi.

https://github.com/bgruening/docker-galaxy-stable

---

### Galaxy Workflow Format

Doesn't Galaxy have a JSON workflow format?

```json
"tool_state": "{\"__page__\": 0, \"__rerun_remap_job_id__\": null,
\"input1\": \"null\", \"chromInfo\": \"\\\"/home/john/workspace/galaxy-central/tool-data/shared/ucsc/chrom/?.len\\\"\",
\"queries\": \"[{\\\"input2\\\": null, \\\"__index__\\\": 0}]\"}",
```

- Neither human writable, nor human readable.
- JSON doesn't allow comments.
  - One shouldn't have to describe a configuration file in JSON,
    let alone write a program in it.

---

### Format 2 Workflows - Example

TODO: SHOW EXAMPLE

---

### Format 2 Workflows

- Started as a way to build test workflows for Galaxy testing framework.
- All steps can be labeled, connections described by ID.
- ``gxformat2`` Python library for conversion and loading into Galaxy.
  - Pypi @ https://pypi.python.org/pypi/gxformat2/
  - Github @ https://github.com/jmchilton/gxformat2



