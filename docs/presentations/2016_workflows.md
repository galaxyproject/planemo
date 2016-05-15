## Planemo: A Scientific Workflow SDK

---

### A Galaxy Philosophy

* #### The most important Galaxy user is the biologist using the GUI, they come first!
  * No one wants to inconvenience a bioinformatician or developer, but we will if absolutely required for the biologist.
* #### Galaxy workflows will never require an SDK.
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

---

### GUI Enhancements - Workflow Run Form

TODO: UPDATE IMAGE

![Workflow Run Form](images/gx_new_workflow_run.png)

???


---

### Workflows are Programs

When I write programs...

* ... I write tests (and write them first)!
* ... I store them on Github!
* ... use a text editor - my text editor!

???

The artisinal, free range, organic text editor of my choice. It is a go 
lang rewrite of a node rewrite of an Erlang editor from 1987 - it is the hottest editor you've never heard of.

---

### Galaxy Workflow Format

Doesn't Galaxy have a JSON workflow format?

```json
"tool_state": "{\"__page__\": 0, \"__rerun_remap_job_id__\": null,
\"input1\": \"null\", \"chromInfo\": \"\\\"/home/john/workspace/galaxy-central/tool-data/shared/ucsc/chrom/?.len\\\"\",
\"queries\": \"[{\\\"input2\\\": null, \\\"__index__\\\": 0}]\"}",
```

Neither human writable, nor human readable.

JSON doesn't allow comments - one shouldn't have to describe a
configuration file in JSON let alone write a program in it.

---

### Format 2 Workflows

TODO: SHOW EXAMPLE
