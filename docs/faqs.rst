=============================
FAQs
=============================

Below are a list of frequently asked questions (and answers) for users new to Planemo.

What is Planemo?
=============================

Planemo is a command-line application which is primarily used for developing Galaxy and CWL tools,
workflows and training materials. It's also used to deploying the tools and workflows once developed
(tools to the Galaxy ToolShed; workflows to the Intergalactic Workflow Commission) and executing
Galaxy-based analyses via the command-line, if you prefer or need to use that rather than Galaxy's
graphical interface.

Is Planemo the right tool for me?
=================================

Yes, if you are comfortable using the command-line and one of the following is true:
  - You want to use a piece of scientific software not currently available in Galaxy and therefore need to wrap it yourself.
  - You want to develop a workflow from individual Galaxy tools and submit it to a community repository such as the Intergalactic Workflow Commission.
  - You want to develop CWL tools or workflows.
  - You're writing a tutorial for the `Galaxy Training Network`_ and would like to use a pre-prepared template.
  - You want to run Galaxy workflows from the command line.

What is the difference between Galaxy and CWL? Which should I use?
==================================================================

Galaxy is a web-based scientific analysis platform which provides access to scientific software
via a graphical interface. It is also a scientific workflow management system and allows chaining
multiple individual tools together to form a complex pipeline, which can be executed as a single
tool. If you place a high value on a graphical interface, usability and developing flexible, reusable
tools which can easily be adapted for other purposes, Galaxy could be a good choice for you.

CWL (Common Workflow Language) is a workflow specification, rather than a workflow management system.
It aims to allow the description of scientific pipelines in a way that is independent of the particular
workflow manager used; if this is important to you, using CWL would be an excellent idea. CWL workflows
are supported by a variety of workflow managers, including Toil, cwltool, Arvados and Galaxy. 

Both Galaxy and CWL place a high value on community, open-source code and FAIR data analysis.

Can I develop Galaxy tools without Planemo?
===========================================

Yes, you can, but it is more complicated to do so and will take more time. Planemo encourages best practices in
software development and follows guidelines agreed on by the Galaxy community, particularly
test-driven development.

Where can I learn to write tools?
=================================

We have a `tutorial`_ in this documentation and `another <https://training.galaxyproject.org/training-material/topics/dev/tutorials/tool-from-scratch/tutorial.html>`__ incorporated into the
Galaxy Training Network, which should get you off to a good start.

If you use the Visual Studio Code editor, we recommend using Planemo in conjunction with the
`Galaxy Language Server`_, which encourages development best practices and should greatly
improve the speed at which you can develop new Galaxy tools.

How do I set up deployment to the ToolShed with Planemo?
========================================================

We recommend using continuous integration (e.g. GitHub Actions) to achieve this and provide
a `template`_ which you can use to set up a GitHub repository which automatically deploys tools
to the ToolShed when a pull request is merged.

If you don't want to handle tool deployment yourself, you can also submit your tool wrappers to
a community repository such as the `IUC`_. This also has the advantage that your wrappers will be
reviewed by experienced Galaxy users and developers.

How can I automate workflow execution with Planemo?
===================================================

Have a look at this `section`_ of the documentation, or this `training`_ provided by the Galaxy
Training Network.

How can I contribute to the project?
====================================

We would love to see new contributions to Planemo! Please have a look `here`_ to see the different
ways you can help.


.. _tutorial: https://planemo.readthedocs.io/en/latest/writing_standalone.html
.. _Galaxy Training Network: https://training.galaxyproject.org/
.. _Galaxy Language Server: https://github.com/galaxyproject/galaxy-language-server
.. _template: https://github.com/galaxyproject/galaxy-tool-repository-template
.. _IUC: https://github.com/galaxyproject/tools-iuc
.. _section: https://planemo.readthedocs.io/en/latest/running.html
.. _training: https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/workflow-automation/tutorial.html
.. _here: https://planemo.readthedocs.io/en/latest/contributing.html
