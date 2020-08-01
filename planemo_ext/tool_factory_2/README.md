*WARNING before you start*

 Install this tool on a private Galaxy ONLY
 Please NEVER on a public or production instance
 
Updated august 2014 by John Chilton adding citation support

Updated august 8 2014 to fix bugs reported by Marius van den Beek

Please cite the resource at
http://bioinformatics.oxfordjournals.org/cgi/reprint/bts573?ijkey=lczQh1sWrMwdYWJ&keytype=ref
if you use this tool in your published work.

**Short Story**

This is an unusual Galaxy tool capable of generating new Galaxy tools.
It works by exposing *unrestricted* and therefore extremely dangerous scripting
to all designated administrators of the host Galaxy server, allowing them to
run scripts in R, python, sh and perl over multiple selected input data sets,
writing a single new data set as output.

*You have a working r/python/perl/bash script or any executable with positional or argparse style parameters*

It can be turned into an ordinary Galaxy tool in minutes, using a Galaxy tool.


**Automated generation of new Galaxy tools for installation into any Galaxy**

A test is generated using small sample test data inputs and parameter settings you supply.
Once the test case outputs have been produced, they can be used to build a
new Galaxy tool. The supplied script or executable is baked as a requirement
into a new, ordinary Galaxy tool, fully workflow compatible out of the box.
Generated tools are installed via a tool shed by an administrator
and work exactly like all other Galaxy tools for your users.

**More Detail**

To use the ToolFactory, you should have prepared a script to paste into a
text box, or have a package in mind and a small test input example ready to select from your history
to test your new script.

```planemo test rgToolFactory2.xml --galaxy_root ~/galaxy --test_data ~/galaxy/tools/tool_makers/toolfactory/test-data``` works for me

There is an example in each scripting language on the Tool Factory form. You
can just cut and paste these to try it out - remember to select the right
interpreter please. You'll also need to create a small test data set using
the Galaxy history add new data tool.

If the script fails somehow, use the "redo" button on the tool output in
your history to recreate the form complete with broken script. Fix the bug
and execute again. Rinse, wash, repeat.

Once the script runs sucessfully, a new Galaxy tool that runs your script
can be generated. Select the "generate" option and supply some help text and
names. The new tool will be generated in the form of a new Galaxy datatype
*toolshed.gz* - as the name suggests, it's an archive ready to upload to a
Galaxy ToolShed as a new tool repository.

Once it's in a ToolShed, it can be installed into any local Galaxy server
from the server administrative interface.

Once the new tool is installed, local users can run it - each time, the script
that was supplied when it was built will be executed with the input chosen
from the user's history. In other words, the tools you generate with the
ToolFactory run just like any other Galaxy tool,but run your script every time.

Tool factory tools are perfect for workflow components. One input, one output,
no variables.

*To fully and safely exploit the awesome power* of this tool,
Galaxy and the ToolShed, you should be a developer installing this
tool on a private/personal/scratch local instance where you are an
admin_user. Then, if you break it, you get to keep all the pieces see
https://bitbucket.org/fubar/galaxytoolfactory/wiki/Home

**Installation**
This is a Galaxy tool. You can install it most conveniently using the
administrative "Search and browse tool sheds" link. Find the Galaxy Main
toolshed at https://toolshed.g2.bx.psu.edu/ and search for the toolfactory
repository. Open it and review the code and select the option to install it.

If you can't get the tool that way, the xml and py files here need to be
copied into a new tools
subdirectory such as tools/toolfactory Your tool_conf.xml needs a new entry
pointing to the xml
file - something like::

  <section name="Tool building tools" id="toolbuilders">
    <tool file="toolfactory/rgToolFactory.xml"/>
  </section>

If not already there,
please add:
<datatype extension="toolshed.gz" type="galaxy.datatypes.binary:Binary"
mimetype="multipart/x-gzip" subclass="True" />
to your local data_types_conf.xml.


**Restricted execution**

The tool factory tool itself will then be usable ONLY by admin users -
people with IDs in admin_users in universe_wsgi.ini **Yes, that's right. ONLY
admin_users can run this tool** Think about it for a moment. If allowed to
run any arbitrary script on your Galaxy server, the only thing that would
impede a miscreant bent on destroying all your Galaxy data would probably
be lack of appropriate technical skills.

**What it does** 

This is a tool factory for simple scripts in python, R and
perl currently. Functional tests are automatically generated. How cool is that.

LIMITED to simple scripts that read one input from the history. Optionally can
write one new history dataset, and optionally collect any number of outputs
into links on an autogenerated HTML index page for the user to navigate -
useful if the script writes images and output files - pdf outputs are shown
as thumbnails and R's bloated pdf's are shrunk with ghostscript so that and
imagemagik need to be available.

Generated tools can be edited and enhanced like any Galaxy tool, so start
small and build up since a generated script gets you a serious leg up to a
more complex one.

**What you do**

You paste and run your script, you fix the syntax errors and
eventually it runs. You can use the redo button and edit the script before
trying to rerun it as you debug - it works pretty well.

Once the script works on some test data, you can generate a toolshed compatible
gzip file containing your script ready to run as an ordinary Galaxy tool in
a repository on your local toolshed. That means safe and largely automated
installation in any production Galaxy configured to use your toolshed.

**Generated tool Security**

Once you install a generated tool, it's just
another tool - assuming the script is safe. They just run normally and their
user cannot do anything unusually insecure but please, practice safe toolshed.
Read the code before you install any tool. Especially this one - it is really scary.

**Send Code**

Patches and suggestions welcome as bitbucket issues please?

**Attribution**

Creating re-usable tools from scripts: The Galaxy Tool Factory
Ross Lazarus; Antony Kaspi; Mark Ziemann; The Galaxy Team
Bioinformatics 2012; doi: 10.1093/bioinformatics/bts573

http://bioinformatics.oxfordjournals.org/cgi/reprint/bts573?ijkey=lczQh1sWrMwdYWJ&keytype=ref

**Licensing**

Copyright Ross Lazarus 2010
ross lazarus at g mail period com

All rights reserved.

Licensed under the LGPL

**Obligatory screenshot**

http://bitbucket.org/fubar/galaxytoolmaker/src/fda8032fe989/images/dynamicScriptTool.png

