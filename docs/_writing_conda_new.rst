
----------------------------------------------------------------
Building New Conda Packages
----------------------------------------------------------------

Frequently packages your tool will require are not found in Bioconda_
or conda-forge yet. In these cases, it is likely best to contribute
your package to one of these projects. Unless the tool is exceedingly
general Bioconda_ is usually the correct starting point.

.. note:: Many things that are not strictly or even remotely "bio" have
    been accepted into Bioconda_ - including tools for image analysis,
    natural language processing, and cheminformatics.

To get quickly learn to write Conda_ recipes for typical Galaxy tools,
please read the following pieces of external documentation.

- `Contributing to Bioconda <https://bioconda.github.io/contributing.html>`__ in particular focusing on

  - `One time setup <https://bioconda.github.io/contrib-setup.html>`__
  - `Contributing a recipe <https://bioconda.github.io/contribute-a-recipe.html>`__ (through "Write a Recipe")
- `Building conda packages <https://conda.io/docs/building/bpp.html#>`__ in particular

  - `Building conda packages with conda skeleton <https://conda.io/docs/build_tutorials/pkgs.html>`__ (the best approach for common scripting languages such as R and Python)
  - `Building conda packages from scratch <https://conda.io/docs/build_tutorials/pkgs2.html>`__
  - `Building conda packages for general code projects <https://conda.io/docs/build_tutorials/postgis.html>`__
  - `Using conda build <https://conda.io/docs/building/recipe.html>`__
- Then return to the Bioconda documentation and read

  - The rest of "Contributing a recipe" continuing from `Testing locally <https://bioconda.github.io/contribute-a-recipe.html#test-locally>`__
  - And finally `Guidelines for bioconda recipes <https://bioconda.github.io/guidelines.html>`__

These guidelines in particular can be skimmed depending on your recipe type, for
instance that document provides specific advice for:

- `Python <https://bioconda.github.io/guidelines.html#python>`__
- `R (CRAN) <https://bioconda.github.io/guidelines.html#r-cran>`__
- `R (Bioconductor) <https://bioconda.github.io/guidelines.html#r-bioconductor>`__
- `Perl <https://bioconda.github.io/guidelines.html#perl>`__
- `C/C++ <https://bioconda.github.io/guidelines.html#c-c>`__

To go a little deeper, you may want to read: 

- `Specification for meta.yaml <https://conda.io/docs/building/meta-yaml.html>`__
- `Environment variables <https://conda.io/docs/building/environment-vars.html>`__
- `Custom channels <https://conda.io/docs/custom-channels.html>`__

And finally to debug problems the `Bioconda troubleshooting <https://bioconda.github.io/troubleshooting.html>`__
documentation may prove useful.

----------------------------------------------------------------
Exercise - Build a Recipe
----------------------------------------------------------------

If you have just completed the exercise above - this exercise can be found in parent folder. Get
there with ``cd ../exercise2``. If not, the exercise can be downloaded with
