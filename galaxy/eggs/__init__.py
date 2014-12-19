PLANEMO_WARNING = ("It looks like you are trying to use Galaxy from within a "
                   "Python environment that has planemo installed. Please "
                   "either deactivate planemo's virtual environment, set up "
                   "one for Galaxy, or uninstall planemo from this "
                   "environment.")

raise ImportError(PLANEMO_WARNING)
