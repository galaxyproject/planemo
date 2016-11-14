"""Module parses R scripts and sends a yaml file to cmd_bioc_tool_init."""
import os


def read_rscript(path):
    """Read the rscript."""
    try:
        with open(os.path.expanduser(path), 'r') as f:
            rscript = f.readlines()
    except Exception as e:
        print(e)
    return rscript


def parse_rscript(script, example_command):
    """Parse script."""
    rscript = read_rscript(script)
    data = {}

    # Find libraries
    lib = Library(rscript)
    library_list = lib.find_library()
    data['library'] = library_list

    # Find inputs
    inputs = Input(rscript, example_command)
    input_list = inputs.find_inputs()
    data['inputs'] = input_list

    # Find outputs
    outputs = Output(rscript, example_command)
    output_list = outputs.find_outputs()
    data['outputs'] = output_list
    return data


def parse_example_command(example_command):
    """Parse example_command to get inputs.

    Each input stored as element in a dictionary list.
    """
    cmd = example_command.replace("\n", " ")
    opts = [i.strip() for i in cmd.split("--")]
    opt_dict = {}
    for opt in opts:
        opt = opt.split(" ")
        if not opt[0] in opt_dict.keys():
            opt_dict[opt[0]] = [opt[1]]
        else:
            opt_dict[opt[0]].append(opt[1])
    return opt_dict


class Library(object):
    """Library class for parsing R scripts."""

    def __init__(self, script):
        """Initialize class Library."""
        self.script = script
        self.searchtext = "library"

    def _prune_library(self, line):
        """Prune line to get the names in library."""
        import re
        split_words = re.compile('\w+').findall(line)
        lib = [w for w in split_words if w != "library"]
        return lib[0]

    def find_library(self):
        """Parse library, to find and check requirements."""
        lib = []
        for i, line in enumerate(self.script):
            line = line.strip()
            if (self.searchtext in line) and (not line.startswith("#")):
                lib_value = self._prune_library(line)
                # if lib_value != "getopt":  # getopt already exists
                lib.append(lib_value)
        return lib


class Input(object):
    """Input class for parsing inputs."""

    def __init__(self, script, example_command):
        """Initialize Input with searchtext = input."""
        self.script = script
        self.example_command = example_command
        self.searchtext = "input"

    def find_inputs(self):
        """Find inputs in example command.

        This parses the R script and has NOTHING TO DO WITH kwds
        """
        opt_dict = parse_example_command(self.example_command)
        inputs = {}
        for key, value in opt_dict.iteritems():
            if self.searchtext in key:  # key here is "input"
                for i, line in enumerate(self.script):
                    line = line.strip()
                    if (key in line) and (not line.startswith("#")):
                        # print >> sys.stderr, 'Line: %s\nkey: %s\nvalue: %s' % (line,key,value)
                        inputs[key] = value
                    else:
                        continue
        # print >> sys.stderr, 'INPUTS: %s' % inputs
        if not bool(inputs):  # if inputs are empty
            print("No inputs found in the Rscript, please specify inputs.")
        return inputs


class Output(object):
    """Output class for parsing outputs."""

    def __init__(self, script, example_command):
        """Initialize Input with searchtext - output."""
        self.script = script
        self.example_command = example_command
        self.searchtext = "output"

    def find_outputs(self):
        """Find outputs in example command."""
        opt_dict = parse_example_command(self.example_command)
        outputs = {}
        for key, value in opt_dict.iteritems():
            if self.searchtext in key:
                for i, line in enumerate(self.script):
                    line = line.strip()
                    if (key in line) and (not line.startswith("#")):
                        outputs[key] = value
                    else:
                        continue
        # if not bool(outputs):  # if outputs are empty
        #     print("No explicit outputs found, please specify outputs.")
        return outputs


if __name__ == "__main__":
    # TODO : Make sure tools with configfile are not used, this is not supported yet
    test_file1 = "/Users/nturaga/Documents/galaxyproject/bioc-galaxy-integration/my_r_tool/my_r_tool.R"
    test_file2 = "/Users/nturaga/Documents/galaxyproject/bioc-galaxy-integration/my_r_tool/my_r_tool_verbose.R"
    test_file3 = "/Users/nturaga/Documents/galaxyproject/bioc-galaxy-integration/my_r_tool/my_r_tool_multi_inputs_outputs.R"
    test_file4 = "/Users/nturaga/Documents/galaxyproject/bioc-galaxy-integration/my_r_tool/my_r_tool_fail_case.R"

    # Test case with explicit input and outputs
    print(" \n ===  Tool test 1 ==== \n ")
    parse_rscript(test_file1, "Rscript my_r_tool.R --input input.csv --output output.csv")

    # Test case with NO EXPLICIT OUTPUT
    print(" \n ===  Tool test 2 ==== \n")
    parse_rscript(test_file2, "Rscript my_r_tool_verbose.R --verbose TRUE --input intput.csv")

    # Test case with tool which has to fail
    print("\n == Tool test 4: Fail case ==== \n")
    parse_rscript(test_file4, "Rscript my_r_tool_fail_case.R")
