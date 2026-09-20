#!/bin/bash
wget https://raw.githubusercontent.com/galaxyproject/tools-iuc/main/tools/bedtools/test-data/unionBedGraphs1.bg -o unionBedGraphs1.bedgraph
wget https://raw.githubusercontent.com/galaxyproject/tools-iuc/main/tools/bedtools/test-data/unionBedGraphs2.bg -o unionBedGraphs2.bedgraph
wget https://raw.githubusercontent.com/galaxyproject/tools-iuc/main/tools/bedtools/test-data/unionBedGraphs3.bg -o unionBedGraphs3.bedgraph
unionBedGraphs -i unionBedGraphs1.bedgraph -i unionBedGraphs2.bedgraph -i unionBedGraphs3.bedgraph > output.bedgraph

# planemo tool_init --force --id union_bedgraphs --requirement bedtools --name 'Union Bedgraphs' --example_command 'unionBedGraphs -i unionBedGraphs1.bg -i unionBedGraphs2.bg -i unionBedGraphs3.bg > output.bedgraph' --example_input 'unionBedGraphs1.bg' --example_input 'unionBedGraphs2.bg' --example_input 'unionBedGraphs3.bg' --example_output 'output.bedgraph' --test_case --doi '10.1093/bioinformatics/btq033' --help_text 'Union BedGraphs'

# planemo t union_bedgraphs.xml ## Fails

# planemo t --conda_dependency_resolution  --conda_auto_install --conda_auto_init union_bedgraphs.xml ## works...
