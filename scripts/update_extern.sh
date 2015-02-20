#!/bin/bash

PLANEMO_SCRIPTS_DIR=`dirname $0`
PLANEMO_ROOT="$PLANEMO_SCRIPTS_DIR/.."

EXTERN_DIR="$PLANEMO_ROOT/planemo_ext"
cd "$PLANEMO_ROOT"


TOOL_FACTORY_PATH="$EXTERN_DIR/tool_factory_2"
mkdir -p $TOOL_FACTORY_PATH
for file in 'rgToolFactory2.xml' 'rgToolFactory2.py' 'getlocalrpackages.py' 'LICENSE';
do
    wget "https://raw.githubusercontent.com/galaxyproject/tools-iuc/master/tools/tool_factory_2/$file"  --output-document "$TOOL_FACTORY_PATH/$file"
done

XSD_PATH="$PLANEMO_ROOT/planemo/xml/xsd/tool"
mkdir -p $XSD_PATH
for file in 'galaxy.xsd' 'citations.xsd' 'citation.xsd' 'LICENSE';
do
    wget "https://raw.githubusercontent.com/JeanFred/Galaxy-XSD/master/$file" --output-document "$XSD_PATH/$file"
done
