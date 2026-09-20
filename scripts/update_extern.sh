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
# Make sure planemo doesn't pick up the tool factory tests. (see github #161)
sed -i -e s/tests/tests_skip/ $TOOL_FACTORY_PATH/rgToolFactory2.xml
