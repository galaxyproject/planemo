#!/bin/bash

usage() {
cat << EOF
Usage: ${0##*/} [-i] /path/to/galaxy...
Sync Planemo shared modules to those same modules in Galaxy directory (or vice versa if -i).

EOF
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

invert=0
OPTIND=1
while getopts ":i" opt; do
    case "$opt" in
        h)
            usage
            exit 0
            ;;
        i)
            invert=1
            ;;
        '?')
            usage >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

PLANEMO_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GALAXY_DIRECTORY=$1
GALAXY_LIB_DIR=$GALAXY_DIRECTORY/lib/galaxy

if [ "$invert" -ne "1" ];
then

    rm -rf $GALAXY_LIB_DIR/objectstore 
    cp -r $PLANEMO_DIRECTORY/planemo_ext/galaxy/objectstore $GALAXY_LIB_DIR

    rm -rf $GALAXY_LIB_DIR/tools/deps
    cp -r $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools/deps $GALAXY_LIB_DIR/tools

    rm -rf $GALAXY_LIB_DIR/jobs/metrics
    cp -r $PLANEMO_DIRECTORY/planemo_ext/galaxy/jobs/metrics $GALAXY_LIB_DIR/jobs

    rm -rf $GALAXY_LIB_DIR/tools/linters
    cp -r $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools/linters $GALAXY_LIB_DIR/tools

    cp $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools/lint.py $GALAXY_LIB_DIR/tools/lint.py
    cp $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools/loader.py $GALAXY_LIB_DIR/tools/loader.py
    cp $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools/loader_directory.py $GALAXY_LIB_DIR/tools/loader_directory.py

    cp $PLANEMO_DIRECTORY/planemo_ext/galaxy/util/plugin_config.py $GALAXY_LIB_DIR/util
    cp $PLANEMO_DIRECTORY/planemo_ext/galaxy/util/xml_macros.py $GALAXY_LIB_DIR/util

else

    rm -rf $PLANEMO_DIRECTORY/planemo_ext/galaxy/objectstore
    cp -r $GALAXY_LIB_DIR/objectstore $PLANEMO_DIRECTORY/planemo_ext/galaxy

    rm -rf $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools/deps
    cp -r $GALAXY_LIB_DIR/tools/deps $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools

    rm -rf $PLANEMO_DIRECTORY/planemo_ext/galaxy/jobs/metrics
    cp -r $GALAXY_LIB_DIR/jobs/metrics $PLANEMO_DIRECTORY/planemo_ext/galaxy/jobs

    rm -rf $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools/linters
    cp -r $GALAXY_LIB_DIR/tools/linters $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools

    cp $GALAXY_LIB_DIR/tools/lint.py $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools
    cp $GALAXY_LIB_DIR/tools/loader.py $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools
    cp $GALAXY_LIB_DIR/tools/loader_directory.py $PLANEMO_DIRECTORY/planemo_ext/galaxy/tools
    cp $GALAXY_LIB_DIR/util/plugin_config.py $PLANEMO_DIRECTORY/planemo_ext/galaxy/util/
    cp $GALAXY_LIB_DIR/util/xml_macros.py $PLANEMO_DIRECTORY/planemo_ext/galaxy/util/

fi
