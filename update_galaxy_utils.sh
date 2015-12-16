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
GALAXY_LIB_DIR=$GALAXY_DIRECTORY/lib

UTIL_FILES=(__init__.py aliaspickler.py bunch.py directory_hash.py docutils_template.txt expressions.py heartbeat.py heartbeat.py inflection.py json.py jstree.py lazy_process.py odict.py object_wrapper.py plugin_config.py properties.py simplegraph.py sleeper.py sockets.py specs.py sqlite.py submodules.py topsort.py topsort.py xml_macros.py)

# Set of files in Galaxy (along with UTIL_FILES above that have minimal dependencies,
# and are Python 2/3 compat.)
GALAXY_LIB=(galaxy/objectstore galaxy/tools/deps galaxy/tools/parser galaxy/jobs/metrics galaxy/tools/linters galaxy/tools/loader_directory.py galaxy/tools/loader.py galaxy/tools/lint.py galaxy/tools/lint_util.py galaxy/tools/deps galaxy/tools/toolbox)

if [ "$invert" -ne "1" ];
then

    for f in "${UTIL_FILES[@]}"
    do
        cp $PLANEMO_DIRECTORY/planemo_ext/galaxy/util/$f $GALAXY_LIB_DIR/galaxy/util
    done

    for f in "${GALAXY_LIB[@]}"
    do
        rm -rf $GALAXY_LIB_DIR/$f
        cp -r $PLANEMO_DIRECTORY/planemo_ext/$f $GALAXY_LIB_DIR/$f
    done

else

    rm -rf $PLANEMO_DIRECTORY/planemo_ext/galaxy/util/*
    for f in "${UTIL_FILES[@]}"
    do
        cp -r $GALAXY_LIB_DIR/galaxy/util/$f $PLANEMO_DIRECTORY/planemo_ext/galaxy/util
    done

    for f in "${GALAXY_LIB[@]}"
    do
        rm -rf $PLANEMO_DIRECTORY/planemo_ext/$f
        cp -r $GALAXY_LIB_DIR/$f $PLANEMO_DIRECTORY/planemo_ext/$f
    done

fi
