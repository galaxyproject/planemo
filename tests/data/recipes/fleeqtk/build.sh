#!/bin/bash

export C_INCLUDE_PATH=${PREFIX}/include
export LIBRARY_PATH=${PREFIX}/lib

make all
mkdir -p $PREFIX/bin
cp -f fleeqtk $PREFIX/bin/
