#!/bin/bash

# This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

if [ -d build ]
then
rm -r build
fi
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ../core
make
ctest
make install
cd ..
rm -r build
