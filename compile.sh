#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/game

if [ ! -d build ]
then
mkdir build
fi

cd build

cmake -DCMAKE_INSTALL_PREFIX=../output ..
make

cd ..
