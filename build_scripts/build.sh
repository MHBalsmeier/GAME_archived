#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

aim_dir=~/compiled/game

### END OF USUAL INPUT SECTION

rm -r $aim_dir/core
rm -r $aim_dir/input
rm -r $aim_dir/grids
if [ -d build ]
then
rm -r build
fi
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$aim_dir ../core
make
ctest
cd ..
