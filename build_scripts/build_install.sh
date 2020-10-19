#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

aim_dir=~/compiled/game

### END OF USUAL INPUT SECTION

if [ -d $aim_dir/core ]
then
rm -r $aim_dir/core
fi
if [ -d ../build ]
then
rm -r ../build
fi
mkdir ../build && cd ../build
cmake -DCMAKE_INSTALL_PREFIX=$aim_dir ../core
make
ctest
make install
chmod +x $aim_dir/core/run.sh
cd - > /dev/null
if [ ! -d $aim_dir/output ]
then
mkdir $aim_dir/output
fi
