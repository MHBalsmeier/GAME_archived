#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

aim_dir=~/compiled/game

### END OF USUAL INPUT SECTION

echo installing grids ...
if [ -d $aim_dir/grids ]
then
rm -r $aim_dir/grids
fi
cp -r grid_generator/grids $aim_dir/grids
echo Grids installed.
