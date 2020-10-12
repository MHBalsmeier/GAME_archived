#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

aim_dir=~/compiled/game_dev

### END OF USUAL INPUT SECTION

echo installing output ...
if [ -d $aim_dir/output ]
then
rm -r $aim_dir/output
fi
mkdir $aim_dir/output
echo output installed.
