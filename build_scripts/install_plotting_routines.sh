#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

aim_dir=~/compiled/game

### END OF USUAL INPUT SECTION

echo installing plotting routines ...
if [ -d $aim_dir/plotting ]
then
rm -r $aim_dir/plotting
fi
cp -r ../plotting $aim_dir/plotting
echo Plotting routines installed.
