#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

aim_dir=~/compiled/game_dev

### END OF USUAL INPUT SECTION

echo installing test initialization states ...
if [ -d $aim_dir/input ]
then
rm -r $aim_dir/input
fi
cp -r ../test_generator/test_states $aim_dir/input
echo Test initialization states installed.
