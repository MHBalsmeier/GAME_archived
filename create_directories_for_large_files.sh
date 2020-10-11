#!/bin/bash

# This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "Creating directories for large files ..."
if [ ! -d grid_generator/grids ]
then
mkdir grid_generator/grids
fi
if [ ! -d orography_generator/orographies ]
then
mkdir orography_generator/orographies
fi
if [ ! -d test_generator/test_states ]
then
mkdir test_generator/test_states
fi
echo "Complete."
