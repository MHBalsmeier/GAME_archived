#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

echo "Creating directories for large files ..."
if [ ! -d grid_generator/grids ]
then
mkdir grid_generator/grids
fi
if [ ! -d grid_generator/statistics ]
then
mkdir grid_generator/statistics
fi
if [ ! -d surface_generator/surface_files ]
then
mkdir surface_generator/surface_files
fi
if [ ! -d test_generator/test_states ]
then
mkdir test_generator/test_states
fi
if [ ! -d output ]
then
mkdir output
fi
if [ ! -d nwp_init ]
then
mkdir nwp_init
fi
if [ ! -d figs ]
then
mkdir figs
fi
echo "Complete."
