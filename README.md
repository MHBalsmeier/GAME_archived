The Global Atmospheric Modeling Framework (GAME) is a global non-hydrostatic hexagonal C grid dynamical core with the possibility to advect a variable number of tracers. For the source terms, it makes use of the atmostracers library. For radiation, it is coupled to RTE+RRTMGP via the C binding rte-rrtmgp-c (to be done).

## Documents

The documentation can be found in the subdirectory doc. It contains an overview of numerics and references to the literature as well as a description of the code.

The handbook of the model can be found in the subdirectory handbook. It contains information on how to generate necessary files like grid files and input data and explains how to compile, configure and run the model.

## Prerequisites

* geos95 (https://github.com/MHBalsmeier/geos95, only for grid generator)
* netcdf library
* OpenMPI
* OpenMP
* eccodes library (installation manual: https://mhbalsmeier.github.io/tutorials/eccodes_on_ubuntu.html)
* CMake
* atmostracers (https://github.com/MHBalsmeier/atmostracers)
* rte-rrtmgp-c (https://github.com/MHBalsmeier/rte-rrtmgp-c)

## Installing the model



### Download



### Build and install






