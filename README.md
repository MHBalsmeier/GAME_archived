The Global Atmospheric Modeling Framework (GAME) is a global non-hydrostatic hexagonal C grid dynamical core with the possibility to advect a variable number of tracers. For the source terms, it makes use of the atmostracers library. For radiation, it is coupled to RTE+RRTMGP via the C binding rte-rrtmgp-c (not yet).

## Documents

The documentation can be found in the subdirectory doc. It contains a full description of the model with a focus on technical aspects (not yet). Scientific fundamentals are cited from the literature.

The handbook of the model can be found in the subdirectory handbook. It contains information on how to configure and run the model.

## Prerequisites

* atmostracers (https://github.com/MHBalsmeier/atmostracers)
* rte-rrtmgp-c (https://github.com/MHBalsmeier/rte-rrtmgp-c)
* geos95 (https://github.com/MHBalsmeier/geos95, only for grid generator)
* eccodes library (installation manual: https://mhbalsmeier.github.io/tutorials/eccodes_on_ubuntu.html)
* netcdf library
* CMake
* OpenMPI
* OpenMP

## Installing the model

### Download

### Build and install

For further information see the handbook.
