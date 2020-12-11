The Geophysical Fluids Modeling Framework (GAME) is a non-hydrostatic hexagonal C grid dynamical core with the possibility to advect a variable number of constituents. For radiation, it is coupled to RTE+RRTMGP.

## Overview

It is known that the forecast skill of a NWP model depends more on physics and data assimilation than on the dynamical core. However, all dynamical cores I know of have inconsistencies even in the most fundamental dynamical quantities (mass, energy forms and entropy). That is why the aim of this project is to develop a next generation dynamical core with the following properties:

* Stability.
* The numerical dispersion relation shall contain no unphysical branches.
* The numerical dispersion relation shall contain all relevant physical branches, including a geostrophic mode.
* Strong scalability on massively parallel computer architectures.
* Mass conservation to machine precision.
* Energy conversions shall be based on a spatial discretization of Poisson brackets.
* Total energy conservation, apart from non-cancelling errors through explicit time stepping.
* Consistent local dissipation and entropy production.
* Satisfaction of the Second Law of Thermodynamics including entropy conservation to machine precision in an adiabatic setup.
* Absence of unphysical numerical stabilizers like divergence damping, fixers, filters and so on.
* No problems with terrain following coordinates. For example: A resting atmosphere around steep orography shall remain at rest.
* Ellipsoidal grid geometry.
* A capable and flexible framework for coupling to physics and to other components of an Earth system model.
* Consistency also in the presence of multiple constituents and radiation.

According to my understanding, a hexagonal C grid is the only discretization where all this can be achieved.

### GAME's principles a.k.a. why a new model is necessary

What GAME does what other models do not do and why:

* It uses the entropy as a prognostic variable. Usually, models use the potential temperature as a prognostic variable which is a conserved quantity and therefore the only forcings are the diabatic forcings rendering it a suitable variable for modeling. However, the same is true for the real entropy (connected to the density times the logarithm of the potential temperature), and this last quantitiy is the much more fundamental physical property.
* It employs the modified TRSK scheme.
* It can assign individual densities (instead of mixing ratios) to constituents as well as individual temperatures and sink velocities.
* It has different options for the time stepping as well as the complexity of the physics to make it useful for a wide range of applications.

### Things to be done

* Improvement of the momentum diffusion operator including dissipation.
* A regional mode.
* Implementation of SLEVE.
* Implementation of MPI parallelization.
* A way to construct Voronoi meshes on an ellipsoid.
* A largely implicit 3D solver for efficiency (larger time step) and better energy conservation properties.
* A nesting option.
* Implementation of ocean dynamics and physics.

## Documents

### Scientific derivations

The derivations of both the continuous equations as well as of the discretization techniques can be found in my textbook on theoretical meteorology (in German): [Kompendium Theoretische Meteorologie](https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf). The most fundamental numerical techniques have firstly been published in the literature cited below.

### Documentation

The documentation can be found in the subdirectory doc. It does not contain scientific derivations but rather describes the code.

### Handbook

The handbook of the model can be found in the subdirectory handbook. It is to be understood as a user manual and contains information on how to generate necessary files like grid files, geospatial data files (orography, for example) and input data and explains how to configure, compile and run the model.

## Installing the model

It is recommended to run the model on Linux. We will not help people who have problems trying to install the model on other operating systems.

### Dependencies

Everything is easy and quick to install.

* [geos95](https://github.com/MHBalsmeier/geos95)
* netcdf library (Ubuntu: sudo apt-get libnetcdf-dev, sudo apt-get libnetcdff-dev)
* eccodes library (installation manual: https://mhbalsmeier.github.io/tutorials/eccodes_on_ubuntu.html)
* CMake (Ubuntu: sudo apt-get install cmake)
* [atmostracers](https://github.com/MHBalsmeier/atmostracers)
* OpenMPI (Ubuntu: sudo apt-get install mpich)
* clone the DCMIP2016 repository: git clone https://github.com/ClimateGlobalChange/DCMIP2016.git
* clone the RTE+RRTMGP repository: git clone https://github.com/earth-system-radiation/rte-rrtmgp

#### For using the plotting routines

The following packages are additionally required if you want to make use of the plotting routines:

* Python and the visualization library scitools-iris (installation manual: https://mhbalsmeier.github.io/tutorials/iris_on_ubuntu.html)
* FFMPEG (Ubuntu: sudo apt-get install ffmpeg)

#### For developing

* Valgrind (Ubuntu: sudo apt-get install valgrind, for doing checks)

### Download

	git clone https://github.com/MHBalsmeier/game.git
	cd game
	./create_directories_for_large_files.sh

### Build and install

Modify the file core/CMakeLists.txt (read the comments). If you want to use radiation, modify rrtmgp_coefficients_file_sw and rrtmgp_coefficients_file_lw in the file core/src/physics/mo_radiation.f90.

In the shell scripts controlling the build process (residing in the directory build\_scripts) change the variable aim\_dir to a place of your choice, then run the scripts. The files with the suffix \_dev are meant to install to a location where new versions can be tested. You also need to install the run scripts in order to have the run scripts of the model where they belong. Install the plotting routines if you want to make use of them.

## Fundamental literature

* Staniforth, A. and Thuburn, J. (2012), Horizontal grids for global weather and climate prediction models: a review. Q.J.R. Meteorol. Soc., 138: 1-26. doi:10.1002/qj.958
* Thuburn, John. (2008). Numerical wave propagation on the hexagonal C-grid. Journal of Computational Physics. 227. 5836-5858. 10.1016/j.jcp.2008.02.010. 
* Gassmann, Almut & Herzog, Hans-Joachim. (2008). Towards a consistent numerical compressible non‐hydrostatic model using generalized Hamiltonian tools. Quarterly Journal of the Royal Meteorological Society. 134. 1597 - 1613. 10.1002/qj.297.
* Thuburn, John et al. “Numerical representation of geostrophic modes on arbitrarily structured C-grids.” J. Comput. Phys. 228 (2009): 8321-8335.
* Ringler, Todd & Thuburn, John & Klemp, J. & Skamarock, W.C.. (2010). A unified approach to energy conservation and potential vorticity dynamics on arbitrarily structured C-grids. J. Comput. Physics. 229. 3065-3090. 10.1016/j.jcp.2009.12.007.
* Gassmann, A. (2013), A global hexagonal C‐grid non‐hydrostatic dynamical core (ICON‐IAP) designed for energetic consistency. Q.J.R. Meteorol. Soc., 139: 152-175. doi:10.1002/qj.1960
* Gassmann, A. Discretization of generalized Coriolis and friction terms on the deformed hexagonal C‐grid. Q J RMeteorol Soc. 2018; 144: 2038– 2053. https://doi.org/10.1002/qj.3294



