The Geophysical Fluids Modeling Framework (GAME) is a non-hydrostatic hexagonal C-grid dynamical core with the possibility to advect a variable number of constituents. For radiation, it is coupled to RTE+RRTMGP.

## Overview

It is known that the forecast skill of a NWP model depends more on physics and data assimilation than on the dynamical core. However, all dynamical cores I know of have inconsistencies even in the most fundamental dynamical quantities (mass, energy forms and entropy). That is why the aim of this project is to develop a next generation dynamical core with the following properties:

* stability
* the numerical dispersion relation shall contain no unphysical branches
* the numerical dispersion relation shall contain all relevant physical branches, including a geostrophic mode
* strong scalability on massively parallel computer architectures
* mass conservation to machine precision
* energy conversions shall be based on a spatial discretization of Poisson brackets
* energy conservation in the spatial discretization
* consistent friction and entropy production
* absence of unphysical numerical stabilizers like divergence damping, fixers and filters
* no stability problems with terrain following coordinates
* ellipsoidal grid geometry
* a capable and flexible framework for coupling to physics and to other components of an Earth system model
* consistency also in the presence of multiple constituents and radiation

According to my understanding, a hexagonal C-grid is the only discretization where all this can be achieved.

### GAME's principles a.k.a. why a new model is necessary

What GAME does what other models do not do and why:


* It employs the modified TRSK scheme.
* It can assign individual densities (instead of mixing ratios) to constituents as well as individual temperatures and sink velocities.
* It has different options for the complexity of the physics to make it useful for a wide range of applications.

### Things to be done

* improvement of the momentum diffusion operator
* a regional mode
* SLEVE
* MPI parallelization
* a way to construct Voronoi meshes on an ellipsoid
* a nesting option
* ocean dynamics and physics

## What I would do differently if I developed the model again

I only have my free time to do NWP model development, and this model architecture requires much more. Therefore I do not continue the work on this model. However, these are the two most important things I would change immediately if I had time or if I could develop the model again:

* I would use the potential temperature as a prognostic variable and not the entropy. This is because the potential temperature is not less fundamental than the entropy and has an easier Poisson bracket structure in an ideal gas.
* I would introduce a hydrostatic background state to achieve more stability in terrain following coordinates, which is necessary for realistic NWP scenarios. This diminishes the self-consistency properties of the Poisson brackets a bit, but this is true only for the spatial discretization. Through the decreased truncation errors in the horizontal gradients, however, the accuracy and possibly even the self-consistency will be increased by a hydrostatic background state when looking at the combination of spatial and temporal discretization.

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

* eccodes library (installation manual: https://gist.github.com/MHBalsmeier/a01ad4e07ecf467c90fad2ac7719844a)
* [geos95](https://github.com/OpenNWP/geos95)
* [atmostracers](https://github.com/OpenNWP/atmostracers)
* Clone the DCMIP2016 repository: git clone https://github.com/ClimateGlobalChange/DCMIP2016.git
* Clone the RTE+RRTMGP repository: `git clone https://github.com/earth-system-radiation/rte-rrtmgp`; then replace all `!$omp` with `!!$omp` in the subdirectories `rte`, `rrtmgp` as well as `extensions`. Use the netcdf files in `rrtmgp/data` from the `develop` branch.

#### For using the plotting routines

The following packages are additionally required if you want to make use of the plotting routines:

* Python and the visualization library scitools-iris (installation manual: https://scitools-iris.readthedocs.io/en/latest/installing.html#installing-from-source-without-conda-on-debian-based-linux-distros-developers)
* FFMPEG (Ubuntu: `sudo apt-get install ffmpeg`)

#### For developing

* Valgrind (Ubuntu: `sudo apt-get install valgrind`, for doing checks)

### Download

	git clone https://github.com/OpenNWP/game.git
	cd game
	./create_directories_for_large_files.sh

### Build and install

Modify the file core/CMakeLists.txt (read the comments). If you want to use radiation, modify rrtmgp_coefficients_file_sw and rrtmgp_coefficients_file_lw in the file core/src/radiation/rterrtmgp_coupler.f90.

In the shell scripts controlling the build process (residing in the directory build\_scripts) change the variable aim\_dir to a place of your choice, then run the scripts. The files with the suffix \_dev are meant to install to a location where new versions can be tested. You also need to install the run scripts in order to have the run scripts of the model where they belong. Install the plotting routines if you want to make use of them.

## Fundamental literature

* Staniforth, A. and Thuburn, J. (2012), Horizontal grids for global weather and climate prediction models: a review. Q.J.R. Meteorol. Soc., 138: 1-26. doi:10.1002/qj.958
* Thuburn, John. (2008). Numerical wave propagation on the hexagonal C-grid. Journal of Computational Physics. 227. 5836-5858. 10.1016/j.jcp.2008.02.010. 
* Gassmann, Almut & Herzog, Hans-Joachim. (2008). Towards a consistent numerical compressible non‐hydrostatic model using generalized Hamiltonian tools. Quarterly Journal of the Royal Meteorological Society. 134. 1597 - 1613. 10.1002/qj.297.
* Thuburn, John et al. “Numerical representation of geostrophic modes on arbitrarily structured C-grids.” J. Comput. Phys. 228 (2009): 8321-8335.
* Ringler, Todd & Thuburn, John & Klemp, J. & Skamarock, W.C.. (2010). A unified approach to energy conservation and potential vorticity dynamics on arbitrarily structured C-grids. J. Comput. Physics. 229. 3065-3090. 10.1016/j.jcp.2009.12.007.
* Gassmann, A. (2013), A global hexagonal C‐grid non‐hydrostatic dynamical core (ICON‐IAP) designed for energetic consistency. Q.J.R. Meteorol. Soc., 139: 152-175. doi:10.1002/qj.1960
* Gassmann, A. Discretization of generalized Coriolis and friction terms on the deformed hexagonal C‐grid. Q J RMeteorol Soc. 2018; 144: 2038– 2053. https://doi.org/10.1002/qj.3294



