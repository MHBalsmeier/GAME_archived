# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

cmake_minimum_required(VERSION 3.10)
set(version_major 0)
set(version_minor 1)
set(version_patch 7)
set(version ${version_major}.${version_minor}.${version_patch})
project(game VERSION ${version})
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpifort)
enable_language(Fortran)
                 # insert the root directory of RTE+RRTMGP here
set(RTERRTMGPDIR /home/max/code/rte-rrtmgp)
add_executable(
game
src/coordinator.c
src/settings.c
src/integrators/rkhevi/manage_rkhevi.c
src/integrators/rkhevi/manage_pressure_gradient.c
src/integrators/rkhevi/forward_tendencies.c
src/integrators/rkhevi/integrate_momentum.c
src/integrators/rkhevi/backward_tendencies.c
src/integrators/rkhevi/integrate_generalized_densities.c
src/integrators/rkhevi/vertical_slice_solvers.c
src/spatial_operators/vorticity_flux.c
src/spatial_operators/vorticities.c
src/spatial_operators/momentum_diff_diss.c
src/spatial_operators/scalar_times_scalar.c
src/spatial_operators/divergences.c
src/spatial_operators/scalar_times_vector.c
src/spatial_operators/gradient_operators.c
src/spatial_operators/inner_product.c
src/diagnostics/diffusion_diagnostics.c
src/diagnostics/remap_horpri2hordual_vector.c
src/diagnostics/remap_verpri2horpri_vector.c
src/diagnostics/ricci.c
src/diagnostics/tangential_wind.c
src/diagnostics/thermodynamic_diagnostics.c
src/diagnostics/vorticity_flux_horizontal.c
src/diagnostics/vorticity_flux_vertical.c
src/spatial_operators/linear_combine_two_states.c
src/io/init_state_setter.c
src/io/write_output.c
src/io/set_grid_props_and_dt.c
src/io/spatial_ops_for_output.c
src/physics/mo_radiation.f90
${RTERRTMGPDIR}/rrtmgp/mo_gas_optics.F90
${RTERRTMGPDIR}/rrtmgp/mo_gas_optics_rrtmgp.F90
${RTERRTMGPDIR}/rrtmgp/mo_rrtmgp_constants.F90
${RTERRTMGPDIR}/rrtmgp/mo_rrtmgp_util_string.F90
${RTERRTMGPDIR}/rrtmgp/mo_rrtmgp_util_reorder.F90
${RTERRTMGPDIR}/rrtmgp/mo_gas_concentrations.F90
${RTERRTMGPDIR}/rrtmgp/kernels/mo_gas_optics_kernels.F90
${RTERRTMGPDIR}/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
${RTERRTMGPDIR}/rte/mo_rte_lw.F90
${RTERRTMGPDIR}/rte/mo_rte_sw.F90
${RTERRTMGPDIR}/rte/mo_fluxes.F90
${RTERRTMGPDIR}/rte/mo_optical_props.F90
${RTERRTMGPDIR}/rte/mo_rte_config.F90
${RTERRTMGPDIR}/rte/mo_rte_kind.F90
${RTERRTMGPDIR}/rte/mo_source_functions.F90
${RTERRTMGPDIR}/rte/mo_rte_util_array.F90
${RTERRTMGPDIR}/rte/kernels/mo_fluxes_broadband_kernels.F90
${RTERRTMGPDIR}/rte/kernels/mo_optical_props_kernels.F90
${RTERRTMGPDIR}/rte/kernels/mo_rte_solver_kernels.F90
${RTERRTMGPDIR}/extensions/mo_fluxes_byband.F90
${RTERRTMGPDIR}/extensions/mo_fluxes_byband_kernels.F90
${RTERRTMGPDIR}/examples/mo_load_coefficients.F90
${RTERRTMGPDIR}/examples/mo_simple_netcdf.F90
${RTERRTMGPDIR}/extensions/mo_rrtmgp_clr_all_sky.F90
)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -Wall")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS} -lgfortran -O3")
                                                # insert output of 'nf-config --fflags --flibs' here
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm")
target_link_libraries(game eccodes m netcdf netcdff geos95 atmostracers)
install(TARGETS game DESTINATION core)
install(FILES run.sh DESTINATION core)
enable_testing()




