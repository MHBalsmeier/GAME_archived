cmake_minimum_required(VERSION 3.10)
project(core VERSION 0.0)
set(CMAKE_C_COMPILER mpicc)
add_executable(
game src/coordinator.c
src/manage_time_stepping/manage_time_stepping.c
src/manage_time_stepping/manage_interpolation.c
src/integrators/explicit_integrators/forward_tendencies.c
src/integrators/explicit_integrators/backward_tendencies.c
src/integrators/explicit_integrators/entropy_density_gas.c
src/integrators/explicit_integrators/integrate
src/integrators/explicit_integrators/integrate_continuity_dry.c
src/integrators/explicit_integrators/integrate_temp_gas.c
src/integrators/explicit_integrators/integrate_tracers.c
src/integrators/implicit_integrators/vertical_slice_solver.c
src/spatial_operators/coriolis_gen.c src/spatial_operators/potential_vorticity.c
src/spatial_operators/dissipation.c
src/spatial_operators/scalar_times_scalar.c
src/spatial_operators/divergence_operators.c
src/spatial_operators/scalar_times_vector.c
src/spatial_operators/gradient_operators.c
src/spatial_operators/
src/spatial_operators/kinetic_energy.c
src/diagnostics/thermodynamic_diagnostics/thermodynamic_diagnostics.c
src/diagnostics/dynamic_diagnostics/recovery/recov_hor_par_pri.c
src/diagnostics/dynamic_diagnostics/recovery/recov_hor_ver_pri.c
src/diagnostics/dynamic_diagnostics/recovery/recov_ver_0_curl.c
src/diagnostics/dynamic_diagnostics/recovery/recov_ver_0_pri.c
src/diagnostics/dynamic_diagnostics/recovery/recov_ver_1_curl.c
src/diagnostics/dynamic_diagnostics/recovery/recov_ver_1_pri.c
src/diagnostics/dynamic_diagnostics/recovery/trsk_modified.c
src/diagnostics/dynamic_diagnostics/ricci.c
src/diagnostics/dynamic_diagnostics/lower_boundary_solver.c
src/diagnostics/turbulence_diagnostics/diffusion_diagnostics.c
src/diagnostics/output_diagnostics/global_scalar_integral.c
src/diagnostics/output_diagnostics/linear_combine_two_states.c
src/io/init_state_setter.c
src/io/write_output.c
src/io/interpolation_t.c
src/io/set_grid_props_and_dt.c
)
add_definitions(-Wall)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
target_link_libraries(game eccodes m netcdf geos95 rte-rrtmgp-c atmostracers)
install(TARGETS game DESTINATION core)
install(FILES run.sh DESTINATION core)
install(DIRECTORY DESTINATION output)
install(DIRECTORY ../test_generator/nc_files/ DESTINATION input)
install(DIRECTORY ../grid_generator/nc_files/ DESTINATION grids)
enable_testing()
