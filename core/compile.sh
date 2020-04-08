gcc src/coordinator.c src/io/set_grid_props_and_dt.c src/io/init_state_setter.c src/io/write_output.c src/time_stepping/runge_kutta_third_order.c src/time_stepping/runge_kutta_fourth_order.c src/io/interpolation_t.c src/r_operators/curl.c src/r_operators/inner.c src/r_operators/add_vector_fields.c src/r_operators/scalar_times_vector.c src/r_operators/determine_tendency.c src/r_operators/cross_product.c src/r_operators/div.c src/r_operators/laplace.c src/r_operators/grad.c src/r_operators/laplace_vec.c src/r_operators/curl_dual.c src/diagnostics/temperature_diagnostics.c src/diagnostics/global_scalar_integral.c src/diagnostics/exner_pressure_diagnostics.c src/r_operators/recover_components/recov_hor_par_dual.c src/r_operators/recover_components/recov_hor_par_pri.c src/r_operators/recover_components/recov_hor_ver_dual.c src/r_operators/recover_components/recov_hor_ver_pri.c src/r_operators/recover_components/recov_ver_0_dual.c src/r_operators/recover_components/recov_ver_0_pri.c src/r_operators/recover_components/recov_ver_1_dual.c src/r_operators/recover_components/recov_ver_1_pri.c src/r_operators/linear_combine_two_states.c -leccodes -lnetcdf -lm -Wl,-rpath=/lib /lib/indextools/libindextools.so /lib/time/libtime.so /lib/geos/libgeos.so -Wall -o game
