/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../integrators.h"
#include "rte-rrtmgp-c.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../../diagnostics/diagnostics.h"

int backward_tendencies(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
	integrate_continuity_dry(state_old, state_new, interpolation, state_tendency, grid, dualgrid, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
    // Radiation is updated here.
    if (config_info -> rad_update == 1)
    {
        calc_rad_heating(radiation_tendency, NO_OF_SCALARS);
    }
    // Temperature diffusion gets updated here, but only at the last RK step and if heat conduction is switched on.
    if (no_rk_step == 2 && (config_info -> temperature_diff_h == 1 || config_info -> temperature_diff_v == 1))
    {
        calc_temp_diffusion_coeffs(state_old, config_info, diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v);
		grad(state_old -> temperature_gas, diagnostics -> temp_gradient, grid);
        scalar_times_vector_scalar_h_v(diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v, diagnostics -> temp_gradient, diagnostics -> temperature_flux_density, grid);
        divv_h(diagnostics -> temperature_flux_density, diffusion_info -> temp_diffusion_heating, grid);
    }
    if (config_info -> tracers_on == 1)
    {
		integrate_tracers(state_old, state_new, interpolation, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
	}
	integrate_entropy_density_dry(state_old, state_new, interpolation, state_tendency, grid, dualgrid, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
    return 0;
}
