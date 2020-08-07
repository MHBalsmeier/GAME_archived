/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../integrators.h"
#include "rte-rrtmgp-c.h"
#include "../../diagnostics/diagnostics.h"
#include "../../spatial_operators/spatial_operators.h"

int backward_tendencies(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
	integrate_continuity_dry(state_old, state_new, interpolation, state_tendency, grid, dualgrid, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
    // Radiation is updated here.
    if (config_info -> rad_update == 1)
    {
        calc_rad_heating(radiation_tendency, NO_OF_SCALARS);
    }
    double rho_h, c_h_v;
    // Diffusion gets updated here.
    if (config_info -> scalar_diffusion_on == 1 && no_rk_step == 2)
    {
        calc_temp_diffusion_coeffs(state_old -> temp_gas, state_old -> density_dry, diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v);
		grad(state_old -> temp_gas, diagnostics -> temp_gradient, grid);
        scalar_times_vector_scalar_h_v(diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v, diagnostics -> temp_gradient, diagnostics -> temperature_flux_density, grid);
        divv_h(diagnostics -> temperature_flux_density, diffusion_info -> temp_diffusion_heating, grid);
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (config_info -> tracers_on == 1)
			{
				rho_h = state_old -> density_dry[i] + state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
				c_h_v = spec_heat_cap_diagnostics_v(state_old -> density_dry[i], state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
				diffusion_info -> temp_diffusion_heating[i] = rho_h*c_h_v*diffusion_info -> temp_diffusion_heating[i];
			}
			else
				diffusion_info -> temp_diffusion_heating[i] = state_old -> density_dry[i]*C_D_V*diffusion_info -> temp_diffusion_heating[i];
		}
    }
    if (config_info -> tracers_on == 1)
    {
		integrate_tracers(state_old, state_new, interpolation, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
	}
	integrate_entropy_density_dry(state_old, state_new, interpolation, state_tendency, grid, dualgrid, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
    return 0;
}
