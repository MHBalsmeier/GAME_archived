/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../integrators.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../../diagnostics/diagnostics.h"
#include <stdio.h>
#include <stdlib.h>

int backward_tendencies(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int no_rk_step, double time_coordinate)
{
    // Radiation is updated here.
    if (config_info -> rad_on == 1 && config_info -> rad_update == 1 && no_rk_step == 0)
    {
    	printf("Starting update of radiative fluxes ...\n");
    	// Fortran needs pointers, this is why this is necessary
    	int no_of_scalars = NO_OF_SCALARS;
    	int no_of_vectors = NO_OF_VECTORS;
    	int no_of_vectors_per_layer = NO_OF_VECTORS_PER_LAYER;
    	int no_of_constituents = NO_OF_CONSTITUENTS;
    	int no_of_condensed_constituents = NO_OF_CONDENSED_CONSTITUENTS;
    	int no_of_layers = NO_OF_LAYERS;
		calc_radiative_flux_convergence(grid -> latitude_scalar, grid -> longitude_scalar, grid -> z_scalar, grid -> z_vector,
		state -> mass_densities, state -> temperature_gas, radiation_tendency, &no_of_scalars, &no_of_vectors, &no_of_vectors_per_layer, &no_of_layers, &no_of_constituents, 
		&no_of_condensed_constituents, &time_coordinate);
    	printf("Update of radiative fluxes completed.\n");
    }
    // Temperature diffusion gets updated here, but only at the first RK step and if heat conduction is switched on.
    if (no_rk_step == 0 && (config_info -> temperature_diff_h == 1 || config_info -> temperature_diff_v == 1))
    {
    	// Now we need to calculate the scalar diffusion coefficients.
        calc_temp_diffusion_coeffs(state, config_info, irreversible_quantities -> scalar_diffusion_coeff_numerical_h, irreversible_quantities -> scalar_diffusion_coeff_numerical_v);
        // The diffusion of the temperature depends on its gradient.
		grad(state -> temperature_gas, diagnostics -> temperature_gradient, grid);
		// Now the diffusive temperature flux density can be obtained.
        scalar_times_vector_scalar_h_v(irreversible_quantities -> scalar_diffusion_coeff_numerical_h, irreversible_quantities -> scalar_diffusion_coeff_numerical_v, diagnostics -> temperature_gradient, diagnostics -> flux_density, grid);
        // The divergence of the diffusive temperature flux density is the diffusive temperature heating.
        divv_h(diagnostics -> flux_density, irreversible_quantities -> temperature_diffusion_heating, grid);
        add_vertical_divv(diagnostics -> flux_density, irreversible_quantities -> temperature_diffusion_heating, grid);
    }
    
	integrate_generalized_densities(state, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, irreversible_quantities, config_info, no_rk_step);
    return 0;
}
