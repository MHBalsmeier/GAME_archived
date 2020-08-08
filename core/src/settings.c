/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/



// This is for Klemp (2008).
// This is where the damping starts in relation to the TOA. 0.75, for example, means that the upper 25 % of the atmosphere are affected by Klemp damping.
const double DAMPING_START_HEIGHT_OVER_TOA = 0.75;
// The maximum damping coefficient (the damping coefficient increases towards the TOA).
const double DAMPING_COEFF_MAX = 0.2;

int get_damping_layer_properties(double *damping_start_height_over_toa, double *damping_coeff_max)
{
	 *damping_start_height_over_toa = DAMPING_START_HEIGHT_OVER_TOA;
	 *damping_coeff_max = DAMPING_COEFF_MAX;
	return 0;
}
