/*
This source file is part of the General Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/



// This is for Klemp (2008).
// This is where the damping starts in relation to the TOA. 0.75, for example, means that the upper 25 % of the atmosphere are affected by Klemp damping.
const double DAMPING_START_HEIGHT_OVER_TOA = 0.75;
// The maximum damping coefficient (the damping coefficient increases towards the TOA).
const double DAMPING_COEFF_MAX = 0.2;
// Wether or not horizontal wind divergence shall be written out.
const int WRITE_OUT_DIVV_H = 1;

int get_damping_layer_properties(double *damping_start_height_over_toa, double *damping_coeff_max)
{
	 *damping_start_height_over_toa = DAMPING_START_HEIGHT_OVER_TOA;
	 *damping_coeff_max = DAMPING_COEFF_MAX;
	return 0;
}

int get_write_settings(int *write_out_divv_h)
{
	*write_out_divv_h = WRITE_OUT_DIVV_H;
	return 0;
}
