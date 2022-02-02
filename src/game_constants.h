/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file collects physical constants that are hardly ever modified.
*/

// calculating the Earth radius
#define SEMIMAJOR 6378137.0
#define SEMIMINOR 6356752.314
#define RADIUS pow(SEMIMAJOR*SEMIMAJOR*SEMIMINOR, 1.0/3)

// other constants
#define OMEGA (7.292115e-5)
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define P_0 100000.0
#define GRAVITY_MEAN_SFC_ABS 9.80616
#define EPSILON_SECURITY (1e-10)
#define SECONDS_PER_HOUR 3600
#define RHO_WATER 1024.0
#define ROUGHNESS_LENGTH_GRASS 0.02
