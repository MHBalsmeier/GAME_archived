/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "enum_and_typedefs.h"
#include "io/io.h"
#include "spatial_operators/spatial_operators.h"
#include "diagnostics/diagnostics.h"
#include "time_stepping/time_stepping.h"
#include "rte-rrtmgp-c.h"
#include <mpi.h>
