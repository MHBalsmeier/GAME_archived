/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "grid_generator.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "geos95.h"
#include <netcdf.h>
#include <math.h>
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int set_f_vec(double latitude_vector[], double direction_dual[], double latitude_vector_dual[], double f_vec[])
{
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H; ++i)
    {
        if (i >= NUMBER_OF_DUAL_VECTORS_H)
        	f_vec[i] = 2*OMEGA*sin(latitude_vector[i - NUMBER_OF_DUAL_VECTORS_H]);
        else
        	f_vec[i] = 2*OMEGA*cos(latitude_vector_dual[i])*sin(direction_dual[i]);
    }
 	return 0;   
}


int set_orography(int RES_ID, int ORO_ID, double z_surface[])
{
    if (ORO_ID != 0)
    {
    	int retval, z_surface_id, ncid;
		int ORO_FILE_LENGTH = 100;
		char *ORO_FILE_PRE = malloc((ORO_FILE_LENGTH + 1)*sizeof(char));
		sprintf(ORO_FILE_PRE, "../orography_generator/nc_files/B%d_O%d.nc", RES_ID, ORO_ID);
		ORO_FILE_LENGTH = strlen(ORO_FILE_PRE);
		free(ORO_FILE_PRE);
		char *ORO_FILE = malloc((ORO_FILE_LENGTH + 1)*sizeof(char));
		sprintf(ORO_FILE, "../orography_generator/nc_files/B%d_O%d.nc", RES_ID, ORO_ID);	    
		if ((retval = nc_open(ORO_FILE, NC_NOWRITE, &ncid)))
		    ERR(retval);
		if ((retval = nc_inq_varid(ncid, "z_surface", &z_surface_id)))
		    ERR(retval);		
		if ((retval = nc_get_var_double(ncid, z_surface_id, &z_surface[0])))
		    ERR(retval);
		if ((retval = nc_close(ncid)))
		    ERR(retval);
    }
    else
    {
    	for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    	{
    		z_surface[i] = 0;
    	}
    }
	return 0;
}








