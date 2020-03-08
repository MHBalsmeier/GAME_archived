#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include <string.h>
#include "/home/max/my_code/game/core/src/enum_and_typedefs.h"
#define FILE_NAME "nc_files/oro_0_res_2.nc"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define UNIT "GREAT_RADIAL_COORDINATE"

int main(int argc, char *argv[])
{
   int ncid, scalar_dimid, var_dimid, oro_id;
   const int NUMBER_OF_PENTAGONS = 12;
   const int NUMBER_OF_HEXAGONS = (int)(10*(pow(4, RES_ID) - 1));
   const int NUMBER_OF_SCALARS_H = 12 + NUMBER_OF_HEXAGONS;
   double *oro;
   oro = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
   int dimids[1];
   int scalar_index,  retval;
   for (int i = 0; i < NUMBER_OF_SCALARS_H; i++)
      oro[i] = 0;
    if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "scalar_index", NUMBER_OF_SCALARS_H, &scalar_dimid)))
      ERR(retval);
   dimids[0] = scalar_dimid;
   if ((retval = nc_def_var(ncid, "R", NC_DOUBLE, 1, dimids, &oro_id)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, oro_id, "units", strlen(UNIT), UNIT)))
      ERR(retval);
   if ((retval = nc_enddef(ncid)))
      ERR(retval);
   if ((retval = nc_put_var_double(ncid, oro_id, &oro[0])))
      ERR(retval);
   if ((retval = nc_close(ncid)))
      ERR(retval);
   free(oro);
   return 0;
}
