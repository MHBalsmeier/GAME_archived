#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "/home/max/my_code/game/core/src/enum_and_typedefs.h"
#include "/home/max/custom_builds/eccodes/include/eccodes.h"
#define FILE_NAME "grib_files/scalar_field_blueprint_res_id_2.grb2"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}
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
   int error = 0;
   FILE *oro_file;
   oro_file = fopen(FILE_NAME, "r");
   codes_handle *handle = NULL;
   handle = codes_handle_new_from_file(NULL, oro_file, PRODUCT_GRIB, &error);
   if(error != 0)
       ECCERR(error);
   long n_lat = 100;
   codes_set_long(handle, "Nj", n_lat);
   codes_get_long(handle, "Nj", &n_lat);
   printf("%ld\n", n_lat);
   codes_set_double_array(handle, "values", &oro, NUMBER_OF_SCALARS_H);
   free(oro);
   return 0;
}
