#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "/home/max/my_code/game/core/src/enum_and_typedefs.h"
#include "/home/max/custom_builds/eccodes/include/eccodes.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

int main(int argc, char *argv[])
{
   const int START_SECTION = 4;
   double *oro;
   oro = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
   int dimids[1];
   int scalar_index,  retval;
   for (int i = 0; i < NUMBER_OF_SCALARS_H; i++)
      oro[i] = 0;
   int error = 0;
   FILE *oro_file;
   if(error != 0)
       ECCERR(error);
   char *SAMPLE_FILE = "grib_files/scalar_field_blueprint_res_id_2.grb2";
   FILE *SAMPLE;
   SAMPLE = fopen(SAMPLE_FILE,"r");
   int err = 0;
   FILE *OUT_GRIB;
   char *OUTPUT_FILE = "grib_files/orography_0_res_id_2.grb2";
   OUT_GRIB = fopen(OUTPUT_FILE, "w");
   grib_multi_handle *multi_handle = NULL;
   multi_handle = codes_grib_multi_handle_new(NULL);
   codes_handle *handle_orography = NULL;
   handle_orography = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
   if (retval = codes_set_long(handle_orography, "discipline", 2))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "subCentre", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "centre", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "significanceOfReferenceTime", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "dataDate", 9999))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "dataTime", 9999))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "productionStatusOfProcessedData", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "typeOfProcessedData", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "grib2LocalSectionNumber", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "shapeOfTheEarth", 7))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "productDefinitionTemplateNumber", 0))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "typeOfGeneratingProcess", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "indicatorOfUnitOfTimeRange", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "stepUnits", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "forecastTime", 0))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "stepRange", 0))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "generatingProcessIdentifier", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "numberOfGridUsed", 20))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "numberOfGridInReference", 1))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "parameterCategory", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "typeOfFirstFixedSurface", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "scaleFactorOfFirstFixedSurface", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "scaledValueOfFirstFixedSurface", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "typeOfSecondFixedSurface", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "scaleFactorOfSecondFixedSurface", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "scaledValueOfSecondFixedSurface", 255))
       ECCERR(retval);
   if (retval = codes_set_long(handle_orography, "level", 0))
       ECCERR(retval);
   if (retval = codes_set_double_array(handle_orography, "values", oro, NUMBER_OF_SCALARS_H))
       ECCERR(retval);
   codes_grib_multi_handle_append(handle_orography, START_SECTION, multi_handle);
   codes_handle_delete(handle_orography);
   codes_grib_multi_handle_write(multi_handle, OUT_GRIB);
   codes_grib_multi_handle_delete(multi_handle);
   fclose(SAMPLE);
   fclose(OUT_GRIB);
   free(oro);
   return 0;
}
