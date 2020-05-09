import numpy as np;
import math as mat;
from scipy.io import netcdf;

def trafo_scalar(RES_ID, geo_props_file, variable_name_latitude, variable_name_longitude):
    NUMBER_OF_PENTAGONS = 12;
    NUMBER_OF_HEXAGONS = int(10*(mat.pow(2, 2*RES_ID) - 1));
    NUMBER_OF_SCALARS_H = NUMBER_OF_PENTAGONS + NUMBER_OF_HEXAGONS;
    lat = np.zeros([NUMBER_OF_SCALARS_H]);
    lon = np.zeros([NUMBER_OF_SCALARS_H]);
    nc_file = netcdf.netcdf_file(geo_props_file);
    latitude_scalar_pre = nc_file.variables[variable_name_latitude];
    longitude_scalar_pre = nc_file.variables[variable_name_longitude];
    latitude_scalar = latitude_scalar_pre[:].copy();
    longitude_scalar = longitude_scalar_pre[:].copy();
    del latitude_scalar_pre;
    del longitude_scalar_pre;
    nc_file.close();
    return latitude_scalar, longitude_scalar;
