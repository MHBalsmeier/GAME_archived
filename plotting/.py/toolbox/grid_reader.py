import numpy as np;
import math as mat;
from scipy.io import netcdf;

def read_grid_props(geo_props_file, variable_name_latitude, variable_name_longitude):
    nc_file = netcdf.netcdf_file(geo_props_file);
    latitude_scalar_pre = nc_file.variables[variable_name_latitude];
    longitude_scalar_pre = nc_file.variables[variable_name_longitude];
    latitude_scalar = latitude_scalar_pre[:].copy();
    longitude_scalar = longitude_scalar_pre[:].copy();
    del latitude_scalar_pre;
    del longitude_scalar_pre;
    nc_file.close();
    return latitude_scalar, longitude_scalar;
