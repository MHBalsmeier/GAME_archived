# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

import eccodes as ec;
import numpy as np;
from colorama import Fore;
from colorama import Style;

def gid_return_extended(filename, short_name, step, level):
    filee = open(filename, "rb");
    for j in np.arange(0, ec.codes_count_in_file(filee)):
        gid = ec.codes_grib_new_from_file(filee, headers_only = False);
        if ec.codes_get(gid, "shortName") == short_name and ec.codes_get(gid, "forecastTime") == step and ec.codes_get(gid, "level") == level:
            filee.close();
            return gid;
        else:
            ec.codes_release(gid);
    filee.close();
    exit(1);

def gid_read_values(filename, short_name, step, level):
    gid = gid_return_extended(filename, short_name, step, level);
    return_array = ec.codes_get_array(gid, "values");
    ec.codes_release(gid);
    return return_array;
    
def fetch_model_output(input_file, step, short_name, level):
    file = open(input_file, "rb");
    gid = ec.codes_grib_new_from_file(file);
    file.close();
    values = gid_read_values(input_file, short_name, step, level);
	lat_vector_deg = np.deg2rad(ec.codes_get_array(gid, "latitudes"));
	lon_vector_deg = np.deg2rad(ec.codes_get_array(gid, "longitudes"));
    return lat, lon, values;
    
    
    
