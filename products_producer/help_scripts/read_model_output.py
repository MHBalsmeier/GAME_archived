import help_scripts.game_grid_generator as gametrafo;
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
    
def fetch_model_output(input_file, step, short_name, level, res_id, grid_props_file):
    file = open(input_file, "rb");
    gid = ec.codes_grib_new_from_file(file);
    file.close();
    values = gid_read_values(input_file, short_name, step, level);
    variable_name_suffix = "scalar";
    if short_name == "u" or short_name == "v" or short_name == "10u" or short_name == "10v":
        variable_name_suffix = "vector";
    lat, lon = gametrafo.trafo_scalar(res_id, grid_props_file, "latitude_" + variable_name_suffix, "longitude_" + variable_name_suffix);
    return lat, lon, values;
    
    
    
