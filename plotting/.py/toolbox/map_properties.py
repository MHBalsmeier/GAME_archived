# This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

def return_central_point(scope):
	if scope == "CEU":
		central_lat_deg = 50;
		central_lon_deg = 10;
		height_map = 3400e3;
		width_map = 3800e3;
	if scope == "CONUS":
		central_lat_deg = 38;
		central_lon_deg = -94;
		height_map = 4100e3;
		width_map = 5800e3;
	return central_lat_deg, central_lon_deg, height_map, width_map;
