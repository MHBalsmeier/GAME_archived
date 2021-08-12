# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

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
	if scope == "CHINA":
		central_lat_deg = 36;
		central_lon_deg = 104;
		height_map = 4500e3;
		width_map = 5800e3;
	if scope == "OZ":
		central_lat_deg = -27;
		central_lon_deg = 135;
		height_map = 4500e3;
		width_map = 5800e3;
	if scope == "CARIB":
		central_lat_deg = 24;
		central_lon_deg = -67;
		height_map = 4500e3;
		width_map = 5800e3;
	if scope == "WRUS":
		central_lat_deg = 59;
		central_lon_deg = 44;
		height_map = 4500e3;
		width_map = 5800e3;
	if scope == "GULF":
		central_lat_deg = 26;
		central_lon_deg = 46;
		height_map = 4500e3;
		width_map = 5800e3;
	if scope == "OCEAN":
		central_lat_deg = 6;
		central_lon_deg = 118;
		height_map = 5000e3;
		width_map = 5800e3;
	if scope == "INDIA":
		central_lat_deg = 23;
		central_lon_deg = 77;
		height_map = 4500e3;
		width_map = 5800e3;
	if scope == "CAFR":
		central_lat_deg = 9;
		central_lon_deg = 21;
		height_map = 4500e3;
		width_map = 5800e3;
	if scope == "SAFR":
		central_lat_deg = -17;
		central_lon_deg = 24;
		height_map = 4500e3;
		width_map = 5800e3;
	if scope == "SAM":
		central_lat_deg = -15;
		central_lon_deg = -60;
		height_map = 4500e3;
		width_map = 5800e3;
	return central_lat_deg, central_lon_deg, height_map, width_map;




