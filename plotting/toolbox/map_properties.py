def return_central_point(scope):
	if scope == "EU":
		central_lat_deg = 50;
		central_lon_deg = 8;
		height_map = 4000e3;
		width_map = 4000e3;
	if scope == "US":
		central_lat_deg = 38;
		central_lon_deg = -98;
		height_map = 5000e3;
		width_map = 5000e3;
	return central_lat_deg, central_lon_deg, height_map, width_map;
