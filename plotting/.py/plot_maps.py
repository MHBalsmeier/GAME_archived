# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

import numpy as np;
import sys;
import toolbox.read_model_output as rmo;
import toolbox.map_properties as mp;
from matplotlib.colors import BoundaryNorm;
import cartopy.feature as cfeature;
import iris;
import toolbox.dist_stuff as ds;
import toolbox.conversions as conv;
import toolbox.time_coord_stuff as tcs;
import matplotlib.pyplot as plt;
import matplotlib as mpl;
import cartopy.crs as ccrs;
import matplotlib.offsetbox as offsetbox;
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER;
import iris.coord_systems as cs;
import iris.plot as iplt;
from scipy.interpolate import griddata;

max_interval = int(sys.argv[1]);
plot_interval = int(sys.argv[2]);
level = int(sys.argv[3]);
short_name = sys.argv[4];
grid_props_file = sys.argv[5];
save_directory = sys.argv[6];
grib_dir_name = sys.argv[7];
projection = sys.argv[8];
run_id = sys.argv[9];
uniform_range = int(sys.argv[10]);
scope = sys.argv[11];
on_pressure_bool = int(sys.argv[12]);
synoptical_time_mode = int(sys.argv[13]);
init_year = int(sys.argv[14]);
init_month = int(sys.argv[15]);
init_day = int(sys.argv[16]);
init_hour = int(sys.argv[17]);

start_timestamp = tcs.find_time_coord(init_year, init_month, init_day, init_hour, 0, 0, 0);

# default values
shift = 0;
rescale = 1;
colormap = "jet";
show_level_on = 1;
contourf_plot = 1;
gravity_mean = 9.80616;

surface_bool = 0;
if short_name == "gh":
	variable_name = "Geopotential height";
	unit_string = "gpdam";
	rescale = 1/gravity_mean;
	contourf_plot = 0;
if short_name == "t":
	variable_name = "Temperature";
	unit_string = "°C";
	shift = -conv.c2k(0);
if short_name == "pt":
	variable_name = "Potential temperature";
	unit_string = "K";
if short_name == "prmsl":
	variable_name = "MSLP / hPa";
	rescale = 0.01;
	shift = -1000;
	unit_string = "hPa";
	show_level_on = 0;
	contourf_plot = 0;
	surface_bool = 1;
if short_name == "sp":
	variable_name = "Surface pressure - 1000 hPa";
	rescale = 0.01;
	shift = -1000;
	unit_string = "hPa";
	show_level_on = 0;
	surface_bool = 1;
if short_name == "cape":
	variable_name = "CAPE";
	unit_string = "J / kg";
	show_level_on = 0;
	surface_bool = 1;
if short_name == "pres":
	variable_name = "Pressure";
	unit_string = "Pa";
if short_name == "r":
	variable_name = "Relative humidity";
	unit_string = "%";
	colormap = "Blues";
if short_name == "u":
	variable_name = "Zonal wind";
	unit_string = "m/s";
if short_name == "v":
	variable_name = "Meridional wind";
	unit_string = "m/s";
if short_name == "2t":
	variable_name = "2 m temperature";
	unit_string = "°C";
	shift = -conv.c2k(0);
	show_level_on = 0;
	surface_bool = 1;
if short_name == "vo":
	variable_name = "Relative vorticity";
	unit_string = "10^-5/s";
	rescale = 1e5;
if short_name == "10u":
	variable_name = "10 m zonal wind";
	unit_string = "m/s";
	show_level_on = 0;
	surface_bool = 1;
if short_name == "10v":
	variable_name = "10 m meridional wind";
	unit_string = "m/s";
	show_level_on = 0;
	surface_bool = 1;
if short_name == "gust":
	variable_name = "10 m gusts";
	unit_string = "m/s";
	show_level_on = 0;
	surface_bool = 1;
if short_name == "rprate":
	variable_name = "Precipitation rate (rain)";
	unit_string = "mm/h";
	rescale = conv.kgm_2s_12mmh_1(1);
	colormap = "Blues";
	show_level_on = 0;
	surface_bool = 1;
if short_name == "sprate":
	variable_name = "Precipitation rate (snow)";
	unit_string = "mm/h";
	rescale = conv.kgm_2s_12mmh_1(1);
	colormap = "Greys";
	show_level_on = 0;
	surface_bool = 1;
if short_name == "tcc":
	variable_name = "Total cloud cover";
	unit_string = "%";
	colormap = "Greys";
	show_level_on = 0;
	surface_bool = 1;
if short_name == "den":
	variable_name = "Dry air density";
	unit_string = "g/m^3";
	rescale = 1000;
if short_name == "wz":
	variable_name = "Vertical velocity";
	unit_string = "m/s";
	rescale = 1;
if short_name == "d":
	variable_name = "Horizontal divergence";
	unit_string = "1/s";
	rescale = 1;

unit_string_for_iris = unit_string;
if short_name == "gh":
    unit_string_for_iris = "dam";

disp_time_in_hr = 0;
time_unit_string = "s";

if max_interval > 24*3600 or synoptical_time_mode == 1:
    disp_time_in_hr = 1;
    time_unit_string = "hr";

savename = run_id + "_" + short_name + "_" + str(level) + "_" + scope;

if on_pressure_bool == 0:
	if surface_bool == 1:
		input_file = grib_dir_name + "/" + run_id + "+0s_surface.grb2";
	else:
		input_file = grib_dir_name + "/" + run_id + "+0s.grb2";
else:
	input_file = grib_dir_name + "/" + run_id + "+0s_pressure_levels.grb2";

	
lat, lon, values_pre = rmo.fetch_model_output(input_file, 0, short_name, level, grid_props_file);

values = np.zeros([len(lat), int(max_interval/plot_interval) + 1]);
values[:, 0] = rescale*values_pre + shift;

for i in np.arange(1, int(max_interval/plot_interval) + 1):
	time_after_init = i*plot_interval;
	if on_pressure_bool == 0:
		if surface_bool == 1:
			input_file = grib_dir_name + "/" + run_id + "+" + str(time_after_init) + "s_surface.grb2";
		else:
			input_file = grib_dir_name + "/" + run_id + "+" + str(time_after_init) + "s.grb2";
	else:
		input_file = grib_dir_name + "/" + run_id + "+" + str(time_after_init) + "s_pressure_levels.grb2";
	lat, lon, values[:, i] = rmo.fetch_model_output(input_file, time_after_init, short_name, level, grid_props_file);
	values[:, i] = rescale*values[:, i] + shift;

scope_bool_vector = np.zeros([len(values[:, 0])], dtype = bool);
if projection == "Gnomonic":
	desired_lat_deg, desired_lon_deg, height_map, width_map = mp.return_central_point(scope);
	for i in range(len(scope_bool_vector)):
		if ds.calc_distance(desired_lat_deg, desired_lon_deg, np.rad2deg(lat[i]), np.rad2deg(lon[i])) < height_map or ds.calc_distance(desired_lat_deg, desired_lon_deg, np.rad2deg(lat[i]), np.rad2deg(lon[i])) < width_map:
			scope_bool_vector[i] = True;

if uniform_range == 1:
	if projection == "Gnomonic":
		total_min = np.min(values[scope_bool_vector, :]);
		total_max = np.max(values[scope_bool_vector, :]);
	else:
		total_min = np.min(values);
		total_max = np.max(values);
	values_range_for_plot = total_max - total_min;
	if (values_range_for_plot > 0.5):
		total_min = np.floor(total_min);
		total_max = np.ceil(total_max);
		values_range_for_plot = total_max - total_min;
	if (values_range_for_plot == 0):
		values_range_for_plot = 0.1;
	color_plot_dist = values_range_for_plot/500;
	bounds = np.arange(total_min, total_max + color_plot_dist, color_plot_dist);
	color_bar_dist = values_range_for_plot/5;
	cmap = plt.get_cmap(colormap);
	norm = BoundaryNorm(bounds, ncolors = cmap.N, clip = True);

points = np.zeros([len(values[:, 0]), 2]);
fig_size = 10;
for i in range(int(max_interval/plot_interval) + 1):
	if uniform_range == 0:
		if projection == "Gnomonic":
			total_min = np.min(values[scope_bool_vector, i]);
			total_max = np.max(values[scope_bool_vector, i]);
		else:
			total_min = np.min(values[:, i]);
			total_max = np.max(values[:, i]);
		if total_min == total_max:
			total_max = total_min + 0.001;
		values_range_for_plot = total_max - total_min;
		if (values_range_for_plot > 0.5):
			total_min = np.floor(total_min);
			total_max = np.ceil(total_max);
			values_range_for_plot = total_max - total_min;
		color_plot_dist = values_range_for_plot/500;
		bounds = np.arange(total_min, total_max + color_plot_dist, color_plot_dist);
		color_bar_dist = values_range_for_plot/5;
		cmap = plt.get_cmap(colormap);
		norm = BoundaryNorm(bounds, ncolors = cmap.N, clip = True);
	time_after_init = i*plot_interval;
	print("plotting " + short_name + " at level " + str(level) + " for t - t_init = " + str(time_after_init) + " s ...");
	if (projection == "Orthographic"):
		fig = plt.figure(figsize = (fig_size, fig_size));
		coord_sys = cs.GeogCS(6371229);
	if (projection == "Mollweide"):
		fig = plt.figure(figsize = (fig_size, 0.5*fig_size));
		coord_sys = cs.GeogCS(6371229);
	if (projection == "EckertIII"):
		fig = plt.figure(figsize = (fig_size, 0.5*fig_size));
		coord_sys = cs.GeogCS(6371229);
	if (projection == "Orthographic"):
		fig = plt.figure(figsize = (fig_size, fig_size));
		ax = plt.axes(projection = ccrs.Orthographic(central_latitude = 0, central_longitude = 0));
		coord_sys = cs.GeogCS(6371229);
	if (projection == "Mollweide"):
		fig = plt.figure(figsize = (fig_size, fig_size));
		ax = plt.axes(projection = ccrs.Mollweide());
		coord_sys = cs.GeogCS(6371229);
	if (projection == "EckertIII"):
		fig = plt.figure(figsize = (fig_size, 0.5*fig_size));
		ax = plt.axes(projection = ccrs.EckertIII());
		coord_sys = cs.GeogCS(6371229);
	if (projection == "Gnomonic"):
		fig = plt.figure(figsize = (fig_size, fig_size));
		proj = ccrs.Gnomonic(central_latitude = desired_lat_deg, central_longitude = desired_lon_deg, globe = None);
		ax = plt.axes(projection = proj);
		ax.set_extent([-width_map/2, width_map/2, -height_map/2, height_map/2], crs = proj);
	lat_plot_deg = np.linspace(-90, 90, 1000);
	lon_plot_deg = np.linspace(-180, 180, 1000);
	Xi, Yi = np.meshgrid(lon_plot_deg, lat_plot_deg);
	points[:, 0] = np.rad2deg(lat);
	points[:, 1] = np.rad2deg(lon);
	values_interpolated = griddata(points, values[:, i], (Yi, Xi), method = "linear");
	if (projection != "Gnomonic"):
		lat_coord = iris.coords.DimCoord(lat_plot_deg, standard_name = "latitude", units = "degrees", coord_system = coord_sys);
		lon_coord = iris.coords.DimCoord(lon_plot_deg, standard_name = "longitude", units = "degrees", coord_system = coord_sys);
	else:
		lat_coord = iris.coords.DimCoord(lat_plot_deg, standard_name = "latitude", units = "degrees");
		lon_coord = iris.coords.DimCoord(lon_plot_deg, standard_name = "longitude", units = "degrees");
	lat_coord.guess_bounds();
	lon_coord.guess_bounds();
	gl = ax.gridlines(draw_labels = True);
	gl.xformatter = LONGITUDE_FORMATTER;
	gl.yformatter = LATITUDE_FORMATTER;
	if scope == "World":
		gl.left_labels = False;
	new_cube = iris.cube.Cube(values_interpolated, units = unit_string_for_iris, dim_coords_and_dims = [(lat_coord, 0), (lon_coord, 1)]);
	if contourf_plot == 1:
		mesh = iplt.pcolormesh(new_cube, cmap = cmap, norm = norm);
		cbar = plt.colorbar(mesh, fraction = 0.02, pad = 0.04, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist));
		cbar.ax.tick_params(labelsize = 12)
		cbar.set_label(unit_string, fontsize = 16);
	else:
		if short_name == "gh":
			levels_vector = np.arange(100, 2055, 8);
		else:
			levels_vector = np.arange(100, 2055, 4);
		big_index = np.where(levels_vector == 1012);
		basic_width = 1;
		linewidths_vector = basic_width*np.ones(np.size(levels_vector));
		linewidths_vector[big_index] = 1.5*basic_width;
		mesh = iplt.contour(new_cube, levels_vector, linewidths = linewidths_vector, colors = "black");
		plt.clabel(mesh, inline = True, fmt = "%1.0f", fontsize = 12, colors = "k");
	if (scope != "World"):
		ax.add_feature(cfeature.LAND);
		ax.add_feature(cfeature.OCEAN);
		countries = cfeature.NaturalEarthFeature(category = "cultural", name = "admin_0_countries", scale = "10m", facecolor = "none");
		ax.add_feature(countries, edgecolor = "black");
	if (scope == "CONUS" or scope == "INDIA" or scope == "CARIB" or scope == "CHINA"):
		states_provinces = cfeature.NaturalEarthFeature(category = "cultural", name = "admin_1_states_provinces_lines", scale = "10m", facecolor = "none");
		ax.add_feature(states_provinces, edgecolor = "gray");
	time_after_init_title = time_after_init;
	if disp_time_in_hr == 1:
		time_after_init_title = int(time_after_init/3600);
	if synoptical_time_mode == 0:
		time_string = "init + " + str(time_after_init_title) + " " + time_unit_string;
	if synoptical_time_mode == 1:
		time_string = "init: " + str(init_year) + "-" + str(init_month) + "-" + str(init_day) + ", " + str(init_hour) + " UTC\n";
		valid_year, valid_month, valid_day, valid_hour, dump, dump, dump = tcs.return_date(start_timestamp + time_after_init);
		time_string = time_string + "valid: " + str(valid_year) + "-" + str(valid_month) + "-" + str(valid_day) + ", " + str(valid_hour) + " UTC (+ " + str(time_after_init_title) + " hrs)";
	if show_level_on == 1:
		if on_pressure_bool == 1:
			if contourf_plot == 0:
				textstr = str(level) + " hPa " + variable_name + " / " + unit_string + " (GAME) \n" + time_string;
			if contourf_plot == 1:
				textstr = str(level) + " hPa " + variable_name + " (GAME) \n" + time_string;
		else:
			textstr = variable_name + " at level " + str(level) + "\n" + time_string;
	else:
		textstr = variable_name + " (GAME) \n" + time_string;
	ob = offsetbox.AnchoredText(textstr, loc = 3);
	ax.add_artist(ob);
	fig.savefig(save_directory + "/" + savename + "+" + str(time_after_init) + "s.png", dpi = 500, bbox_inches = "tight");
	plt.close();
	print("done");



