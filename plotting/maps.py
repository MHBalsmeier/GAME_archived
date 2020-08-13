import numpy as np;
import sys;
import toolbox.read_model_output as rmo;
from matplotlib.colors import BoundaryNorm;
import cartopy.feature as cfeature;
import iris;
import toolbox.game_grid_generator as ggg;
import toolbox.dist_stuff as ds;
import toolbox.conversions as conv;
import matplotlib.pyplot as plt;
import matplotlib as mpl;
import cartopy.crs as ccrs;
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER;
from iris.coord_systems import GeogCS;
import iris.plot as iplt;
from scipy.interpolate import griddata;

max_interval = int(sys.argv[1]);
plot_interval = int(sys.argv[2]);
level = int(sys.argv[3]);
short_name = sys.argv[4];
grid_props_file = sys.argv[5];
save_folder = sys.argv[6];
grib_dir_name = sys.argv[7];
projection = sys.argv[8];
run_id = sys.argv[9];

# default values
shift = 0;
rescale = 1;
colormap = "jet";

if short_name == "pt":
    variable_name = "Potential temperature";
    unit_string = "K";
if short_name == "prmsl":
    variable_name = "MSLP - 1000 hPa";
    rescale = 0.01;
    shift = -1000;
    unit_string = "hPa";
if short_name == "sp":
    variable_name = "Surface pressure - 1000 hPa";
    rescale = 0.01;
    shift = -1000;
    unit_string = "hPa";
if short_name == "cape":
    variable_name = "CAPE";
    unit_string = "J / kg";
if short_name == "pres":
    variable_name = "Pressure";
    unit_string = "Pa";
if short_name == "rh":
    variable_name = "Relative humidity";
    unit_string = "%";
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
if short_name == "vo":
    variable_name = "Relative vorticity";
    unit_string = "10^-5/s";
    rescale = 1e5;
if short_name == "10u":
    variable_name = "10 m zonal wind";
    unit_string = "m/s";
if short_name == "10v":
    variable_name = "10 m meridional wind";
    unit_string = "m/s";
if short_name == "gust":
    variable_name = "10 m gusts";
    unit_string = "m/s";
if short_name == "rprate":
    variable_name = "Liquid precipitation rate";
    unit_string = "mm/h";
    rescale = conv.kgm_2s_12mmh_1(1);
    colormap = "Greys";
if short_name == "sprate":
    variable_name = "Solid precipitation rate";
    unit_string = "mm/h";
    rescale = conv.kgm_2s_12mmh_1(1);
    colormap = "Greys";
if short_name == "tcc":
    variable_name = "Total cloud cover";
    unit_string = "%";
    colormap = "Greys";
if short_name == "den":
    variable_name = "Dry air density";
    unit_string = "g/m^3";
    rescale = 1000;
if short_name == "wz":
    variable_name = "Vertical velocity";
    unit_string = "m/s";
    rescale = 1;

disp_time_in_hr = 0;
time_unit_string = "s";

if max_interval > 24*3600:
    disp_time_in_hr = 1;
    time_unit_string = "hr";

savename = run_id + "_" + short_name + "_" + str(level);

input_file = grib_dir_name + "/" + run_id + "+0s.grb2";
lat, lon, values_pre = rmo.fetch_model_output(input_file, 0, short_name, level, 4, grid_props_file);

values = np.zeros([len(lat), int(max_interval/plot_interval) + 1]);
values[:, 0] = rescale*values_pre + shift;

for i in np.arange(1, int(max_interval/plot_interval) + 1):
    time_after_init = i*plot_interval;
    input_file = grib_dir_name + "/" + run_id + "+" + str(time_after_init) + "s.grb2";
    lat, lon, values[:, i] = rmo.fetch_model_output(input_file, time_after_init, short_name, level, 4, grid_props_file);
    values[:, i] = rescale*values[:, i] + shift;

total_min = np.min(values);
total_max = np.max(values);
values_range_for_plot = total_max - total_min;
if (values_range_for_plot > 0.5):
	total_min = np.floor(np.min(values));
	total_max = np.ceil(np.max(values));
	values_range_for_plot = total_max - total_min;
color_plot_dist = values_range_for_plot/500;
bounds = np.arange(total_min, total_max + color_plot_dist, color_plot_dist);
color_bar_dist = values_range_for_plot/5;

cmap = plt.get_cmap(colormap);
norm = BoundaryNorm(bounds, ncolors = cmap.N, clip = True);
points = np.zeros([len(values[:, 0]), 2]);
fig_size = 10;
for i in range(int(max_interval/plot_interval) + 1):
	time_after_init = i*plot_interval;
	print("plotting " + short_name + " at level " + str(level) + " for t - t_init = " + str(time_after_init) + " s ...");
	if (projection == "Orthographic"):
		fig = plt.figure(figsize = (fig_size, fig_size));
	if (projection == "Mollweide"):
		fig = plt.figure(figsize = (fig_size, 0.5*fig_size));
	if (projection == "EckertIII"):
		fig = plt.figure(figsize = (fig_size, 0.5*fig_size));
	if (projection == "Orthographic"):
		ax = plt.axes(projection=ccrs.Orthographic(central_latitude = 0, central_longitude = 0));
	if (projection == "Mollweide"):
		ax = plt.axes(projection=ccrs.Mollweide());
	if (projection == "EckertIII"):
		ax = plt.axes(projection=ccrs.EckertIII());
	gl = ax.gridlines(draw_labels = True);
	gl.xformatter = LONGITUDE_FORMATTER;
	gl.yformatter = LATITUDE_FORMATTER;
	lat_plot_deg = np.linspace(-90, 90, 1000);
	lon_plot_deg = np.linspace(-180, 180, 1000);
	Xi, Yi = np.meshgrid(lon_plot_deg, lat_plot_deg);
	points[:, 0] = np.rad2deg(lat);
	points[:, 1] = np.rad2deg(lon);
	values_interpolated = griddata(points, values[:, i], (Yi, Xi), method = "linear");
	cs = GeogCS(6371229);
	lat_coord = iris.coords.DimCoord(lat_plot_deg, standard_name = "latitude", units = "degrees", coord_system = cs);
	lon_coord = iris.coords.DimCoord(lon_plot_deg, standard_name = "longitude", units = "degrees", coord_system = cs);
	lat_coord.guess_bounds();
	lon_coord.guess_bounds();
	new_cube = iris.cube.Cube(values_interpolated, units = unit_string, dim_coords_and_dims = [(lat_coord, 0), (lon_coord, 1)]);
	mesh = iplt.pcolormesh(new_cube, cmap = cmap, norm = norm);
	cbar = plt.colorbar(mesh, fraction = 0.02, pad = 0.04, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist));
	cbar.ax.tick_params(labelsize = 12)
	cbar.set_label(unit_string, fontsize = 16);
	time_after_init_title = time_after_init;
	if disp_time_in_hr == 1:
		time_after_init_title = int(time_after_init/3600);
	plt.title(variable_name + "; " + str(time_after_init_title) + " " + time_unit_string + " after init", fontsize = 20);
	fig.savefig(save_folder + "/" + savename + "+" + str(time_after_init) + "s.png", dpi = 500);
	plt.close();
	print("done");



