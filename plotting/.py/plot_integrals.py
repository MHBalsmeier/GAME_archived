# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

import numpy as np;
import sys;
import matplotlib.pyplot as plt;
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

fig_save_path = sys.argv[1];
output_dir = sys.argv[2];
write_out_dry_mass_integral = int(sys.argv[3]);
write_out_entropy_integral = int(sys.argv[4]);
write_out_energy_integral = int(sys.argv[5]);
run_id = sys.argv[6];
write_out_linearized_entropy_integral = int(sys.argv[7]);
delta_t = float(sys.argv[8]);

fig_size = 6;
if write_out_dry_mass_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size));
	plt.title("Dry mass");
	plt.xlabel("time since init / hr");
	plt.ylabel("change relative to init value / %");
	data = np.genfromtxt(output_dir + "/dry_mass");
	time_vector = 1/3600*delta_t*(data[:, 0] - data[0, 0]);
	plt.xlim([min(time_vector), max(time_vector)]);
	dry_mass_vector = np.genfromtxt(output_dir + "/dry_mass")[:, 1];
	plt.plot(time_vector, 100*(dry_mass_vector/dry_mass_vector[0] - 1));
	plt.grid();
	fig.savefig(fig_save_path + "/" + run_id + "_dry_mass_integral.png", dpi = 500);
	plt.close();
	
if write_out_entropy_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size));
	ax = plt.axes();
	ax.grid();
	plt.title("Entropy");
	plt.xlabel("time since init / hr");
	plt.ylabel("change relative to init value / %");
	data = np.genfromtxt(output_dir + "/entropy");
	time_vector = 1/3600*delta_t*(data[:, 0] - data[0, 0]);
	plt.xlim([min(time_vector), max(time_vector)]);
	entropy_vector = data[:, 1];
	plt.plot(time_vector, 100*(entropy_vector/entropy_vector[0] - 1));
	if (write_out_linearized_entropy_integral == 1):
		data = np.genfromtxt(output_dir + "/linearized_entropy");
		linearized_entropy_vector = data[:, 1];
		plt.plot(time_vector, 100*(linearized_entropy_vector/linearized_entropy_vector[0] - 1));
		plt.legend(["entropy", "linearized entropy"], fontsize = 16);
	fig.savefig(fig_save_path + "/" + run_id + "_entropy_integral.png", dpi = 500);
	plt.close();

if write_out_energy_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size));
	plt.title("Energy");
	plt.xlabel("time since init / hr");
	plt.ylabel("change relative to init value of total energy / %");
	data = np.genfromtxt(output_dir + "/energy");
	time_vector = 1/3600*delta_t*(data[:, 0] - data[0, 0]);
	plt.xlim([min(time_vector), max(time_vector)]);
	kinetic_vector = data[:, 1];
	potential_vector = data[:, 2];
	internal_vector = data[:, 3];
	total_begin = kinetic_vector[0] + potential_vector[0] + internal_vector[0];
	plt.plot(time_vector, 100*(kinetic_vector - kinetic_vector[0])/total_begin);
	plt.plot(time_vector, 100*(potential_vector - potential_vector[0])/total_begin);
	plt.plot(time_vector, 100*(internal_vector - internal_vector[0])/total_begin);
	plt.plot(time_vector, 100*(kinetic_vector + potential_vector + internal_vector - total_begin)/total_begin);
	print("relative energy change: " + str(100*(kinetic_vector[-1] + potential_vector[-1] + internal_vector[-1] - total_begin)/total_begin) + " %");
	plt.legend(["kinetic", "potential", "internal", "total"]);
	plt.grid();
	fig.savefig(fig_save_path + "/" + run_id + "_energy_integrals.png", dpi = 500);
	plt.close();
	
	
	
	
	
	
	
	
	
	
	

