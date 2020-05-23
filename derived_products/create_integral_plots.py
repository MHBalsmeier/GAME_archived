import numpy as np;
import sys;
import matplotlib.pyplot as plt;

fig_save_path = sys.argv[1];
output_dir = sys.argv[2];
write_out_dry_mass_integral = int(sys.argv[3]);
write_out_entropy_integral = int(sys.argv[4]);
write_out_energy_integral = int(sys.argv[5]);

fig_size = 6;
if write_out_dry_mass_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size));
	plt.title("Dry mass");
	plt.xlabel("time since init / s");
	plt.ylabel("relative to init value");
	data = np.genfromtxt(output_dir + "/dry_mass");
	time_vector = data[:, 0] - data[0, 0];
	plt.xlim([min(time_vector), max(time_vector)]);
	dry_mass_vector = np.genfromtxt(output_dir + "/dry_mass")[:, 1];
	plt.plot(time_vector, dry_mass_vector/dry_mass_vector[0]);
	fig.savefig(fig_save_path + "/dry_mass_integral.png", dpi = 500);
	plt.close();
	
if write_out_entropy_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size));
	plt.title("Entropy");
	plt.xlabel("time since init / s");
	plt.ylabel("relative to init value");
	data = np.genfromtxt(output_dir + "/entropy");
	time_vector = data[:, 0] - data[0, 0];
	plt.xlim([min(time_vector), max(time_vector)]);
	entropy_vector = data[:, 1];
	plt.plot(time_vector, entropy_vector/entropy_vector[0]);
	fig.savefig(fig_save_path + "/entropy_integral.png", dpi = 500);
	plt.close();

if write_out_energy_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size));
	plt.title("Energy");
	plt.xlabel("time since init / s");
	plt.ylabel("relative to init value");
	data = np.genfromtxt(output_dir + "/energy");
	time_vector = data[:, 0] - data[0, 0];
	plt.xlim([min(time_vector), max(time_vector)]);
	kinetic_vector = data[:, 1];
	potential_vector = data[:, 2];
	internal_vector = data[:, 3];
	total_begin = kinetic_vector[0] + potential_vector[0] + internal_vector[0];
	plt.plot(time_vector, kinetic_vector/total_begin);
	plt.plot(time_vector, potential_vector/total_begin);
	plt.plot(time_vector, internal_vector/total_begin);
	plt.plot(time_vector, (kinetic_vector + potential_vector + internal_vector)/total_begin);
	plt.legend(["kinetic", "potential", "internal", "total"]);
	fig.savefig(fig_save_path + "/energy_integrals.png", dpi = 500);
	plt.close();
	
	
	
	
	
	
	
	
	
	
	

