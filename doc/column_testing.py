import numpy as np;
import math;
import matplotlib.pyplot as plt;

# This is for testing vertical column solvers.

# top of atmosphere
TOA = 80e3;
# number of layers
no_of_layers = 100;
# background T
T_0 = 273.15;
# time step
delta_t = 0.1;
# number of steps you want to integrate
no_of_steps = 2000;
#name of the time integration scheme
time_integration = "Euler explicit";
# turn avection on or off
adv = 0;
# turn gravity on or off
grav_switch = 0;
# initial surface pressure
p_surf = 101325;

### END OF INPUT SECTION

# constants:

c_p = 1005;
c_v = 717.942189;
R_d = c_p - c_v;
kappa = c_p/c_v;
# gravitational acceleration
g = 9.8;

c_s = math.sqrt(kappa*R_d*T_0);

courant_number = c_s/(TOA/no_of_layers/delta_t);

print("Courant number: " + str(np.round(courant_number, 2)));

no_of_levels = no_of_layers + 1;

# grid generating procedure
z_levels = np.linspace(TOA, 0, no_of_levels);

z_layers = np.zeros([no_of_layers]);

for i in range(no_of_layers):
	z_layers[i] = 0.5*(z_levels[i] + z_levels[i + 1]);

# pracing the levels in the middle between the layers
for i in range(no_of_layers - 1):
	z_levels[i + 1] = 0.5*(z_layers[i] + z_layers[i + 1]);

# the state variables
T_old = np.zeros([no_of_layers]);
w_old = np.zeros([no_of_levels]);
rho_old = np.zeros([no_of_levels]);
T_new = np.zeros([no_of_layers]);
w_new = np.zeros([no_of_levels]);
rho_new = np.zeros([no_of_levels]);

# initializing the T field
for i in range(no_of_layers):
	T_old[i] = T_0;
	# adding a small perurbation
	if (i == 50):
		T_old[i] = T_old[i] + 0.1;
	T_new[i] = T_old[i];
# initializing the velocity field
for i in range(no_of_levels):
	w_old[i] = 0;
	w_new[i] = w_old[i];
# initializing the density field
for i in range(no_of_layers):
	rho_old[i] = p_surf/(R_d*T_old[i]);
	rho_new[i] = rho_old[i];

# loop over all time steps
e_int = np.zeros([no_of_steps]);
e_kin = np.zeros([no_of_steps]);
e_pot = np.zeros([no_of_steps]);
e_tot = np.zeros([no_of_steps]);
for i in range(no_of_steps):
	if (time_integration == "Euler explicit"):
		for j in range(no_of_layers):
			div_w = (w_old[j] - w_old[j + 1])/(z_levels[j] - z_levels[j + 1]);
			T_new[j] = T_old[j] - delta_t*R_d*T_old[j]/c_v*div_w;
		for j in np.arange(1, no_of_levels - 1):
			grad_T = (T_old[j - 1] - T_old[j])/(z_layers[j - 1] - z_layers[j]);
			w_new[j] = w_old[j] - delta_t*c_p*grad_T;
		for j in range(no_of_layers):
			T_old[j] = T_new[j];
		for j in range(no_of_levels):
			w_old[j] = w_new[j];
	for j in range(no_of_layers):	
		e_int[i] = e_int[i] + rho_old[j]*c_v*T_old[j];
		e_pot[i] = e_pot[i] + grav_switch*rho_old[j]*g*z_layers[j];
		e_kin[i] = e_kin[i] + 0.5*rho_old[j]*(0.5*w_old[j]**2 + 0.5*w_old[j + 1]**2);
	e_tot[i] = e_int[i] + e_kin[i] + e_pot[i];
	
w_rescale = 100;
z_rescale = 1e-3;
fig = plt.figure();
plt.title("Vertical velocity (last time step)");
plt.ylabel("Height / km");
plt.xlabel("Vertical velocity / cm / s");
plt.plot(w_rescale*w_new, z_rescale*z_levels);
plt.xlim([1.1*np.min(w_rescale*w_new), 1.1*np.max(w_rescale*w_new)]);
plt.ylim([np.min(z_rescale*z_levels), np.max(z_rescale*z_levels)]);
plt.show();

fig = plt.figure();
plt.title("Energy evolution");
plt.xlabel("time / s");
plt.ylabel("percentage change rel. to total init");
plt.plot(delta_t*np.arange(0, no_of_steps), (e_int- e_int[0])/e_tot[0]);
plt.plot(delta_t*np.arange(0, no_of_steps), (e_pot- e_pot[0])/e_tot[0]);
plt.plot(delta_t*np.arange(0, no_of_steps), (e_kin- e_kin[0])/e_tot[0]);
plt.plot(delta_t*np.arange(0, no_of_steps), (e_tot - e_tot[0])/e_tot[0]);
plt.xlim([0, np.max(delta_t*np.arange(0, no_of_steps))]);
plt.legend(["internal", "potential", "kinetic", "total"]);
plt.show();











