# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

# This file is for plotting zonal means.

import matplotlib.pyplot as plt;

save_directory = "/home/max/held_suarez";
short_name = "u";

fig = plt.figure();
plt.title("Zonal mean of " + short_name);
plt.xlabel("latitude / deg");
plt.ylabel("height above MSL / m");
plt.savefig(save_directory + "/" + short_name);
