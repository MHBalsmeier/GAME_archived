echo "Starting to compile grid generator ..."
mpicc -O3 grid_generator.c -lnetcdf -lm -lgeos95 -Wall -o grid_generator
echo "Grid generator compiled successfully."
