echo "Starting to compile grid generator ..."
mpicc -O3 grid_generator.c -lnetcdf -lm -lgeos95 -Wall -o grid_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Grid generator compilation failed.$NC
else
echo "Grid generator compiled successfully."
fi
