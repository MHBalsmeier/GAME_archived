oro_id=2
echo "Starting grid file creation."
mpirun -np 1 ./grid_generator $oro_id
if [ $? -ne 0 ]
then
echo -e ${RED}Grid file creation failed.$NC
else
echo "Grid file created sucessfully."
fi
