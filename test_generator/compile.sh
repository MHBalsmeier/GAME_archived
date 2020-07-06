echo "Starting to compile test generator."
gcc test_generator.c -leccodes -lnetcdf -lm -lgeos95 -latmostracers -Wl,-rpath=/lib -Wall -o test_generator
echo "Test generator compiled succesfully."
