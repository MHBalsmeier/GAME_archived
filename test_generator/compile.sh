echo "Starting to compile test generator."
gcc test_generator.c -leccodes -lnetcdf -lm -lgeos95 -latmostracers -Wl,-rpath=/lib -Wall -o test_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Test generator compilation failed.$NC
else
echo "Test generator compiled sucessfully."
fi
