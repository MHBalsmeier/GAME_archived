# test_ids:
# 0:	standard atmosphere without orography
# 1:	standard atmosphere with Gaussian mountain
# 2:	JW test, dry, balanced
# 3:	JW test, dry, perturbed
# 4:	JW test, moist, perturbed
test_id=3
echo "Starting the test state generation ..."
./test_generator $test_id
if [ $? -ne 0 ]
then
echo -e ${RED}Test state file creation failed.$NC
else
echo "Test state file created sucessfully."
fi
