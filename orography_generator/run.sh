oro_id=2
echo "Starting orography generation ..."
./orography_generator $oro_id
if [ $? -ne 0 ]
then
echo -e ${RED}Orography file creation failed.$NC
else
echo "Orography file created sucessfully."
fi
