# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

oro_id=1
echo "Starting orography generation ..."
./orography_generator $oro_id
if [ $? -ne 0 ]
then
echo -e ${RED}Orography file creation failed.$NC
else
echo "Orography file created sucessfully."
fi
