# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "Starting to compile orography generator ..."
gcc orography_generator.c -lnetcdf -lgeos95 -lm -o orography_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Orography generator compilation failed.$NC
else
echo "Orography generator compiled successfully."
fi
