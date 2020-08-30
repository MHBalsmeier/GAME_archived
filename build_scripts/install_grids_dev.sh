aim_dir=~/compiled/game_dev

### END OF USUAL INPUT SECTION

echo installing grids ...
if [ -d $aim_dir/grids ]
then
rm -r $aim_dir/grids
fi
cp -r grid_generator/grids $aim_dir/grids
echo Grids installed.
