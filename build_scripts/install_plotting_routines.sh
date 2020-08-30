aim_dir=~/compiled/game

### END OF USUAL INPUT SECTION

echo installing plotting routines ...
if [ -d $aim_dir/plotting ]
then
rm -r $aim_dir/plotting
fi
cp -r plotting $aim_dir/plotting
echo Plotting routines installed.
