echo installing plotting routines ...
aim_dir=~/compiled/game_dev/plotting
if [ -d $aim_dir ]
then
rm -r $aim_dir
fi
cp -r plotting $aim_dir
echo Plotting routines installed.
