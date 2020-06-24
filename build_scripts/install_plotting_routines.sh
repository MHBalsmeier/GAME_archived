echo installing plotting_routines ...
aim_dir=~/compiled/game/plotting_routines
if [ -d $aim_dir ]
then
rm -r $aim_dir
fi
cp -r plotting_routines $aim_dir
echo plotting_routines installed.
