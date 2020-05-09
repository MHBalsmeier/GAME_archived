echo installing output ...
aim_dir=~/compiled/game/output
if [ -d $aim_dir ]
then
rm -r $aim_dir
fi
mkdir $aim_dir
echo output installed.
