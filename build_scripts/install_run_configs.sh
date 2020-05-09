echo installing run_configs ...
aim_dir=~/compiled/game/run_configs
if [ -d $aim_dir ]
then
rm -r $aim_dir
fi
cp -r run_configs $aim_dir
echo run_configs installed.
