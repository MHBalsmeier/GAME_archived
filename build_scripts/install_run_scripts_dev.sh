echo installing run scripts ...
aim_dir=~/compiled/game_dev/run_scripts
if [ -d $aim_dir ]
then
rm -r $aim_dir
fi
cp -r run_scripts $aim_dir
echo run scripts installed.
