aim_dir=~/compiled/game

### END OF USUAL INPUT SECTION

echo installing run scripts ...
if [ -d $aim_dir/run_scripts ]
then
rm -r $aim_dir/run_scripts
fi
cp -r run_scripts $aim_dir/run_scripts
echo run scripts installed.
