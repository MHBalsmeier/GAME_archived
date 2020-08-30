aim_dir=~/compiled/game

### END OF USUAL INPUT SECTION

echo installing test initialization states ...
if [ -d $aim_dir/input ]
then
rm -r $aim_dir/input
fi
cp -r test_generator/test_states $aim_dir/input
echo Test initialization states installed.
