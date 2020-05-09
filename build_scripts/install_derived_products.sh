echo installing derived_products ...
aim_dir=~/compiled/game/derived_products
if [ -d $aim_dir ]
then
rm -r $aim_dir
fi
cp -r derived_products $aim_dir
echo derived_products installed.
