aim_dir=~/compiled/game

### END OF USUAL INPUT SECTION

if [ -d $aim_dir/core ]
then
rm -r $aim_dir/core
fi
if [ -d build ]
then
rm -r build
fi
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$aim_dir ../core
make
ctest
make install
chmod +x $aim_dir/core/run.sh
cd ..
