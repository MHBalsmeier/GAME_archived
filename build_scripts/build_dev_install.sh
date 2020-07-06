aim_dir=~/compiled/game_dev
rm -r $aim_dir/core
rm -r $aim_dir/input
rm -r $aim_dir/grids
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
