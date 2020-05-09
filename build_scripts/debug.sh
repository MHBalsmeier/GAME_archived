if [ -d build ]
then
rm -r build
fi
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ../core
make
ctest
make install
cd ..
rm -r build
