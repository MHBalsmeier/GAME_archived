LD_LIBRARY_PATH=/home/max/custom_builds/eccodes/lib/
export LD_LIBRARY_PATH
gcc orography_generator.c $LD_LIBRARY_PATH"libeccodes.so" -o orography_generator.exe
