LD_LIBRARY_PATH=/home/max/source/eccodes/build/lib/
export LD_LIBRARY_PATH
gcc orography_generator.c /home/max/source/eccodes/build/lib/libeccodes.so -lm -o orography_generator.exe
