LD_LIBRARY_PATH=/home/max/source/eccodes/build/lib/
export LD_LIBRARY_PATH
gcc test_generator.c /home/max/source/eccodes/build/lib/libeccodes.so -lm -o test_generator.exe
