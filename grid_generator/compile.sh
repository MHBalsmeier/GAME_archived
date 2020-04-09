gcc grid_generator.c -lnetcdf -lm /lib/geos/libgeos.so /lib/conv/libconv.so /lib/indextools/libindextools.so -Wl,-rpath=/lib -Wall -o grid_generator
