gcc grid_generator.c -lnetcdf -lm /lib/geos/libgeos.so /lib/conv/libconv.so -Wl,-rpath=/lib -o grid_generator
