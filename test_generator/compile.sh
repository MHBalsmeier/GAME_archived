gcc test_generator.c -leccodes -lnetcdf -lm /lib/geos/libgeos.so /lib/indextools/libindextools.so /lib/addcomp/libaddcomp.so -Wl,-rpath=/lib -Wall -o test_generator
