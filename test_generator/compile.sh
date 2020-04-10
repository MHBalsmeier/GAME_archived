gcc test_generator.c -leccodes -lnetcdf -lm /lib/indextools/libindextools.so /lib/geos/libgeos.so /lib/addcomp/libaddcomp.so -Wl,-rpath=/lib -Wall -o test_generator
