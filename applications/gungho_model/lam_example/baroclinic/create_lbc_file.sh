INPUT=lfric_diag.nc
OUTPUT=lbc.nc

# Change the file from a time sequence of 3d fields(axes=space),
# to a single 4d field (axes=space+time) by changing the time dimension
# from unlimited to time=NN.
ncks --fix_rec_dmn time -O $INPUT ${INPUT}_tmp

# Change the file so that the time starts at 0, not 2000
ncap2 -O -s 'time-=2100; time_bounds-=2100' ${INPUT}_tmp ${OUTPUT}
rm ${INPUT}_tmp



