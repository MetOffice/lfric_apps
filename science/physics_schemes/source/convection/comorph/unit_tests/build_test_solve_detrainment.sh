
#!/bin/bash

# Code Owner: Please refer to the UM file CodeOwners.txt
# This file belongs in section: convection_comorph

# 1st argument $1 is path to compile directory

# Find full path to comorph directory containing this script
comorph=$(readlink -f $(dirname $(readlink -f "$0"))"/..")

echo "Compiling unit test at: "$comorph"/unit_tests/test_solve_detrainment.F90"

echo "Compile directory: "$1

cd $1

# ifort -O0 -g -debug all -check all -warn interface -traceback \
gfortran -O0 -g -Wall -ffpe-trap=invalid,zero -fbounds-check -Warray-bounds \
         -fcheck-array-temporaries -finit-real=nan -fimplicit-none \
         -std=f2008ts -Wtabs -fbacktrace \
      -o test_solve_detrainment.exe \
      \
      $comorph"/control/comorph_constants_mod.F90" \
      \
      $comorph"/interface/standalone/raise_error.F90" \
      $comorph"/util/brent_dekker_mod.F90" \
      $comorph"/plume_model/solve_detrainment.F90" \
      \
      $comorph"/unit_tests/test_solve_detrainment.F90" \
1> std_out.txt  2> std_err.txt

echo "Done.  Compile output is in:"
echo $1/"std_out.txt"
echo $1/"std_err.txt"
