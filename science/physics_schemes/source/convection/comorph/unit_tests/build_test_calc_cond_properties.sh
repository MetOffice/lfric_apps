
#!/bin/bash

# Code Owner: Please refer to the UM file CodeOwners.txt
# This file belongs in section: convection_comorph

# 1st argument $1 is path to compile directory

# Find full path to comorph directory containing this script
comorph=$(readlink -f $(dirname $(readlink -f "$0"))"/..")

echo "Compiling unit test at: "$comorph"/unit_tests/test_calc_cond_properties.F90"

echo "Compile directory: "$1

cd $1

ifort -O0 -g -debug all -check all -warn interface -traceback \
      -o test_calc_cond_properties.exe \
      \
      $comorph"/interface/standalone/raise_error.F90" \
      $comorph"/util/brent_dekker_mod.F90" \
      \
      $comorph"/control/comorph_constants_mod.F90" \
      $comorph"/control/set_dependent_constants.F90" \
      \
      $comorph"/moist_proc/fall_out.F90" \
      \
      $comorph"/microphysics/set_cond_radius.F90" \
      $comorph"/microphysics/fall_speed.F90" \
      $comorph"/microphysics/solve_wf_cond.F90" \
      $comorph"/microphysics/calc_kqkt.F90" \
      $comorph"/microphysics/calc_cond_properties.F90" \
      \
      $comorph"/unit_tests/test_calc_cond_properties.F90" \
1> std_out.txt  2> std_err.txt

echo "Done.  Compile output is in:"
echo $1/"std_out.txt"
echo $1/"std_err.txt"

