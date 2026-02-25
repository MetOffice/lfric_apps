
#!/bin/bash

# Code Owner: Please refer to the UM file CodeOwners.txt
# This file belongs in section: convection_comorph

# 1st argument $1 is path to compile directory

# Find full path to comorph directory containing this script
comorph=$(readlink -f $(dirname $(readlink -f "$0"))"/..")

echo "Compiling unit test at: "$comorph"/unit_tests/test_moist_proc.F90"

echo "Compile directory: "$1

cd $1

ifort -O0 -g -debug all -check all -warn interface -traceback \
      -o test_moist_proc.exe \
      \
      $comorph"/control/comorph_constants_mod.F90" \
      $comorph"/interface/standalone/qsat_data.F90" \
      $comorph"/interface/standalone/set_qsat.F90" \
      $comorph"/interface/standalone/raise_error.F90" \
      $comorph"/control/set_dependent_constants.F90" \
      $comorph"/control/diag_type_mod.F90" \
      $comorph"/control/cmpr_type_mod.F90" \
      $comorph"/util/check_bad_values.F90" \
      $comorph"/util/brent_dekker_mod.F90" \
      \
      $comorph"/moist_thermo/lat_heat_mod.F90" \
      $comorph"/moist_thermo/set_dqsatdt.F90" \
      $comorph"/moist_thermo/calc_q_tot.F90" \
      $comorph"/moist_thermo/set_cp_tot.F90" \
      $comorph"/moist_thermo/dry_adiabat.F90" \
      $comorph"/moist_thermo/calc_virt_temp_dry.F90" \
      $comorph"/moist_thermo/calc_rho_dry.F90" \
      $comorph"/moist_thermo/sat_adjust.F90" \
      \
      $comorph"/moist_proc/moist_proc_diags_type_mod.F90" \
      $comorph"/moist_proc/linear_qs_mod.F90" \
      $comorph"/moist_proc/phase_change_coefs_mod.F90" \
      $comorph"/moist_proc/calc_phase_change_coefs.F90" \
      $comorph"/moist_proc/proc_incr.F90" \
      $comorph"/moist_proc/solve_tq.F90" \
      $comorph"/moist_proc/calc_cond_temp.F90" \
      $comorph"/moist_proc/toggle_melt.F90" \
      $comorph"/moist_proc/melt_ctl.F90" \
      $comorph"/moist_proc/modify_coefs_liq.F90" \
      $comorph"/moist_proc/modify_coefs_ice.F90" \
      $comorph"/moist_proc/check_negatives.F90" \
      $comorph"/moist_proc/moist_proc_consistency.F90" \
      $comorph"/moist_proc/phase_change_solve.F90" \
      $comorph"/moist_proc/fall_in.F90" \
      $comorph"/moist_proc/fall_out.F90" \
      $comorph"/moist_proc/moist_proc_conservation.F90" \
      \
      $comorph"/microphysics/activate_cond.F90" \
      $comorph"/microphysics/ice_nucleation.F90" \
      $comorph"/microphysics/set_cond_radius.F90" \
      $comorph"/microphysics/fall_speed.F90" \
      $comorph"/microphysics/solve_wf_cond.F90" \
      $comorph"/microphysics/calc_kqkt.F90" \
      $comorph"/microphysics/calc_cond_properties.F90" \
      $comorph"/microphysics/collision_rate.F90" \
      $comorph"/microphysics/ice_rain_to_graupel.F90" \
      $comorph"/microphysics/collision_ctl.F90" \
      $comorph"/microphysics/microphysics_1.F90" \
      $comorph"/microphysics/microphysics_2.F90" \
      \
      $comorph"/moist_proc/moist_proc.F90" \
      \
      $comorph"/unit_tests/test_moist_proc.F90" \
1> std_out.txt  2> std_err.txt

echo "Done.  Compile output is in:"
echo $1/"std_out.txt"
echo $1/"std_err.txt"
