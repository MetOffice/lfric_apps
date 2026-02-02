#!/bin/bash

# Code Owner: Please refer to the UM file CodeOwners.txt
# This file belongs in section: convection_comorph

# 1st argument $1 is path to compile directory

# Find full path to comorph directory containing this script
comorph=$(readlink -f "$(dirname "$(readlink -f "$0")")""/..")

echo "Compiling unit test at: ""$comorph""/unit_tests/test_comorph.F90"

echo "Compile directory: ""$1"

cd "$1" || exit

#module swap ifort ifort/19.0_64

export OMP_NUM_THREADS=2

#ifort -O0 -g -debug all -check all -warn interface -traceback \
#      -init=snan,arrays -init=huge -check bounds -check shape \
#      -parallel -qopenmp \
gfortran -O0 -g -Wall -ffpe-trap=invalid,zero -fbounds-check -Warray-bounds \
         -fcheck-array-temporaries -finit-real=nan -fimplicit-none \
         -std=f2008ts -Wtabs -fbacktrace -fopenmp \
    -o test_comorph.exe \
      \
      "$comorph""/control/comorph_constants_mod.F90" \
      "$comorph""/interface/standalone/raise_error.F90" \
      "$comorph""/interface/standalone/qsat_data.F90" \
      "$comorph""/interface/standalone/set_qsat.F90" \
      "$comorph""/control/set_dependent_constants.F90" \
      "$comorph""/control/cmpr_type_mod.F90" \
      \
      "$comorph""/util/brent_dekker_mod.F90" \
      "$comorph""/util/compress.F90" \
      "$comorph""/util/decompress.F90" \
      "$comorph""/util/copy_field.F90" \
      "$comorph""/util/diff_field.F90" \
      "$comorph""/util/init_zero.F90" \
      "$comorph""/util/check_bad_values.F90" \
      \
      "$comorph""/moist_thermo/lat_heat_mod.F90" \
      "$comorph""/moist_thermo/set_dqsatdt.F90" \
      "$comorph""/moist_thermo/calc_q_tot.F90" \
      "$comorph""/moist_thermo/set_cp_tot.F90" \
      "$comorph""/moist_thermo/dry_adiabat.F90" \
      "$comorph""/moist_thermo/sat_adjust.F90" \
      "$comorph""/moist_thermo/calc_virt_temp_dry.F90" \
      "$comorph""/moist_thermo/calc_virt_temp.F90" \
      "$comorph""/moist_thermo/linear_qs_mod.F90" \
      "$comorph""/moist_thermo/calc_qvl_supersat.F90" \
      "$comorph""/moist_thermo/calc_rho_dry.F90" \
      "$comorph""/moist_thermo/calc_layer_mass.F90" \
      \
      "$comorph""/control/grid_type_mod.F90" \
      "$comorph""/control/fields_type_mod.F90" \
      "$comorph""/control/env_half_mod.F90" \
      "$comorph""/control/cloudfracs_type_mod.F90" \
      "$comorph""/control/turb_type_mod.F90" \
      "$comorph""/control/parcel_type_mod.F90" \
      "$comorph""/control/res_source_type_mod.F90" \
      "$comorph""/control/subregion_mod.F90" \
      "$comorph""/control/fields_2d_mod.F90" \
      "$comorph""/control/diag_type_mod.F90" \
      "$comorph""/control/fields_diags_type_mod.F90" \
      "$comorph""/control/parcel_diags_type_mod.F90" \
      "$comorph""/control/subregion_diags_type_mod.F90" \
      "$comorph""/control/diags_2d_type_mod.F90" \
      "$comorph""/control/diags_super_type_mod.F90" \
      "$comorph""/control/set_l_within_bl.F90" \
      "$comorph""/control/force_cloudfrac_consistency.F90" \
      "$comorph""/control/set_cloudfracs_k.F90" \
      \
      "$comorph""/moist_proc/moist_proc_diags_type_mod.F90" \
      "$comorph""/moist_proc/phase_change_coefs_mod.F90" \
      "$comorph""/moist_proc/calc_phase_change_coefs.F90" \
      "$comorph""/moist_proc/proc_incr.F90" \
      "$comorph""/moist_proc/solve_tq.F90" \
      "$comorph""/moist_proc/calc_cond_temp.F90" \
      "$comorph""/moist_proc/toggle_melt.F90" \
      "$comorph""/moist_proc/melt_ctl.F90" \
      "$comorph""/moist_proc/modify_coefs_liq.F90" \
      "$comorph""/moist_proc/modify_coefs_ice.F90" \
      "$comorph""/moist_proc/check_negatives.F90" \
      "$comorph""/moist_proc/moist_proc_consistency.F90" \
      "$comorph""/moist_proc/phase_change_solve.F90" \
      "$comorph""/moist_proc/fall_in.F90" \
      "$comorph""/moist_proc/fall_out.F90" \
      "$comorph""/moist_proc/moist_proc_conservation.F90" \
      \
      "$comorph""/microphysics/activate_cond.F90" \
      "$comorph""/microphysics/ice_nucleation.F90" \
      "$comorph""/microphysics/set_cond_radius.F90" \
      "$comorph""/microphysics/fall_speed.F90" \
      "$comorph""/microphysics/solve_wf_cond.F90" \
      "$comorph""/microphysics/calc_kqkt.F90" \
      "$comorph""/microphysics/calc_cond_properties.F90" \
      "$comorph""/microphysics/collision_rate.F90" \
      "$comorph""/microphysics/ice_rain_to_graupel.F90" \
      "$comorph""/microphysics/collision_ctl.F90" \
      "$comorph""/microphysics/microphysics_1.F90" \
      "$comorph""/microphysics/microphysics_2.F90" \
      "$comorph""/moist_proc/moist_proc.F90" \
      \
      "$comorph""/interface/standalone/tracer_source.F90" \
      \
      "$comorph""/plume_model/plume_model_diags_type_mod.F90" \
      "$comorph""/plume_model/sublevs_mod.F90" \
      "$comorph""/plume_model/calc_core_mean_ratio.F90" \
      "$comorph""/plume_model/calc_env_nsq.F90" \
      "$comorph""/plume_model/init_sublevs.F90" \
      "$comorph""/plume_model/entrain_fields.F90" \
      "$comorph""/plume_model/entdet_res_source.F90" \
      "$comorph""/plume_model/precip_res_source.F90" \
      "$comorph""/plume_model/solve_detrainment.F90" \
      "$comorph""/plume_model/wind_w_eqn.F90" \
      "$comorph""/plume_model/set_ent.F90" \
      "$comorph""/plume_model/set_det.F90" \
      "$comorph""/plume_model/calc_sat_height.F90" \
      "$comorph""/plume_model/set_par_cloudfrac.F90" \
      "$comorph""/plume_model/calc_mean_q_cl_with_core.F90" \
      "$comorph""/plume_model/momentum_eqn.F90" \
      "$comorph""/plume_model/parcel_dyn.F90" \
      "$comorph""/plume_model/update_edge_virt_temp.F90" \
      "$comorph""/plume_model/update_par_radius.F90" \
      "$comorph""/plume_model/interp_diag_conv_cloud_a.F90" \
      "$comorph""/plume_model/set_diag_conv_cloud_a.F90" \
      "$comorph""/plume_model/calc_cape.F90" \
      "$comorph""/plume_model/conv_level_step.F90" \
      "$comorph""/plume_model/normalise_integrals.F90" \
      \
      "$comorph""/plume_genesis/genesis_diags_type_mod.F90" \
      "$comorph""/plume_genesis/test_unstable.F90" \
      "$comorph""/plume_genesis/init_test.F90" \
      "$comorph""/plume_genesis/calc_fields_next.F90" \
      "$comorph""/plume_genesis/calc_turb_perts.F90" \
      "$comorph""/plume_genesis/calc_turb_parcel.F90" \
      "$comorph""/plume_genesis/set_par_fields.F90" \
      "$comorph""/plume_genesis/set_region_cond_fields.F90" \
      "$comorph""/plume_genesis/calc_env_region_tq_nb.F90" \
      "$comorph""/plume_genesis/calc_env_regions.F90" \
      "$comorph""/plume_genesis/calc_init_mass.F90" \
      "$comorph""/plume_genesis/calc_init_par_fields.F90" \
      "$comorph""/plume_genesis/region_parcel_calcs.F90" \
      "$comorph""/plume_genesis/calc_qss_forcing_init.F90" \
      "$comorph""/plume_genesis/cor_init_mass_liq_1.F90" \
      "$comorph""/plume_genesis/cfl_limit_init_mass.F90" \
      "$comorph""/plume_genesis/add_region_parcel.F90" \
      "$comorph""/plume_genesis/normalise_init_parcel.F90" \
      "$comorph""/plume_genesis/init_mass_moist_frac.F90" \
      \
      "$comorph""/control/draft_diags_type_mod.F90" \
      "$comorph""/control/init_diag_array.F90" \
      "$comorph""/control/comorph_diags_type_mod.F90" \
      "$comorph""/control/calc_turb_diags.F90" \
      "$comorph""/control/calc_sum_massflux.F90" \
      "$comorph""/control/calc_delta_tv.F90" \
      "$comorph""/control/par_gen_distinct_layers.F90" \
      "$comorph""/control/calc_pressure_incr_diag.F90" \
      "$comorph""/control/mass_rearrange_calc.F90" \
      "$comorph""/control/mass_rearrange.F90" \
      "$comorph""/control/add_res_source.F90" \
      "$comorph""/control/add_conv_cloud.F90" \
      "$comorph""/control/calc_diag_conv_cloud.F90" \
      "$comorph""/control/conv_genesis_ctl.F90" \
      "$comorph""/control/save_parcel_bl_top.F90" \
      "$comorph""/control/homog_conv_bl.F90" \
      "$comorph""/control/homog_conv_bl_ctl.F90" \
      "$comorph""/control/conv_sweep_compress.F90" \
      "$comorph""/control/conv_sweep_ctl.F90" \
      "$comorph""/control/conv_incr_ctl.F90" \
      "$comorph""/control/find_cmpr_any.F90" \
      "$comorph""/control/cfl_limit_indep.F90" \
      "$comorph""/control/cfl_limit_sum_ent.F90" \
      "$comorph""/control/cfl_limit_scaling.F90" \
      "$comorph""/control/apply_scaling.F90" \
      "$comorph""/control/conv_closure_ctl.F90" \
      "$comorph""/control/comorph_main.F90" \
      "$comorph""/control/comorph_ctl.F90" \
      \
      "$comorph""/unit_tests/set_test_profiles.F90" \
      "$comorph""/interface/standalone/standalone_test.F90" \
      "$comorph""/unit_tests/test_comorph.F90" \
1> std_out.txt  2> std_err.txt

echo "Done.  Compile output is in:"
echo "$1""/std_out.txt"
echo "$1""/std_err.txt"
