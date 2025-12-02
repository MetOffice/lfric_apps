##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Active variables
all: export ACTIVE_tl_poly_advective_kernel_mod            := advective dtdx dtdy v u tracer wind
all: export ACTIVE_tl_poly1d_vert_adv_kernel_mod           := advective wind dpdz tracer
all: export ACTIVE_tl_vorticity_advection_kernel_mod       := r_u wind vorticity vorticity_at_quad \
                                                              u_at_quad j_vorticity vorticity_term \
                                                              res_dot_product cross_product1 cross_product2 mul2
all: export ACTIVE_tl_sample_eos_pressure_kernel_mod       := exner tmp_exner \
                                                              theta_vd_cell theta_vd_e theta \
                                                              moist_dyn_gas rho_cell rho_e rho
all: export ACTIVE_tl_moist_dyn_gas_kernel_mod             := moist_dyn_gas mr_v mr_v_at_dof
all: export ACTIVE_tl_hydrostatic_kernel_mod               := r_u exner theta \
                                                              moist_dyn_gas moist_dyn_tot \
                                                              exner_e theta_v_e grad_theta_v_at_quad \
                                                              exner_at_quad theta_v_at_quad \
                                                              grad_term res_dot_product
all: export ACTIVE_tl_moist_dyn_mass_kernel_mod            := moist_dyn_tot mr_v mr_cl mr_r \
                                                              mr_s mr_g mr_ci \
                                                              mr_v_at_dof mr_cl_at_dof \
                                                              mr_r_at_dof mr_s_at_dof \
                                                              mr_g_at_dof mr_ci_at_dof
all: export ACTIVE_tl_kinetic_energy_gradient_kernel_mod   := r_u u ru_e u_e \
                                                              u_at_quad ke_at_quad \
                                                              res_dot_product mul2
all: export ACTIVE_tl_project_eos_pressure_kernel_mod      := exner exner_e exner_at_quad \
                                                              tmp_exner r_exner \
                                                              rho rho_e rho_at_quad theta \
                                                              theta_vd_e theta_vd_at_quad \
                                                              moist_dyn_gas
all: export ACTIVE_tl_rhs_project_eos_kernel_mod           := rhs_eos exner exner_e exner_quad \
                                                              rho rho_e rho_quad \
                                                              theta theta_vd_e theta_vd_quad \
                                                              moist_dyn_gas eos
all: export ACTIVE_tl_rhs_sample_eos_kernel_mod            := exner exner_e exner_cell \
                                                              rhs_eos rho rho_e rho_cell \
                                                              theta theta_vd_e theta_vd_cell moist_dyn_gas
all: export ACTIVE_stabilise_bl_u_kernel_mod               := u_stabilised u_initial u_final
all: export ACTIVE_apply_mixed_lu_operator_kernel_mod      := wind theta exner lhs_u lhs_t
all: export ACTIVE_apply_mixed_operator_kernel_mod         := u_e t_col lhs_p lhs_w lhs_uv exner wind_w wind_uv
all: export ACTIVE_opt_apply_variable_hx_kernel_mod        := lhs x pressure rhs_p \
                                                              div_u t_e t_e1_vec t_e2_vec
all: export ACTIVE_apply_elim_mixed_lp_operator_kernel_mod := theta exner u lhs_exner \
                                                              lhs_e m3e_pe p3t_te q32_ue \
                                                              p_e t_e u_e
all: export ACTIVE_combine_w2_field_kernel_mod         := uvw w uv
all: export ACTIVE_w2_to_w1_projection_kernel_mod      := v_w1 u_w2 vu res_dot_product wind
all: export ACTIVE_sample_field_kernel_mod             := field_1 field_2 f_at_node
all: export ACTIVE_sample_flux_kernel_mod              := flux u
all: export ACTIVE_split_w2_field_kernel_mod           := uvw w uv
all: export ACTIVE_strong_curl_kernel_mod              := xi res_dot_product curl_u u
all: export ACTIVE_matrix_vector_kernel_mod            := lhs x
all: export ACTIVE_dg_inc_matrix_vector_kernel_mod     := lhs x
all: export ACTIVE_transpose_matrix_vector_kernel_mod  := lhs lhs_e x x_e
all: export ACTIVE_dg_matrix_vector_kernel_mod         := lhs x
all: export ACTIVE_sci_average_w2b_to_w2_kernel_mod    := field_w2 field_w2_broken
all: export ACTIVE_sci_extract_w_kernel_mod             := velocity_w2v u_physics
all: export ACTIVE_sci_combine_multidata_field_kernel_mod := field1_in field2_in field_out
all: export ACTIVE_tl_horizontal_mass_flux_kernel_mod     := mass_flux wind
all: export ACTIVE_tl_vertical_mass_flux_kernel_mod       := mass_flux wind
all: export ACTIVE_w3v_advective_update_kernel_mod        := advective_increment tracer dtdz t_U t_D
all: export ACTIVE_tl_w3v_advective_update_kernel_mod     := advective_increment wind w
all: export ACTIVE_horizontal_mass_flux_kernel_mod        := mass_flux reconstruction
all: export ACTIVE_vertical_mass_flux_kernel_mod          := mass_flux reconstruction
