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
all: export ACTIVE_sci_average_w2b_to_w2_kernel_mod    := field_w2 field_w2_broken
all: export ACTIVE_sci_extract_w_kernel_mod             := velocity_w2v u_physics
all: export ACTIVE_sci_combine_multidata_field_kernel_mod := field1_in field2_in field_out
all: export ACTIVE_tl_horizontal_mass_flux_kernel_mod     := mass_flux wind
all: export ACTIVE_tl_vertical_mass_flux_kernel_mod       := mass_flux wind
all: export ACTIVE_w3v_advective_update_kernel_mod        := advective_increment tracer dtdz t_U t_D
all: export ACTIVE_tl_w3v_advective_update_kernel_mod     := advective_increment wind w
all: export ACTIVE_horizontal_mass_flux_kernel_mod        := mass_flux reconstruction
all: export ACTIVE_vertical_mass_flux_kernel_mod          := mass_flux reconstruction
