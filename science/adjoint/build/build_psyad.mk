##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Wrapper script to build the required PSyAD targets.
# The PSyAD kernels are built using three stages.
# 1) Pre-patch: copies (and patches) tangent linear kernels from their base dir to PSYAD_WDIR.
# 2) PSyAD: generates adjoint kernels and adjoint test algorithms from the pre-patch stage tangent linear kernels.
# 3) Post-patch: copies (and patches) adjoint kernels and adjoint test algorithms from PSYAD_WDIR to WORKING_DIR.

# Kernels which need adjointing...
#
KERNEL_LIST = kernel/linear_physics/stabilise_bl_u_kernel_mod.f90 \
              kernel/solver/apply_mixed_lu_operator_kernel_mod.f90 \
              kernel/solver/apply_mixed_operator_kernel_mod.f90 \
              kernel/solver/opt_apply_variable_hx_kernel_mod.f90 \
              kernel/solver/apply_elim_mixed_lp_operator_kernel_mod.f90

# Variables per kernel...
#
ACTIVE_tl_poly_advective_kernel_mod            := advective dtdx dtdy v u tracer wind
ACTIVE_tl_poly1d_vert_adv_kernel_mod           := advective wind dpdz tracer
ACTIVE_tl_vorticity_advection_kernel_mod       := r_u wind vorticity vorticity_at_quad \
                                                  u_at_quad j_vorticity vorticity_term \
                                                  res_dot_product cross_product1 cross_product2 mul2
ACTIVE_stabilise_bl_u_kernel_mod               := u_stabilised u_initial u_final
ACTIVE_apply_mixed_lu_operator_kernel_mod      := wind theta exner lhs_u lhs_t
ACTIVE_apply_mixed_operator_kernel_mod         := u_e t_col lhs_p lhs_w lhs_uv exner wind_w wind_uv
ACTIVE_opt_apply_variable_hx_kernel_mod        := lhs x pressure rhs_p \
                                                  div_u t_e t_e1_vec t_e2_vec
ACTIVE_apply_elim_mixed_lp_operator_kernel_mod := theta exner u lhs_exner \
                                                  lhs_e m3e_pe p3t_te q32_ue \
                                                  p_e t_e u_e
ACTIVE_combine_w2_field_kernel_mod             := uvw w uv
ACTIVE_w2_to_w1_projection_kernel_mod          := v_w1 u_w2 vu res_dot_product wind
ACTIVE_sample_field_kernel_mod                 := field_1 field_2 f_at_node
ACTIVE_sample_flux_kernel_mod                  := flux u
ACTIVE_split_w2_field_kernel_mod               := uvw w uv
ACTIVE_strong_curl_kernel_mod                  := xi res_dot_product curl_u u
ACTIVE_sci_average_w2b_to_w2_kernel_mod        := field_w2 field_w2_broken
ACTIVE_sci_extract_w_kernel_mod                := velocity_w2v u_physics
ACTIVE_sci_combine_multidata_field_kernel_mod  := field1_in field2_in field_out
ACTIVE_tl_horizontal_mass_flux_kernel_mod      := mass_flux wind
ACTIVE_tl_vertical_mass_flux_kernel_mod        := mass_flux wind
ACTIVE_w3v_advective_update_kernel_mod         := advective_increment tracer dtdz t_U t_D
ACTIVE_tl_w3v_advective_update_kernel_mod      := advective_increment wind w
ACTIVE_horizontal_mass_flux_kernel_mod         := mass_flux reconstruction
ACTIVE_vertical_mass_flux_kernel_mod           := mass_flux reconstruction

.SECONDEXPANSION:

ADJOINT_KERNEL_LIST = $(patsubst %_kernel_mod.f90,%_kernel_mod_adjoint.f90,$(KERNEL_LIST))

# Converting kernel filenames into test algorithm filenames is somewhat
# involved.
#
kernel_type = $(firstword $(subst _, ,$(notdir $(kernel))))
kernel_name = $(patsubst tl_%_kernel_mod.f90,%,$(notdir $(kernel)))
TEST_LIST = $(foreach kernel,$(KERNEL_LIST),$(patsubst kernel/%,algorithm/%,$(dir $(kernel)))$(if $(filter tl,$(kernel_type)),atlt_,adjt_)$(kernel_name)_alg_mod.x90)

#build_adjoint: $(error $(addprefix $(WORKING_DIR)/,$(ADJOINT_KERNEL_LIST))) $(addprefix $(WORKING_DIR)/,$(TEST_LIST))
build_adjoint: $(WORKING_DIR)/kernel/transport/mol/tl_poly1d_vert_adv_kernel_mod_adjoint.f90 $(WORKING_DIR)/kernel/fem/strong_curl_kernel_mod_adjoint.f90

# Generate kernel adjoints and test algorithms.
#
# PSyAd requires test algorithm names to follow certain patterns. This is
# complicated to generate.
#
# kernel_name   - Just the name of the kernel, without any of the trailing
#                 bits.
# linear_kernel - "tl_" or empty.
# alg_name      - "atlt_..." if "tl_" appears, otherwise "adjt_..."
#
# If rule is triggered by kernel target:
#
$(PSYAD_WDIR)/kernel/%_kernel_mod_adjoint.f90 : kernel_name = $(notdir $*)
$(PSYAD_WDIR)/kernel/%_kernel_mod_adjoint.f90 : kernel_file = $@
$(PSYAD_WDIR)/kernel/%_kernel_mod_adjoint.f90 : linear_kernel = $(filter tl,$(subst _,$(SPACE),$(kernel_name)))
$(PSYAD_WDIR)/kernel/%_kernel_mod_adjoint.f90 : algorithm_name = $(if $(linear_kernel),$(patsubst tl_%,atlt_%,$(kernel_name)),$(addprefix adjt_,$(kernel_name)))
$(PSYAD_WDIR)/kernel/%_kernel_mod_adjoint.f90 : algorithm_file = $(PSYAD_WDIR)/algorithm/$(dir $*)$(algorithm_name)_alg_mod.x90
#
# If rule is triggered by algorithm target:
#
$(PSYAD_WDIR)/algorithm/%_alg_mod.x90: algorithm_name = $(notdir $*)
$(PSYAD_WDIR)/algorithm/%_alg_mod.x90: algorithm_file = $@
$(PSYAD_WDIR)/algorithm/%_alg_mod.x90: linear_algorithm = $(filter atlt,$(subst _,$(SPACE),$(algorithm_name)))
$(PSYAD_WDIR)/algorithm/%_alg_mod.x90: kernel_name = $(if $(linear_algorithm),$(patsubst atlt_%,tl_%,$(algorithm_name)),$(patsubst adjt_%,%,$(algorithm_name)))
$(PSYAD_WDIR)/algorithm/%_alg_mod.x90: kernel_file = $(PSYAD_WDIR)/kernel/$(dir $*)$(kernel_name)_kernel_mod_adjoint.f90
$(PSYAD_WDIR)/kernel/%_kernel_mod_adjoint.f90 $(PSYAD_WDIR)/algorithm/%_alg_mod.x90 \
    : $$(PSYAD_WDIR)/kernel/%_kernel_mod.f90 \
    | $$(PSYAD_WDIR)/kernel/$$(dir $$*) $$(PSYAD_WDIR)/algorithm/$$(dir $$*)
	$(call MESSAGE,PSyAd,$<)
	psyad -api lfric \
	      -oad $(kernel_file) \
	      -otest $(algorithm_file) \
	      -c $(PSYAD_CONFIG_FILE) \
	      -a $(ACTIVE_$(basename $(notdir $<))) \
	      -- $<

# If a patch exists, patch it to PSyAd workspace (pre-PSyAd)
#
$(PSYAD_WDIR)/%_kernel_mod.f90: $$(PATCH_DIR)/kernel/$$(notdir $$*)_kernel_mod.patch $(WORKING_DIR)/%_kernel_mod.f90 \
    | $$(dir $$@)
	$(call MESSAGE,Pre-patch,$<)
	patch $(WORKING_DIR)/$*_kernel_mod.f90 $< -o $@

# If no patch exists, just copy to PSyAd workspace (pre-PSyAd)
#
$(PSYAD_WDIR)/%_kernel_mod.f90: $(WORKING_DIR)/%_kernel_mod.f90 \
    | $$(dir $$@)
	$(call MESSAGE,Pre-copy,$<)
	cp $< $@


# If a patch exists, patch kernel to the buildspace (post-PsyAd)
#
$(WORKING_DIR)/%_kernel_mod_adjoint.f90: $$(PATCH_DIR)/kernel/$$(notdir $$*)_kernel_mod_adjoint.patch \
    $(PSYAD_WDIR)/%_kernel_mod_adjoint.f90 \
    | $$(dir $$@)
	$(call MESSAGE,Post-patch,$<)
	patch $(PSYAD_WDIR)/$*_kernel_mod_adjoint.f90 $< -o $@

# If no patch exists, just copy kernel to buildspace (post-PsyAd)
#
$(WORKING_DIR)/%_kernel_mod_adjoint.f90: $(PSYAD_WDIR)/%_kernel_mod_adjoint.f90 \
    | $$(dir $$@)
	$(call MESSAGE,Post-copy,$<)
	cp $< $@

# If a patch exists, patch algorithm to the buildspace (post-PSyAd)
#
# Building the adjoint kernel also generates the test algorithm so the former
# is a prerequisite of this rule
#
$(WORKING_DIR)/%_alg_mod.x90: $$(PATCH_DIR)/algorithm/$$(notdir $$*)_alg_mod.patch \
    $$(PSYAD_WDIR)/$$(dir $$*)$$(patsubst algorithm/adjt_%,%,$$(patsubst atlt_%,kernel/tl_%,$$(notdir $$*)))_kernel_mod_adjoint.f90 \
    | $$(dir $@))
	$(call MESSAGE,Patch,$<)
	patch $(PSYAD_WDIR)/$*_alg_mod.x90 $< -o $@

# If no patch exists, just copy algorithm to buildspace (post-PSyAd)
#
$(WORKING_DIR)/%_alg_mod.x90: $(PSYAD_WDIR)/%_alg_mod.f90 \
    | $$(dir $@))
	$(call MESSAGE,Copy,$<)
	cp $< $@


source_dirs := $(shell find $(WORKING_DIR) -type d -printf %P/\\n)
kernel_wdir := $(addprefix $(PSYAD_WDIR)/,$(source_dirs))
algorithm_wdir := $(patsubst kernel/%,$(PSYAD_WDIR)/algorithm/%,$(filter kernel/%,$(source_dirs)))
$(sort $(kernel_wdir) $(algorithm_wdir)):
	mkdir -p $@


include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
