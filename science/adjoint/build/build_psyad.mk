##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Orchestrates the building of adjoint kernels and test algorithms.
#
# Expects the following macros to be defined:
# BUILD_ADJ_TESTS   - Set if the test algorithms are needed
# PATCH_DIR         - Directory containing pre and post PSyAd patches
# PSYAD_CONFIG_FILE - Configuration file for PSyAd
# PSYAD_WDIR        - Extra scratch space for this process
# WORKING_DIR       - Build's scratch space
##############################################################################
#
# Kernels which need adjointing...
#
KERNEL_LIST = kernel/linear_physics/stabilise_bl_u_kernel_mod.f90 \
              kernel/solver/apply_mixed_lu_operator_kernel_mod.f90 \
              kernel/solver/apply_mixed_operator_kernel_mod.f90 \
              kernel/solver/opt_apply_variable_hx_kernel_mod.f90 \
              kernel/solver/apply_elim_mixed_lp_operator_kernel_mod.f90

# Variables per kernel...
#
ACTIVE_stabilise_bl_u_kernel_mod               := u_stabilised u_initial u_final
ACTIVE_apply_mixed_lu_operator_kernel_mod      := wind theta exner lhs_u lhs_t
ACTIVE_apply_mixed_operator_kernel_mod         := u_e t_col lhs_p lhs_w lhs_uv exner wind_w wind_uv
ACTIVE_opt_apply_variable_hx_kernel_mod        := lhs x pressure rhs_p \
                                                  div_u t_e t_e1_vec t_e2_vec
ACTIVE_apply_elim_mixed_lp_operator_kernel_mod := theta exner u lhs_exner \
                                                  lhs_e m3e_pe p3t_te q32_ue \
                                                  p_e t_e u_e

do_psyad: $(if $(BUILD_ADJ_TESTS),$(WORKING_DIR)/driver/gen_adj_kernel_tests_mod.f90) \
          $(foreach kernel,$(KERNEL_LIST),$(WORKING_DIR)/$(kernel).adjoint)

$(WORKING_DIR)/driver/gen_adj_kernel_tests_mod.f90: | $(WORKING_DIR)/driver
	$(call MESSAGE,PSyAd test driver,$@)
	python $(ADJOINT_BUILD)/psyad_driver $(ADJOINT_BUILD)/gen_adj_kernel_tests_mod.f90.jinja \
                                         $@ $(notdir $(KERNEL_LIST))

%.f90.adjoint: %.f90 | $(PSYAD_WDIR)
	$(call MESSAGE,PSyAd,$<)
	$(ADJOINT_BUILD)/psyad_wrapper $(if $(BUILD_ADJ_TESTS),--algorithm-dir=$(PSYAD_WDIR)) \
                                   $(PSYAD_CONFIG_FILE) \
                                   $< \
                                   $(PATCH_DIR) \
                                   $(WORKING_DIR) \
                                   $(PSYAD_WDIR) \
                                   $(ACTIVE_$(patsubst %.f90.adjoint,%,$(notdir $*)))

$(PSYAD_WDIR) $(WORKING_DIR)/adjoint $(WORKING_DIR)/driver:
	$(call MESSAGE,Create,$@)
	mkdir -p $@
