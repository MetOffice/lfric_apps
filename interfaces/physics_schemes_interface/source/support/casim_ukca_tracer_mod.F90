!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!-------------------------------------------------------------------------------
!> @brief Gather the GLOMAP aerosol fields into a UKCA tracer array column
!> @details The UM CASIM aerosol routines locate the GLOMAP mode numbers and
!!          component mass mixing ratios within the UKCA tracer array, using
!!          the indexing set up by ukca_set_conv_indices. This routine builds
!!          that array for a single LFRic column from the individual GLOMAP
!!          prognostic fields. Only the GLOMAP entries are filled, as the
!!          chemistry tracers are never looked at by CASIM.

module casim_ukca_tracer_mod

  use constants_mod, only: i_def, r_def, r_um

  implicit none

  private
  public :: casim_ukca_tracer_column

contains

  !> @brief Fill a column of the UKCA tracer array from the GLOMAP fields
  !> @param[in]  nlayers     Number of layers
  !> @param[in]  undf_wth    Number of unique degrees of freedom for potential
  !!                          temperature space
  !> @param[in]  base_wth    Dofmap entry for the bottom of the column in
  !!                          potential temperature space
  !> @param[in]  n_ait_sol   Soluble Aitken mode number mixing ratio
  !> @param[in]  ait_sol_su  Soluble Aitken mode H2SO4 mass mixing ratio
  !> @param[in]  ait_sol_bc  Soluble Aitken mode black carbon mass mixing ratio
  !> @param[in]  ait_sol_om  Soluble Aitken mode organic mass mixing ratio
  !> @param[in]  n_acc_sol   Soluble accumulation mode number mixing ratio
  !> @param[in]  acc_sol_su  Soluble accumulation mode H2SO4 mass mixing ratio
  !> @param[in]  acc_sol_bc  Soluble accumulation mode black carbon m.m.r.
  !> @param[in]  acc_sol_om  Soluble accumulation mode organic m.m.r.
  !> @param[in]  acc_sol_ss  Soluble accumulation mode sea salt m.m.r.
  !> @param[in]  n_cor_sol   Soluble coarse mode number mixing ratio
  !> @param[in]  cor_sol_su  Soluble coarse mode H2SO4 mass mixing ratio
  !> @param[in]  cor_sol_bc  Soluble coarse mode black carbon m.m.r.
  !> @param[in]  cor_sol_om  Soluble coarse mode organic m.m.r.
  !> @param[in]  cor_sol_ss  Soluble coarse mode sea salt m.m.r.
  !> @param[in]  n_ait_ins   Insoluble Aitken mode number mixing ratio
  !> @param[in]  ait_ins_bc  Insoluble Aitken mode black carbon m.m.r.
  !> @param[in]  ait_ins_om  Insoluble Aitken mode organic m.m.r.
  !> @param[in]  n_acc_ins   Insoluble accumulation mode number mixing ratio
  !> @param[in]  acc_ins_du  Insoluble accumulation mode dust m.m.r.
  !> @param[in]  n_cor_ins   Insoluble coarse mode number mixing ratio
  !> @param[in]  cor_ins_du  Insoluble coarse mode dust m.m.r.
  !> @param[out] tracer_ukca_col  UKCA tracer array for this column
  subroutine casim_ukca_tracer_column( nlayers, undf_wth, base_wth,            &
                                       n_ait_sol, ait_sol_su, ait_sol_bc,      &
                                       ait_sol_om,                             &
                                       n_acc_sol, acc_sol_su, acc_sol_bc,      &
                                       acc_sol_om, acc_sol_ss,                 &
                                       n_cor_sol, cor_sol_su, cor_sol_bc,      &
                                       cor_sol_om, cor_sol_ss,                 &
                                       n_ait_ins, ait_ins_bc, ait_ins_om,      &
                                       n_acc_ins, acc_ins_du,                  &
                                       n_cor_ins, cor_ins_du,                  &
                                       tracer_ukca_col )

    use nlsizes_namelist_mod, only: tr_ukca
    use um_ukca_init_mod,     only: tracer_names,                              &
                                    fldname_n_ait_sol,                         &
                                    fldname_ait_sol_su,                        &
                                    fldname_ait_sol_bc,                        &
                                    fldname_ait_sol_om,                        &
                                    fldname_n_acc_sol,                         &
                                    fldname_acc_sol_su,                        &
                                    fldname_acc_sol_bc,                        &
                                    fldname_acc_sol_om,                        &
                                    fldname_acc_sol_ss,                        &
                                    fldname_n_cor_sol,                         &
                                    fldname_cor_sol_su,                        &
                                    fldname_cor_sol_bc,                        &
                                    fldname_cor_sol_om,                        &
                                    fldname_cor_sol_ss,                        &
                                    fldname_n_ait_ins,                         &
                                    fldname_ait_ins_bc,                        &
                                    fldname_ait_ins_om,                        &
                                    fldname_n_acc_ins,                         &
                                    fldname_acc_ins_du,                        &
                                    fldname_n_cor_ins,                         &
                                    fldname_cor_ins_du

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in) :: base_wth

    real(kind=r_def), intent(in), dimension(undf_wth) :: n_ait_sol
    real(kind=r_def), intent(in), dimension(undf_wth) :: ait_sol_su
    real(kind=r_def), intent(in), dimension(undf_wth) :: ait_sol_bc
    real(kind=r_def), intent(in), dimension(undf_wth) :: ait_sol_om
    real(kind=r_def), intent(in), dimension(undf_wth) :: n_acc_sol
    real(kind=r_def), intent(in), dimension(undf_wth) :: acc_sol_su
    real(kind=r_def), intent(in), dimension(undf_wth) :: acc_sol_bc
    real(kind=r_def), intent(in), dimension(undf_wth) :: acc_sol_om
    real(kind=r_def), intent(in), dimension(undf_wth) :: acc_sol_ss
    real(kind=r_def), intent(in), dimension(undf_wth) :: n_cor_sol
    real(kind=r_def), intent(in), dimension(undf_wth) :: cor_sol_su
    real(kind=r_def), intent(in), dimension(undf_wth) :: cor_sol_bc
    real(kind=r_def), intent(in), dimension(undf_wth) :: cor_sol_om
    real(kind=r_def), intent(in), dimension(undf_wth) :: cor_sol_ss
    real(kind=r_def), intent(in), dimension(undf_wth) :: n_ait_ins
    real(kind=r_def), intent(in), dimension(undf_wth) :: ait_ins_bc
    real(kind=r_def), intent(in), dimension(undf_wth) :: ait_ins_om
    real(kind=r_def), intent(in), dimension(undf_wth) :: n_acc_ins
    real(kind=r_def), intent(in), dimension(undf_wth) :: acc_ins_du
    real(kind=r_def), intent(in), dimension(undf_wth) :: n_cor_ins
    real(kind=r_def), intent(in), dimension(undf_wth) :: cor_ins_du

    ! Shaped as the UM tracer array is, for a single column, so that it can be
    ! handed straight to the UM CASIM aerosol routines
    real(kind=r_um), intent(out) :: tracer_ukca_col(1, 1, 0:nlayers, tr_ukca)

    integer(kind=i_def) :: k, n

    ! The chemistry tracers are not used by CASIM, and neither are the
    ! soluble accumulation and coarse mode dust masses, for which LFRic has
    ! no prognostic. Everything therefore starts as zero and only the fields
    ! which LFRic holds are copied in.
    tracer_ukca_col(:,:,:,:) = 0.0_r_um

    do n = 1, tr_ukca
      select case (tracer_names(n))

      case (fldname_n_ait_sol)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( n_ait_sol(base_wth + k), r_um )
        end do
      case (fldname_ait_sol_su)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( ait_sol_su(base_wth + k), r_um )
        end do
      case (fldname_ait_sol_bc)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( ait_sol_bc(base_wth + k), r_um )
        end do
      case (fldname_ait_sol_om)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( ait_sol_om(base_wth + k), r_um )
        end do

      case (fldname_n_acc_sol)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( n_acc_sol(base_wth + k), r_um )
        end do
      case (fldname_acc_sol_su)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( acc_sol_su(base_wth + k), r_um )
        end do
      case (fldname_acc_sol_bc)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( acc_sol_bc(base_wth + k), r_um )
        end do
      case (fldname_acc_sol_om)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( acc_sol_om(base_wth + k), r_um )
        end do
      case (fldname_acc_sol_ss)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( acc_sol_ss(base_wth + k), r_um )
        end do

      case (fldname_n_cor_sol)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( n_cor_sol(base_wth + k), r_um )
        end do
      case (fldname_cor_sol_su)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( cor_sol_su(base_wth + k), r_um )
        end do
      case (fldname_cor_sol_bc)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( cor_sol_bc(base_wth + k), r_um )
        end do
      case (fldname_cor_sol_om)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( cor_sol_om(base_wth + k), r_um )
        end do
      case (fldname_cor_sol_ss)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( cor_sol_ss(base_wth + k), r_um )
        end do

      case (fldname_n_ait_ins)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( n_ait_ins(base_wth + k), r_um )
        end do
      case (fldname_ait_ins_bc)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( ait_ins_bc(base_wth + k), r_um )
        end do
      case (fldname_ait_ins_om)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( ait_ins_om(base_wth + k), r_um )
        end do

      case (fldname_n_acc_ins)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( n_acc_ins(base_wth + k), r_um )
        end do
      case (fldname_acc_ins_du)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( acc_ins_du(base_wth + k), r_um )
        end do

      case (fldname_n_cor_ins)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( n_cor_ins(base_wth + k), r_um )
        end do
      case (fldname_cor_ins_du)
        do k = 0, nlayers
          tracer_ukca_col(1,1,k,n) = real( cor_ins_du(base_wth + k), r_um )
        end do

      end select
    end do

  end subroutine casim_ukca_tracer_column

end module casim_ukca_tracer_mod
