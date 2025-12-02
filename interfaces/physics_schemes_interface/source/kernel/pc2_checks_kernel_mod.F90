!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to cloud scheme.

module pc2_checks_kernel_mod

use argument_mod,      only: arg_type,          &
                             GH_FIELD, GH_REAL, &
                             GH_READ, GH_WRITE, &
                             CELL_COLUMN
use fs_continuity_mod, only: WTHETA
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: pc2_checks_kernel_type
  private
  type(arg_type) :: meta_args(17) = (/                &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! mv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! mi_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ms_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! bcf_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! theta_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! exner_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! dtheta_response
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! dmv_response_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! dmcl_response_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! dmci_response_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! dms_response_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! dcfl_response_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! dcff_response_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA)  & ! dbcf_response_wth
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: pc2_checks_code
end type

public :: pc2_checks_code

contains

!> @brief Interface to the pc2 checks code
!> @details Performs some consistency checks as part of the PC2
!>          as a result of combining cloud increments caculated in
!>          parallel without knowledge of each other.
!>          The PC2 cloud scheme is described in UMDP 30.
!> @param[in]     nlayers              Number of layers
!> @param[in]     mv_wth               Vapour mass mixing ratio
!> @param[in]     ml_wth               Liquid cloud mass mixing ratio
!> @param[in]     mi_wth               Ice cloud mass mixing ratio
!> @param[in]     ms_wth               Snow mass mixing ratio
!> @param[in]     cfl_wth              Liquid cloud fraction
!> @param[in]     cff_wth              Ice cloud fraction
!> @param[in]     bcf_wth              Bulk cloud fraction
!> @param[in]     theta_wth            Potential temperature field
!> @param[in]     exner_wth            Exner pressure in potential temperature space
!> @param[in,out] dtheta_response_wth  Change in theta
!> @param[in,out] dmv_response_wth     Change in water vapour
!> @param[in,out] dmcl_response_wth    Change in liquid water content
!> @param[in,out] dmci_response_wth    Change in ice    water content
!> @param[in,out] dms_response_wth     Change in snow   water content
!> @param[in,out] dcfl_response_wth    Change in liquid cloud fraction
!> @param[in,out] dcff_response_wth    Change in ice    cloud fraction
!> @param[in,out] dbcf_response_wth    Change in bulk   cloud fraction
!> @param[in]     ndf_wth              Number of degrees of freedom per cell for
!!                                      potential temperature space
!> @param[in]     undf_wth             Number of unique of degrees of freedom
!!                                      for potential temperature space
!> @param[in]     map_wth              Dofmap for the cell at the base of the column
!!                                      for potential temperature space

subroutine pc2_checks_code( nlayers,                   &
                                           ! Atm fields
                            mv_wth,                    &
                            ml_wth,                    &
                            mi_wth,                    &
                            ms_wth,                    &
                            cfl_wth,                   &
                            cff_wth,                   &
                            bcf_wth,                   &
                            theta_wth,                 &
                            exner_wth,                 &
                                           ! Responses
                            dtheta_response_wth,       &
                            dmv_response_wth,          &
                            dmcl_response_wth,         &
                            dmci_response_wth,         &
                            dms_response_wth,          &
                            dcfl_response_wth,         &
                            dcff_response_wth,         &
                            dbcf_response_wth,         &
                                           ! Other
                            ndf_wth, undf_wth, map_wth )

    use constants_mod, only: r_def, i_def, r_um, i_um

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use nlsizes_namelist_mod,       only: row_length, rows, model_levels
    use pc2_checks_mod,             only: pc2_checks
    use planet_constants_mod,       only: p_zero, kappa
    use gen_phys_inputs_mod,        only: l_mr_physics

    use free_tracers_inputs_mod,    only: n_wtrac
    use water_tracers_mod,          only: wtrac_type

    implicit none

    ! Arguments

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth

    real(kind=r_def), intent(in),  dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mi_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ms_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: bcf_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: theta_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_wth

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    ! The changes to the fields as a result
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dtheta_response_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv_response_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmcl_response_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmci_response_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dms_response_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_response_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcff_response_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dbcf_response_wth

    real(r_um), dimension(row_length,rows,model_levels) :: &
                  qv_work, qcl_work, qcf_work,             &
                  cfl_work, cff_work, bcf_work,            &
                  t_work, theta_work, pressure, qcf2_work

    integer(i_um) :: k

    integer(i_um), parameter :: nSCMDpkgs=15
    logical, parameter :: l_scmdiags(nscmdpkgs)=.false.

    ! Water tracer field which is not currently used but is required by
    ! UM routine
    type(wtrac_type), dimension(n_wtrac) :: wtrac

    !-----------------------------------------------------------------------
    ! Initialisation of prognostic variables and arrays
    !-----------------------------------------------------------------------

    do k = 1, model_levels

      t_work(1,1,k)   = theta_wth(map_wth(1) + k) *            &
                        exner_wth(map_wth(1) + k)

      ! Pressure at centre of theta levels
      pressure(1,1,k) = p_zero*(exner_wth(map_wth(1) + k))     &
                                **(1.0_r_um/kappa)

      ! Moist prognostics
      qv_work(1,1,k)  = mv_wth(map_wth(1) + k)
      qcl_work(1,1,k) = ml_wth(map_wth(1) + k)
      qcf_work(1,1,k) = ms_wth(map_wth(1) + k)
      qcf2_work(1,1,k) = mi_wth(map_wth(1) + k)

      ! Cast LFRic cloud fractions onto cloud fraction work arrays.
      bcf_work(1,1,k) = bcf_wth(map_wth(1) + k)
      cfl_work(1,1,k) = cfl_wth(map_wth(1) + k)
      cff_work(1,1,k) = cff_wth(map_wth(1) + k)

    end do

    call pc2_checks( pressure,                                 &
                     t_work, bcf_work, cfl_work, cff_work,     &
                     qv_work, qcl_work, qcf_work, l_mr_physics,&
                     row_length, rows, model_levels,           &
                     0_i_um, 0_i_um, 0_i_um, 0_i_um, qcf2_work,&
                     wtrac)

    ! Recast back to LFRic space
    do k = 1, model_levels
      ! *_work arrays have been updated

      ! New theta found from new temperature.
      theta_work(1,1,k) = t_work(1,1,k) /             &
                          exner_wth(map_wth(1) + k)

      ! All increments found from difference between the updated *_work values
      ! and the values that were intent in.
      dtheta_response_wth(map_wth(1) + k) = theta_work(1,1,k) &
                                          - theta_wth(map_wth(1) + k)
      !
      dmv_response_wth (map_wth(1)+k) = qv_work (1,1,k) - mv_wth (map_wth(1) + k)
      dmcl_response_wth(map_wth(1)+k) = qcl_work(1,1,k) - ml_wth (map_wth(1) + k)
      dmci_response_wth(map_wth(1)+k) = qcf2_work(1,1,k) - mi_wth (map_wth(1) + k)
      dms_response_wth(map_wth(1)+k) = qcf_work(1,1,k) - ms_wth (map_wth(1) + k)
      dcfl_response_wth(map_wth(1)+k) = cfl_work(1,1,k) - cfl_wth(map_wth(1) + k)
      dcff_response_wth(map_wth(1)+k) = cff_work(1,1,k) - cff_wth(map_wth(1) + k)
      dbcf_response_wth(map_wth(1)+k) = bcf_work(1,1,k) - bcf_wth(map_wth(1) + k)
    end do

    dtheta_response_wth(map_wth(1)+0) = dtheta_response_wth(map_wth(1)+1)
    dmv_response_wth   (map_wth(1)+0) = dmv_response_wth   (map_wth(1)+1)
    dmcl_response_wth  (map_wth(1)+0) = dmcl_response_wth  (map_wth(1)+1)
    dmci_response_wth  (map_wth(1)+0) = dmci_response_wth  (map_wth(1)+1)
    dms_response_wth   (map_wth(1)+0) = dms_response_wth   (map_wth(1)+1)
    dcfl_response_wth  (map_wth(1)+0) = dcfl_response_wth  (map_wth(1)+1)
    dcff_response_wth  (map_wth(1)+0) = dcff_response_wth  (map_wth(1)+1)
    dbcf_response_wth  (map_wth(1)+0) = dbcf_response_wth  (map_wth(1)+1)

end subroutine pc2_checks_code

end module pc2_checks_kernel_mod
