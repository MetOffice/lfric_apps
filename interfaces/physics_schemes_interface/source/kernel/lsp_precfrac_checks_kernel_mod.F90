!------------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief Interface to prognostic precip fraction checks routine

module lsp_precfrac_checks_kernel_mod

use argument_mod,      only: arg_type,          &
                             GH_FIELD, GH_REAL, &
                             GH_READ, GH_WRITE, &
                             CELL_COLUMN
use fs_continuity_mod, only: WTHETA, W3
use kernel_mod,        only: kernel_type

implicit none

private

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: lsp_precfrac_checks_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/                                          &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),     & ! exner_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! mv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! mi_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ms_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! mr_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! mg_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA)  & ! precfrac_wth
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: lsp_precfrac_checks_code
end type

public :: lsp_precfrac_checks_code

contains

!> @brief Interface to prognostic precip fraction checks routine
!> @details Performs sanity-checks on the prognostic fractional coverage of
!>          falling rain and graupel.  Ensure the in-precip-shaft precip
!>          concentration ( m_rain [+m_graupel] )/precfrac
!>          has a plausible value, and ensure 0 <= precfrac <= 1.
!> @param[in]     nlayers              Number of layers
!> @param[in]     exner_w3             Exner pressure in w3 space
!> @param[in]     mv_wth               Vapour mass mixing ratio
!> @param[in]     ml_wth               Liquid cloud mass mixing ratio
!> @param[in]     mi_wth               Ice cloud mass mixing ratio
!> @param[in]     ms_wth               Snow mass mixing ratio
!> @param[in]     mr_wth               Rain mass mixing ratio
!> @param[in]     mg_wth               Graupel mass mixing ratio
!> @param[in,out] precfrac_wth         Prognostic precip fraction
!> @param[in]     ndf_w3               Degrees of freedom per cell; w3 space
!> @param[in]     undf_w3              Unique of degrees of freedom; w3 space
!> @param[in]     map_w3               Dofmap  celofl at column base; w3 space
!> @param[in]     ndf_wth              Degrees of freedom per cell; wth space
!> @param[in]     undf_wth             Unique of degrees of freedom; wth space
!> @param[in]     map_wth              Dofmap  celofl at column base; wth space

subroutine lsp_precfrac_checks_code( nlayers,                                  &
                                     exner_w3,                                 &
                                     mv_wth, ml_wth, mi_wth,                   &
                                     ms_wth, mr_wth, mg_wth,                   &
                                     precfrac_wth,                             &
                                     ndf_w3,  undf_w3,  map_w3,                &
                                     ndf_wth, undf_wth, map_wth )

    use constants_mod,           only: r_def, i_def, r_um, i_um

    ! UM modules
    use nlsizes_namelist_mod,    only: row_length, rows, model_levels
    use atm_fields_bounds_mod,   only: pdims
    use planet_constants_mod,    only: p_zero, kappa
    use lsp_precfrac_checks_mod, only: lsp_precfrac_checks

    implicit none

    ! Arguments

    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    integer(kind=i_def), intent(in) :: nlayers

    real(kind=r_def), intent(in),    dimension(undf_w3)  :: exner_w3
    real(kind=r_def), intent(in),    dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: mi_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: ms_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: mr_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: mg_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: precfrac_wth

    ! Arguments converted to r_um precision for input to lsp_precfrac_checks
    real(r_um), dimension(row_length,rows,model_levels) ::                     &
                  p_rho_levels, q, qcl, qcf, qcf2, qrain, qgraup, precfrac

    ! Loop counter
    integer(i_um) :: k


    ! Copy / convert fields for input to lsp_precfrac_checks
    do k = 1, nlayers

      ! Pressure at layer boundaries
      p_rho_levels(1,1,k) = p_zero * exner_w3(map_w3(1)+k-1)**(1.0_r_um/kappa)

      ! Water species masses
      q(1,1,k)      = mv_wth(map_wth(1)+k)
      qcl(1,1,k)    = ml_wth(map_wth(1)+k)
      qcf2(1,1,k)   = mi_wth(map_wth(1)+k)
      qcf(1,1,k)    = ms_wth(map_wth(1)+k)
      qrain(1,1,k)  = mr_wth(map_wth(1)+k)
      qgraup(1,1,k) = mg_wth(map_wth(1)+k)

      ! Input value of prognostic precip fraction
      precfrac(1,1,k) = precfrac_wth(map_wth(1)+k)

    end do

    ! Call the checking routine
    CALL lsp_precfrac_checks( pdims, p_rho_levels,                             &
                              q, qcl, qcf, qcf2, qrain, qgraup,                &
                              precfrac )

    ! Recast back to LFRic space
    do k = 1, nlayers
      ! Updated value of prognostic precip fraction
      precfrac_wth(map_wth(1)+k) = precfrac(1,1,k)
    end do
    precfrac_wth(map_wth(1)+0) = precfrac_wth(map_wth(1)+1)

end subroutine lsp_precfrac_checks_code

end module lsp_precfrac_checks_kernel_mod
