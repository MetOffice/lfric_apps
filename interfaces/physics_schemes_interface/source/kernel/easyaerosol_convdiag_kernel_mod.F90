!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Calculate fields for easyaerosol diagnostics

module easyaerosol_convdiag_kernel_mod
use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_REAL, GH_READ,                       &
                             GH_INTEGER, GH_SCALAR,                            &
                             GH_WRITE, CELL_COLUMN,                            &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2,                        &
                             ANY_DISCONTINUOUS_SPACE_3,                        &
                             ANY_DISCONTINUOUS_SPACE_4

use kernel_mod,        only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: easyaerosol_convdiag_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                &
       ! n_band
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
       ! mode_dimen
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
       ! aer_prop
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),        &
       ! aer_prop_diag
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2),       &
       ! aer_mix_ratio
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3),        &
       ! rho_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4)         &
       /)

  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: easyaerosol_convdiag_code
end type

public :: easyaerosol_convdiag_code

contains

!> @brief         Convert aerosol absorption or scattering efficiency data
!>                (units m2 kg-1) from the socrates arrays into absorption or
!>                (scattering coefficient (units of m-1) and populate the
!>                diagnostic arrays. The data in the socrates arrays has to be
!>                multiplied by both air density and the aerosol mass mixing
!>                ratio.
!> @param[in]     nlayers          The number of layers
!> @param[in]     n_band           Number of spectral bands (sw or lw)
!> @param[in]     mode_dimen       Number of spectral bands (sw or lw)
!> @param[in]     aer_prop         aerosol property from easyaerosol
!> @param[in,out] aer_diag         aerosol property diagnostic
!> @param[in]     aer_mix_ratio    Aerosol mixing ratio for all modes
!> @param[in]     rho_wtheta       Air density on theta levels
!> @param[in]     ndf_aer          No. of DOFs per cell for aer space
!> @param[in]     undf_aer         No. unique of DOFs for aer space
!> @param[in]     map_aer          Dofmap for aer space column base cell
!> @param[in]     ndf_diag          No. of DOFs per cell for diag space
!> @param[in]     undf_diag         No. unique of DOFs for diag space
!> @param[in]     map_diag          Dofmap for diag space column base cell
!> @param[in]     ndf_mmr          No. of DOFs per cell for mmr space
!> @param[in]     undf_mmr         No. unique of DOFs for mmr space
!> @param[in]     map_mmr          Dofmap for mmr space column base cell
!> @param[in]     ndf_rho          No. DOFs per cell for 2D space
!> @param[in]     undf_rho         No. unique DOFs for rho_wtheta space
!> @param[in]     map_rho          Dofmap for rho_wtheta space
!
subroutine easyaerosol_convdiag_code( nlayers,                                 &
                                      n_band,                                  &
                                      mode_dimen,                              &
                                      aer_prop,                                &
                                      aer_diag,                                &
                                      aer_mix_ratio,                           &
                                      rho_wtheta,                              &
                                      ndf_aer,  undf_aer,  map_aer,            &
                                      ndf_diag,  undf_diag,  map_diag,         &
                                      ndf_mmr,  undf_mmr,  map_mmr,            &
                                      ndf_rho,  undf_rho,  map_rho)

  use constants_mod,       only: r_def, i_def, i_um

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: n_band
  integer(kind=i_def), intent(in) :: mode_dimen

  integer(kind=i_def), intent(in) :: ndf_aer
  integer(kind=i_def), intent(in) :: undf_aer
  integer(kind=i_def), dimension(ndf_aer), intent(in) :: map_aer

  integer(kind=i_def), intent(in) :: ndf_diag
  integer(kind=i_def), intent(in) :: undf_diag
  integer(kind=i_def), dimension(ndf_diag), intent(in) :: map_diag

  integer(kind=i_def), intent(in) :: ndf_mmr
  integer(kind=i_def), intent(in) :: undf_mmr
  integer(kind=i_def), dimension(ndf_mmr), intent(in) :: map_mmr

  integer(kind=i_def), intent(in) :: ndf_rho
  integer(kind=i_def), intent(in) :: undf_rho
  integer(kind=i_def), dimension(ndf_rho), intent(in) :: map_rho

  real(kind=r_def), intent(in),    dimension(undf_aer)  :: aer_prop
  real(kind=r_def), intent(inout), dimension(undf_diag) :: aer_diag
  real(kind=r_def), intent(in),    dimension(undf_mmr)  :: aer_mix_ratio
  real(kind=r_def), intent(in),    dimension(undf_rho)  :: rho_wtheta

  ! Loop counters / array index
  integer(i_um) :: k, i_band, i_bnd_mode, i_mode, i_aer, i_diag, i_rho, i_mmr

  i_bnd_mode = 0
  do i_band = 1, n_band
    do i_mode = 1, mode_dimen

      i_bnd_mode = i_bnd_mode + 1

      k = 0
      i_diag = map_aer(1) + ((i_bnd_mode-1) * (nlayers+1)) + k
      aer_diag (i_diag) = 0.0_r_def

      do k = 1, nlayers

        i_aer  = map_aer(1)  + ((i_bnd_mode-1) * (nlayers+1)) + k
        i_diag = map_diag(1) + ((i_bnd_mode-1) * (nlayers+1)) + k
        i_mmr  = map_mmr(1)  + ((i_mode-1)     * (nlayers+1)) + k
        i_rho  = map_rho(1)  + k

        aer_diag (i_diag) = aer_prop(i_aer)                                    &
                            * rho_wtheta(i_rho)                                &
                            * aer_mix_ratio(i_mmr)
      end do
    end do
  end do

end subroutine easyaerosol_convdiag_code

end module easyaerosol_convdiag_kernel_mod
