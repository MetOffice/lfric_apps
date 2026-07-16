!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint vertical mass flux using an
!!        upwind reconstruction of a field on cell edges.
module atl_vertical_mass_flux_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_REAL,           &
                              GH_READINC, GH_READ,         &
                              GH_INC,                      &
                              CELL_COLUMN,                 &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              ANY_W2
use constants_mod,     only : r_def, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: atl_vertical_mass_flux_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/              &
       arg_type(GH_FIELD,  GH_REAL, GH_READINC, ANY_W2), &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,    ANY_W2), &
       arg_type(GH_FIELD,  GH_REAL, GH_INC,     ANY_W2), &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,    ANY_DISCONTINUOUS_SPACE_1) &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: atl_vertical_mass_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: atl_vertical_mass_flux_code

contains

!> @brief Computes the vertical mass flux: wind*reconstruction.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] mass_flux      Vertical mass flux
!> @param[in]     ls_wind        Linear wind field
!> @param[in,out] wind           Perturbation wind field
!> @param[in]     reconstruction Tracer field reconstructed on cell edges
!> @param[in]     ndf_w2         Number of degrees of freedom per cell
!> @param[in]     undf_w2        Number of unique degrees of freedom for the wind field
!> @param[in]     map_w2         Dofmap for the cell at the base of the column
!> @param[in]     ndf_md         Number of degrees of freedom per cell
!> @param[in]     undf_md        Number of unique degrees of freedom for the
!!                               reconstructed field
!> @param[in]     map_md         Dofmap for the cell at the base of the column
subroutine atl_vertical_mass_flux_code( nlayers,                &
                                        mass_flux,              &
                                        ls_wind,                &
                                        wind,                   &
                                        reconstruction,         &
                                        ndf_w2,                 &
                                        undf_w2,                &
                                        map_w2,                 &
                                        ndf_md,                 &
                                        undf_md,                &
                                        map_md )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_md
  integer(kind=i_def), intent(in)                    :: undf_md
  integer(kind=i_def), dimension(ndf_md), intent(in) :: map_md
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(undf_md), intent(in)    :: reconstruction
  real(kind=r_def), dimension(undf_w2), intent(inout) :: wind
  real(kind=r_def), dimension(undf_w2), intent(in)    :: ls_wind
  real(kind=r_def), dimension(undf_w2), intent(inout) :: mass_flux

  ! Internal variables
  integer(kind=i_def) :: k, df, ijp1, ijp2, offset
  real(kind=r_def)    :: wgt

  ! df of dof on face
  if  ( ndf_w2 == 2 ) then
    ! W2v space, reconstruction has ndata=2
    df = 1
    offset = 0
  else
    ! W2 space, reconstruction has ndata=6
    df = 5
    offset = 4*nlayers
  end if

  mass_flux( map_w2(df) + nlayers ) = 0.0_r_def
  do k = nlayers - 1, 1, -1
    ! Value from bottom face of cell above edge
    ijp2 = map_md(1) + offset + nlayers + k-1
    ! Take value from top face of cell below edge
    ijp1 = map_md(1) + offset + k

    wgt = (0.5_r_def - sign(0.5_r_def, ls_wind(map_w2(df)+k)))
    wind(k + map_w2(df)) = wind(k + map_w2(df)) &
                         + mass_flux(map_w2(df) + k)*( &
                           wgt*reconstruction(ijp1) &
                           + (1.0_r_def-wgt)*reconstruction(ijp2))
                
  end do
  mass_flux( map_w2(df) ) = 0.0_r_def

end subroutine atl_vertical_mass_flux_code

end module atl_vertical_mass_flux_kernel_mod
