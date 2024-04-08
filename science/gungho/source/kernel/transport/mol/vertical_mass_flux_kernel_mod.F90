!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes vertical mass flux using an
!!        upwind reconstruction of a field on cell edges
module vertical_mass_flux_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_REAL,           &
                              GH_WRITE, GH_READ,           &
                              CELL_COLUMN,                 &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              ANY_W2
use constants_mod,     only : r_tran, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: vertical_mass_flux_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/              &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_W2), &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_W2), &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1) &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_mass_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vertical_mass_flux_code

contains

!> @brief Computes the vertical mass flux: wind*reconstruction.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] mass_flux      Vertical mass flux
!> @param[in]     wind           Wind field
!> @param[in]     reconstruction Tracer field reconstructed on cell edges
!> @param[in]     ndf_w2         Number of degrees of freedom per cell
!> @param[in]     undf_w2        Number of unique degrees of freedom for the wind field
!> @param[in]     map_w2         Dofmap for the cell at the base of the column
!> @param[in]     ndf_md         Number of degrees of freedom per cell
!> @param[in]     undf_md        Number of unique degrees of freedom for the
!!                               reconstructed field
!> @param[in]     map_md         Dofmap for the cell at the base of the column
subroutine vertical_mass_flux_code( nlayers,                &
                                    mass_flux,              &
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

  real(kind=r_tran), dimension(undf_md), intent(in)    :: reconstruction
  real(kind=r_tran), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_tran), dimension(undf_w2), intent(inout) :: mass_flux

  ! Internal variables
  integer(kind=i_def) :: k, df, ijp1, ijp2, offset
  real(kind=r_tran)   :: wgt

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

  mass_flux( map_w2(df) ) = 0.0_r_tran
  ! Vertical Flux, loop over edges ignoring the top and bottom face (where flux = 0)
  ! Reconstruction is stored on a layer first multidata field
  ! with index for face f = (0,1) = (Bottom, Top): map_md(1) + f*nlayers + k
  ! each cell contains the values for when it is the upwind cell for each edge
  ! so if u.n > 0 then we set the flux to be the value on the top edge from the cell below
  ! and if u.n < 0 then we set the flux to be the value on the bottom edge from this cell
  do k = 1, nlayers-1
    ! Value from bottom face of cell above edge
    ijp2 = map_md(1) + offset + nlayers + k-1
    ! Take value from top face of cell below edge
    ijp1 = map_md(1) + offset + k

    wgt = (0.5_r_tran - sign(0.5_r_tran, wind(map_w2(df)+k)))

    mass_flux( map_w2(df) + k) = wind(map_w2(df)+k)*( &
                                 wgt*reconstruction(ijp1) &
                               + (1.0_r_tran-wgt)*reconstruction(ijp2))

  end do
  mass_flux( map_w2(df) + nlayers ) = 0.0_r_tran

end subroutine vertical_mass_flux_code

end module vertical_mass_flux_kernel_mod
