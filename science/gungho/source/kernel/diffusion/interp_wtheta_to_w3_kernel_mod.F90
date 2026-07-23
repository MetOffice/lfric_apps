!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Obtains a W3 field from a Wtheta field by sampling at W3 DoFs.
!>
!> @details Takes a Wtheta field and returns a W3 field, whose values are
!>          obtained by interpolating the Wtheta values to the W3 DoFs.
!>
!>          Only designed to work with the lowest order finite elements.
!>
module interp_wtheta_to_w3_kernel_mod

  use argument_mod,            only : arg_type,          &
                                      GH_FIELD, GH_REAL, &
                                      GH_WRITE, GH_READ, &
                                      CELL_COLUMN
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : Wtheta, W3
  use kernel_mod,              only : kernel_type

  implicit none

  private

  type, public, extends(kernel_type) :: interp_wtheta_to_w3_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),     &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),     &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: interp_wtheta_to_w3_code
  end type

public :: interp_wtheta_to_w3_code

contains

!> @brief Obtains a W3 field from a Wtheta field by sampling at W3 DoFs.
!> @param[in]     nlayers    Number of layers
!> @param[in,out] field_w3   Output W3 field
!> @param[in]     field_wth  Input Wtheta field
!> @param[in]     height_w3  Height above surface of W3 DoFs
!> @param[in]     height_wt  Height above surface of Wtheta DoFs
!> @param[in]     ndf_w3     Number of degrees of freedom per cell for W3
!> @param[in]     undf_w3    Total number of degrees of freedom for W3
!> @param[in]     map_w3     Dofmap for the cell at the base of the column for W3
!> @param[in]     ndf_wt     Number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wt    Total number of degrees of freedom for Wtheta
!> @param[in]     map_wt     Dofmap for the cell at the base of the column for Wtheta

subroutine interp_wtheta_to_w3_code( nlayers,   &
                                     field_w3,  &
                                     field_wt,  &
                                     height_w3, &
                                     height_wt, &
                                     ndf_w3,    &
                                     undf_w3,   &
                                     map_w3,    &
                                     ndf_wt,    &
                                     undf_wt,   &
                                     map_wt     )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: undf_w3, ndf_w3
  integer(kind=i_def),                     intent(in)    :: undf_wt, ndf_wt

  integer(kind=i_def), dimension(ndf_w3),  intent(in)    :: map_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in)    :: map_wt

  real(kind=r_def), dimension(undf_w3), intent(inout) :: field_w3
  real(kind=r_def), dimension(undf_wt), intent(in)    :: field_wt
  real(kind=r_def), dimension(undf_w3), intent(in)    :: height_w3
  real(kind=r_def), dimension(undf_wt), intent(in)    :: height_wt

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_def) :: weight_lower, weight_upper, weight_denom

  do k = 0, nlayers - 1
    weight_denom = height_wt(map_wt(1)+k+1) - height_wt(map_wt(1)+k)
    weight_upper = height_w3(map_w3(1)+k) - height_wt(map_wt(1)+k)
    weight_lower = height_wt(map_wt(1)+k+1) - height_w3(map_w3(1)+k)
    field_w3(map_w3(1)+k) =                                                    &
                weight_upper / weight_denom * field_wt(map_wt(1)+k+1)          &
              + weight_lower / weight_denom * field_wt(map_wt(1)+k)
  end do

end subroutine interp_wtheta_to_w3_code

end module interp_wtheta_to_w3_kernel_mod
