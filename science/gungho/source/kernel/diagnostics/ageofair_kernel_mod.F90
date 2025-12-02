!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Zeros certain levels in age-of-air tracer field
!>
!> @details Zeros age-of-air tracer for levels nearest surface.
module ageofair_kernel_mod

  use argument_mod,      only : arg_type, func_type,     &
                                GH_FIELD, GH_REAL,       &
                                GH_READ, GH_WRITE,       &
                                CELL_COLUMN, GH_SCALAR,  &
                                GH_INTEGER
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W0, W2, W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: ageofair_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/             &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)   &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ageofair_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: ageofair_code

contains

!> @brief Zeros certain levels in age-of-air tracer field
!> @param[in]     nlayers  Number of layers
!> @param[in,out] ageofair Age-of-air tracer field
!> @param[in]     ndf_w3   Number of degrees of freedom per cell
!> @param[in]     undf_w3  Total number of degrees of freedom
!> @param[in]     map_w3   Dofmap for the cell at the base of the column for W3
subroutine ageofair_code( nlayers,                                &
                          ageofair,                               &
                          reset_level,                            &
                          ndf_w3,                                 &
                          undf_w3,                                &
                          map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3
  integer(kind=i_def), intent(in) :: undf_w3
  real(kind=r_def), dimension(undf_w3), intent(inout) :: ageofair
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in) :: reset_level

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def) :: df

  do k = 0, reset_level-1
    do df=1,ndf_w3
      ageofair( map_w3(df) + k ) =  0.0_r_def
    end do
  end do

end subroutine ageofair_code

end module ageofair_kernel_mod
