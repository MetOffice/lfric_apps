!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Mitigation to linear model boundary layer instability.

module stabilise_bl_u_kernel_mod

  use argument_mod,      only : arg_type, CELL_COLUMN, GH_FIELD, GH_SCALAR, &
                                GH_REAL, GH_INTEGER, GH_READ, GH_INC
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: stabilise_bl_u_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/             &
         arg_type(GH_FIELD,  GH_REAL, GH_INC,  W2), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, W2), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, W2), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ ), &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ )  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: stabilise_bl_u_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  public :: stabilise_bl_u_code

contains

!> @brief Mitigation to linear model boundary layer instability in column of u.
!> @details Sets u_stabilised(k) equal to:
!!            alpha(k) * u_initial(k) + (1 - alpha(k)) * u_final(k)
!!          where alpha(k) lies between 0 and 1, is largest near the surface,
!!          and decreases linearly to zero over a certain number of levels.
!> @param[in]     nlayers               Number of layers
!> @param[in,out] u_stabilised          Column of u to be stabilised
!> @param[in]     u_initial             Column of u at start of timestep
!> @param[in]     u_final               Column of u at end of timestep
!> @param[in]     n_levels_to_stabilise Number of levels (from base) to modify
!> @param[in]     max_stabilisation     Strength of stabilisation, range [0, 1]
!> @param[in]     ndf1                  Number of degrees of freedom per cell
!!                                      for the output field
!> @param[in]     undf1                 Unique number of degrees of freedom for
!!                                      the output field
!> @param[in]     map1                  Dofmap for the cell at the base of the
!!                                      column for the output field
subroutine stabilise_bl_u_code(nlayers, u_stabilised, u_initial, u_final, &
                               n_levels_to_stabilise, max_stabilisation,  &
                               ndf1, undf1, map1)

  implicit none

  ! Arguments
  integer(kind=i_def),                   intent(in)    :: undf1, ndf1
  integer(kind=i_def),                   intent(in)    :: nlayers
  real(kind=r_def),    dimension(undf1), intent(inout) :: u_stabilised
  real(kind=r_def),    dimension(undf1), intent(in)    :: u_initial, u_final
  integer(kind=i_def),                   intent(in)    :: n_levels_to_stabilise
  real(kind=r_def),                      intent(in)    :: max_stabilisation
  integer(kind=i_def), dimension(ndf1),  intent(in)    :: map1

  ! Internal variables
  integer(kind=i_def) :: df, k
  real(kind=r_def)    :: scalar, alpha(0 : nlayers - 1), level_as_float

  ! Set alpha, vector of stabilisation strength in the vertical
  level_as_float = real( n_levels_to_stabilise, r_def )
  scalar = max_stabilisation / level_as_float
  do k = 0, n_levels_to_stabilise - 1
    level_as_float = real( n_levels_to_stabilise - k, r_def )
    alpha(k) = scalar * level_as_float
  end do
  do k = n_levels_to_stabilise, nlayers - 1
    alpha(k) = 0.0_r_def
  end do

  ! Perform stabilisation
  do k = 0, nlayers - 1
    do df = 1, ndf1
      u_stabilised(map1(df) + k) =                       &
        alpha(k) * u_initial(map1(df) + k)               &
        + (1.0_r_def - alpha(k)) * u_final(map1(df) + k)
    end do
  end do

end subroutine stabilise_bl_u_code

end module  stabilise_bl_u_kernel_mod

