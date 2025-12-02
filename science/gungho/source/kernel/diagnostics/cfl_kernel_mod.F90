!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Limits wind velocity.
!>
module cfl_kernel_mod

  use argument_mod,            only : arg_type,          &
                                      GH_FIELD, GH_REAL, &
                                      GH_READ, GH_INC,   &
                                      CELL_COLUMN
  use constants_mod,           only : i_def, r_def
  use fs_continuity_mod,       only : W2
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: cfl_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/              &
         arg_type(GH_FIELD*3, GH_REAL, GH_INC,  W2), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: cfl_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: cfl_code

contains

!> @brief Calculates components of the advective Courant number on
!>        W2 dofs
!! @param[in]     nlayers Number of layers
!! @param[in,out] cflx CFL calculated on 'x' dofs in 'x' direction
!! @param[in,out] cfly CFL calculated on 'y' dofs in 'y' direction
!! @param[in,out] cflz CFL calculated on 'z' dofs in 'z' direction
!! @param[in]     wind Wind
!! @param[in]     dJ_on_w2 detJ evaluated on w2 points
!! @param[in]     ndf_w2 Number of degrees of freedom per cell
!! @param[in]     undf_w2 Total number of degrees of freedom
!! @param[in]     map_w2 Dofmap for the cell at the base of the column
!! @param[in]     dt The model timestep length
subroutine cfl_code(nlayers, cflx, cfly, cflz, wind, dJ_on_w2, &
                    ndf_w2, undf_w2, map_w2, dt)

  use analytic_density_profiles_mod, only : analytic_density

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf_w2, undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(undf_w2), intent(inout) :: cflx, cfly, cflz
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_w2), intent(in)    :: dJ_on_w2
  real(kind=r_def),                     intent(in)    :: dt

  ! Internal variables
  integer(kind=i_def) :: k

  do k = 0, nlayers-1

      ! Compute the cfl for this dof
      cflx(map_w2(1)+k) = wind(map_w2(1)+k)*dt/dJ_on_w2(map_w2(1)+k)
      cfly(map_w2(2)+k) = wind(map_w2(2)+k)*dt/dJ_on_w2(map_w2(2)+k)
      cflx(map_w2(3)+k) = wind(map_w2(3)+k)*dt/dJ_on_w2(map_w2(3)+k)
      cfly(map_w2(4)+k) = wind(map_w2(4)+k)*dt/dJ_on_w2(map_w2(4)+k)
      cflz(map_w2(5)+k) = wind(map_w2(5)+k)*dt/dJ_on_w2(map_w2(5)+k)
      cflz(map_w2(6)+k) = wind(map_w2(6)+k)*dt/dJ_on_w2(map_w2(6)+k)

  end do

end subroutine cfl_code

end module cfl_kernel_mod
