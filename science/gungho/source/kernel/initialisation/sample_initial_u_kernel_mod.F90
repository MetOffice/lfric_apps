!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the rhs for the initialisation of the wind field.
!>
!> @details The kernel computes the rhs of the equation u = u0 where u0 is the
!>          analytically defined wind field. The analytic wind field is projected
!>          onto u using Galerkin projection.
!>
module sample_initial_u_kernel_mod

  use argument_mod,            only : arg_type,          &
                                      GH_FIELD, GH_REAL, &
                                      GH_SCALAR,         &
                                      GH_INC, GH_READ,   &
                                      ANY_SPACE_9, CELL_COLUMN
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W2
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: sample_initial_u_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                       &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),          &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9), &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),              &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)               &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: sample_initial_u_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: sample_initial_u_code

contains

!> @brief Compute the right hand side to initialise the wind field.
!! @param[in] nlayers Number of layers
!! @param[in,out] wind The velocity vector
!! @param[in] chi_1 X component of the coordinate field
!! @param[in] chi_2 Y component of the coordinate field
!! @param[in] chi_3 Z component of the coordinate field
!! @param[in] u0 Initial wind u component
!! @param[in] v0 Initial wind v component
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column
!! @param[in] ndf_chi Number of dofs per cell for the coordinate field
!! @param[in] undf_chi Total number of degrees of freedom
!! @param[in] map_chi Dofmap for the coordinate field
subroutine sample_initial_u_code(nlayers,                   &
                                 wind,                      &
                                 chi_1, chi_2, chi_3,       &
                                 u0, v0,                    &
                                 ndf, undf, map,            &
                                 ndf_chi, undf_chi, map_chi &
                                 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_chi
  integer(kind=i_def), intent(in) :: undf, undf_chi

  integer(kind=i_def), dimension(ndf),     intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(undf),     intent(inout) :: wind
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), intent(in) :: u0
  real(kind=r_def), intent(in) :: v0

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_def)    :: dx, dy, dz

  do k = 0, nlayers-1
    ! Assumes a linear DG coordinate field and constant horizontal wind: u=U0,v=V0.
    ! The vertical wind component is computed from the horizontal components.

    ! u dofs
    dz = 0.5_r_def*(chi_3(map_chi(5)+k) - chi_3(map_chi(1)+k) + chi_3(map_chi(7)+k) - chi_3(map_chi(3)+k))
    dy = 0.5_r_def*(chi_2(map_chi(3)+k) - chi_2(map_chi(1)+k) + chi_2(map_chi(7)+k) - chi_2(map_chi(5)+k))
    wind(map(1)+k) = U0*dy*dz
    dz = 0.5_r_def*(chi_3(map_chi(6)+k) - chi_3(map_chi(2)+k) + chi_3(map_chi(8)+k) - chi_3(map_chi(4)+k))
    dy = 0.5_r_def*(chi_2(map_chi(4)+k) - chi_2(map_chi(2)+k) + chi_2(map_chi(8)+k) - chi_2(map_chi(6)+k))
    wind(map(3)+k) = U0*dy*dz
    ! v dofs, The basis function has been neglected, but is valued -1 at v dof
    ! points (since it points in the negative chi_hat(2) direction)
    ! therefore a positive value of V0 corresponds to a negative value of
    ! wind(map(2)) & wind(map(4))
    dz = 0.5_r_def*(chi_3(map_chi(5)+k) - chi_3(map_chi(1)+k) + chi_3(map_chi(6)+k) - chi_3(map_chi(2)+k))
    dx = 0.5_r_def*(chi_1(map_chi(2)+k) - chi_1(map_chi(1)+k) + chi_1(map_chi(6)+k) - chi_1(map_chi(5)+k))
    wind(map(2)+k) = -V0*dx*dz
    dz = 0.5_r_def*(chi_3(map_chi(7)+k) - chi_3(map_chi(3)+k) + chi_3(map_chi(8)+k) - chi_3(map_chi(4)+k))
    dx = 0.5_r_def*(chi_1(map_chi(4)+k) - chi_1(map_chi(3)+k) + chi_1(map_chi(8)+k) - chi_1(map_chi(7)+k))
    wind(map(4)+k) = -V0*dx*dz

    ! The vertical component is left unchanged ,i.e. forced to be zero.
    ! Since the vertical component is the normal flux through vertical faces
    ! this means that the truly vertical velocity (in increasing z direction)
    ! will be: ~(u*dz/dx + v*dz/dy)
    wind(map(5)+k) = 0.0_r_def
    wind(map(6)+k) = 0.0_r_def
  end do

end subroutine sample_initial_u_code

end module sample_initial_u_kernel_mod
