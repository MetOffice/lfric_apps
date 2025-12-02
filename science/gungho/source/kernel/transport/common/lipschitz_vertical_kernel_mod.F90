!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the Lipschitz number for the vertical departure distances.
!> @details The Lipschitz number gives the stability condition for a semi-
!!          Lagrangian transport scheme. In 1D, this is simply the dimensionless
!!          departure distance, minus the value immediately upwind of it. If
!!          this number reaches 1, it implies that trajectories of particles
!!          have crossed, and the scheme will become unstable.
module lipschitz_vertical_kernel_mod

  use argument_mod,      only: arg_type, GH_READ, &
                               GH_FIELD, GH_REAL, &
                               GH_WRITE, CELL_COLUMN
  use constants_mod,     only: r_tran, i_def
  use fs_continuity_mod, only: W2V
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: lipschitz_vertical_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2V),                           &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2V)                            &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: lipschitz_vertical_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: lipschitz_vertical_code

contains

!> @brief Adjusts vertical departure points to avoid them breaking the
!!        Lipschitz condition, and trajectories crossing
!> @param[in]     nlayers             Number of layers in the mesh
!> @param[in,out] lipschitz_vertical  Lipschitz number of vertical dep pts
!> @param[in]     dep_pts_z           Vertical departure points
!> @param[in]     ndf_w2v             Num of DoFs per cell for W2V
!> @param[in]     undf_w2v            Num of DoFs per partition for W2V
!> @param[in]     map_w2v             DoFmap for W2V
subroutine lipschitz_vertical_code( nlayers,                    &
                                    lipschitz_vertical,         &
                                    dep_pts_z,                  &
                                    ndf_w2v, undf_w2v, map_w2v  &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_w2v, undf_w2v

  real(kind=r_tran),   dimension(undf_w2v), intent(inout) :: lipschitz_vertical
  real(kind=r_tran),   dimension(undf_w2v), intent(in)    :: dep_pts_z
  integer(kind=i_def), dimension(ndf_w2v),  intent(in)    :: map_w2v

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_tran)   :: courant, upwind_courant

  lipschitz_vertical(map_w2v(1)) = 0.0_r_tran

  do k = 1, nlayers - 1
    courant = dep_pts_z(map_w2v(1)+k)
    ! Work out upwind value
    if (courant > 0.0_r_tran) then
      upwind_courant = dep_pts_z(map_w2v(1)+k-1)
      lipschitz_vertical(map_w2v(1)+k) = courant - upwind_courant
    else if (courant < 0.0_r_tran) then
      upwind_courant = dep_pts_z(map_w2v(1)+k+1)
      lipschitz_vertical(map_w2v(1)+k) = upwind_courant - courant
    else
      lipschitz_vertical(map_w2v(1)+k) = 0.0_r_tran
    end if

  end do

  lipschitz_vertical(map_w2v(1)+nlayers) = 0.0_r_tran

end subroutine lipschitz_vertical_code

end module lipschitz_vertical_kernel_mod
