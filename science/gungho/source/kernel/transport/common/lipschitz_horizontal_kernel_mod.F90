!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the Lipschitz number for the horizontal departure distances.
!> @details The Lipschitz number gives the stability condition for a semi-
!!          Lagrangian transport scheme. In 1D, this is simply the dimensionless
!!          departure distance, minus the value immediately upwind of it. If
!!          this number reaches 1, it implies that trajectories of particles
!!          have crossed, and the scheme will become unstable.
!!
!!          This kernel computes this Lipschitz number for the horizontal
!!          departure distances. It loops over columns, but only calculates the
!!          upwind values in that column.
!!
!!          NB: as this kernel computes only the upwind values in the column,
!!          to compute values at DoFs on the edge of a processor's local domain,
!!          we would need this kernel to loop up to depth one in the halos.
!!
!!          This could be done by setting the metadata for the Lipschitz number
!!          to be GH_INC, but this would force a halo exchange on the departure
!!          distance field. We don't want to do this because (a) it adds cost,
!!          and (b) the departure distance values at panel boundaries are
!!          deliberately different if they have been computed from det(J) on the
!!          extended mesh.
!!
!!          Since the Lipschitz number is only a diagnostic, it is better to
!!          accept that values at local domain edges may not be written on some
!!          processors. This important value from this field is the maximum,
!!          which will not be affected by this (since the value will be computed
!!          on at least one processor).

module lipschitz_horizontal_kernel_mod

  use argument_mod,      only: arg_type, GH_READ, &
                               GH_FIELD, GH_REAL, &
                               GH_WRITE, CELL_COLUMN
  use constants_mod,     only: r_tran, i_def
  use fs_continuity_mod, only: W2H
  use kernel_mod,        only: kernel_type
  use reference_element_mod, only: W, E, N, S

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: lipschitz_horizontal_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2H),                           &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2H)                            &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: lipschitz_horizontal_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: lipschitz_horizontal_code

contains

!> @brief Computes the horizontal Lipschitz number
!> @param[in]     nlayers           Number of layers in the mesh
!> @param[in,out] lipschitz_number  Lipschitz number of horizontal dep pts
!> @param[in]     dep_pts_xy        Horizontal departure distances
!> @param[in]     ndf_w2h           Num of DoFs per cell for W2H
!> @param[in]     undf_w2h          Num of DoFs per partition for W2H
!> @param[in]     map_w2h           DoFmap for W2H
subroutine lipschitz_horizontal_code( nlayers,                    &
                                      lipschitz_number,           &
                                      dep_pts_xy,                 &
                                      ndf_w2h, undf_w2h, map_w2h  &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h

  real(kind=r_tran),   dimension(undf_w2h), intent(inout) :: lipschitz_number
  real(kind=r_tran),   dimension(undf_w2h), intent(in)    :: dep_pts_xy
  integer(kind=i_def), dimension(ndf_w2h),  intent(in)    :: map_w2h

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_tran)   :: courant, upwind_courant

  ! Lipschitz for W ------------------------------------------------------------
  do k = 0, nlayers - 1
    courant = dep_pts_xy(map_w2h(W)+k)
    if (courant < 0.0_r_tran) then
      upwind_courant = dep_pts_xy(map_w2h(E)+k)
      lipschitz_number(map_w2h(W)+k) = upwind_courant - courant
    end if
  end do

  ! Lipschitz for E ------------------------------------------------------------
  do k = 0, nlayers - 1
    courant = dep_pts_xy(map_w2h(E)+k)
    if (courant > 0.0_r_tran) then
      upwind_courant = dep_pts_xy(map_w2h(W)+k)
      lipschitz_number(map_w2h(E)+k) = courant - upwind_courant
    end if
  end do

  ! Lipschitz for S ------------------------------------------------------------
  do k = 0, nlayers - 1
    courant = dep_pts_xy(map_w2h(S)+k)
    if (courant > 0.0_r_tran) then
      upwind_courant = dep_pts_xy(map_w2h(N)+k)
      lipschitz_number(map_w2h(S)+k) = courant - upwind_courant
    end if
  end do

  ! Lipschitz for N ------------------------------------------------------------
  do k = 0, nlayers - 1
    courant = dep_pts_xy(map_w2h(N)+k)
    if (courant < 0.0_r_tran) then
      upwind_courant = dep_pts_xy(map_w2h(S)+k)
      lipschitz_number(map_w2h(N)+k) = upwind_courant - courant
    end if
  end do

end subroutine lipschitz_horizontal_code

end module lipschitz_horizontal_kernel_mod
