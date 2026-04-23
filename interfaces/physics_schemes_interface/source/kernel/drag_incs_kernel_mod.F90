!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Filters the drag increments in W2 to prevent drag from accelerating
!>        the wind
module drag_incs_kernel_mod

  use argument_mod,          only: arg_type, CELL_COLUMN,                      &
                                   GH_FIELD, GH_WRITE, GH_READ,                &
                                   GH_REAL, GH_INTEGER,                        &
                                   ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,         only: i_def, r_def
  use fs_continuity_mod,     only: W2
  use kernel_mod,            only: kernel_type
  use reference_element_mod, only: N

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  ! The type declaration for the kernel. Contains the metadata needed by the
  ! Psy layer.
  !
  type, public, extends(kernel_type) :: drag_incs_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
        arg_type(GH_FIELD, GH_REAL,    GH_WRITE, W2),                          &
        arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                          &
        arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                          &
        arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   &
        arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)    &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: drag_incs_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: drag_incs_code

contains

!> @param[in]  nlayers Number of layers
!> @param[out] du_out            Filtered output increment
!> @param[in]  du_in             Unfiltered input increment
!> @param[in]  u_in              Wind to be incremented
!> @param[in]  face_selector_ew  Face selector for east-west faces
!> @param[in]  face_selector_ns  Face selector for north-south faces
!> @param[in]  ndf               Number of degrees of freedom per cell
!> @param[in]  undf              Total number of degrees of freedom
!> @param[in]  map               Dofmap for the cell at the base of the column
!> @param[in]  ndf_2d            Num DoFs per cell for 2D field
!> @param[in]  undf_2d           Total num DoFs for 2D field in this partition
!> @param[in]  map_2d            Dofmap for the 2D field
subroutine drag_incs_code(nlayers, du_out, du_in, u_in,       &
                          face_selector_ew, face_selector_ns, &
                          ndf, undf, map, ndf_2d, undf_2d, map_2d)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, undf, ndf_2d, undf_2d
  integer(kind=i_def), dimension(ndf), intent(in) :: map
  integer(kind=i_def), dimension(ndf_2d), intent(in) :: map_2d
  real(kind=r_def), dimension(undf), intent(inout) :: du_out
  real(kind=r_def), dimension(undf), intent(in) :: du_in, u_in
  integer(kind=i_def), dimension(undf_2d), intent(in) :: face_selector_ew
  integer(kind=i_def), dimension(undf_2d), intent(in) :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: j, df, k

  ! Only loop over horizontal DoFs
  do j = 1, face_selector_ew(map_2d(1)) + face_selector_ns(map_2d(1))
    df = j
    if (j == 3 .and. face_selector_ns(map_2d(1)) == 2                          &
               .and. face_selector_ew(map_2d(1)) == 1) df = N
    do k = 0, nlayers-1
      ! Ensure drag is always acting in opposite sense to the wind,
      ! i.e. u and du aren't the same sign
      if (sign(1.0_r_def,u_in(map(df)+k)) == sign(1.0_r_def,du_in(map(df)+k))) then
        ! Set output increment to 0
        du_out(map(df)+k) = 0.0_r_def
      ! Ensure drag does not decelerate the wind beyond 0,
      ! i.e. u+du is the same sign as u
      else if ( sign(1.0_r_def,u_in(map(df)+k)+du_in(map(df)+k)) &
            /= sign(1.0_r_def,u_in(map(df)+k)) ) then
        ! Set output increment to negative of input wind
        du_out(map(df)+k) = -u_in(map(df)+k)
      else
        ! Output increment is input increment
        du_out(map(df)+k) = du_in(map(df)+k)
      end if
    end do
  end do

end subroutine drag_incs_code

end module drag_incs_kernel_mod
