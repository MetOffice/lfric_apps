!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the horizontal departure wind.
!> @details Computes the wind used for the horizontal departure points by dividing
!!          the advecting wind by the Det(J) in the upwind cell.

module calc_departure_wind_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 CELL_COLUMN, GH_WRITE, &
                                 GH_READ, STENCIL, CROSS
  use constants_mod,      only : i_def, r_tran
  use fs_continuity_mod,  only : W3, W2, W2h
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: calc_departure_wind_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2h),               & ! dep_wind
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                & ! adv_wind
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(CROSS)) & ! detj
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_departure_wind_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: calc_departure_wind_code

contains

  !> @brief Compute the horizontal departure wind.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] dep_wind       Departure wind on W2h DOFs
  !> @param[in]     adv_wind       Advecting wind on W2 DOFs
  !> @param[in]     detj           Volume factor at W3
  !> @param[in]     stencil_size   Local length of Det(J) at W3 stencil
  !> @param[in]     stencil_map    Dofmap for the Det(J) at W3 stencil
  !> @param[in]     ndf_w2h        Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h       Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h        Map for W2h
  !> @param[in]     ndf_w2         Number of degrees of freedom for W2 per cell
  !> @param[in]     undf_w2        Number of unique degrees of freedom for W2
  !> @param[in]     map_w2         Map for W2
  !> @param[in]     ndf_w3         Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3        Number of unique degrees of freedom for W3
  !> @param[in]     map_w3         Map for W3

  subroutine calc_departure_wind_code( nlayers,      &
                                       dep_wind,     &
                                       adv_wind,     &
                                       detj,         &
                                       stencil_size, &
                                       stencil_map,  &
                                       ndf_w2h,      &
                                       undf_w2h,     &
                                       map_w2h,      &
                                       ndf_w2,       &
                                       undf_w2,      &
                                       map_w2,       &
                                       ndf_w3,       &
                                       undf_w3,      &
                                       map_w3 )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_w2
    integer(kind=i_def), intent(in) :: ndf_w2
    integer(kind=i_def), intent(in) :: stencil_size

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
    integer(kind=i_def), dimension(ndf_w3,stencil_size), intent(in) :: stencil_map

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_w2h), intent(inout) :: dep_wind
    real(kind=r_tran), dimension(undf_w2),  intent(in)    :: adv_wind
    real(kind=r_tran), dimension(undf_w3),  intent(in)    :: detj

    ! Variables
    integer(kind=i_def) :: k

    ! Cross stencil extent 1 has form
    !     | 5 |
    ! | 2 | 1 | 4 |
    !     | 3 |

    ! W2 has horizontal DOF map
    !       4
    !     1   3
    !       2

    if (stencil_size < 5) then
      ! Edge of LAM, set dep_wind to zero
      do k=0, nlayers-1
        dep_wind(map_w2h(1)+k) = 0.0_r_tran
        dep_wind(map_w2h(2)+k) = 0.0_r_tran
        dep_wind(map_w2h(3)+k) = 0.0_r_tran
        dep_wind(map_w2h(4)+k) = 0.0_r_tran
      end do
    else

      do k=0, nlayers-1

        ! For each horizontal W2 dof divide by the upwind Det(J)

        dep_wind(map_w2h(1)+k) = (0.5_r_tran + sign(0.5_r_tran, adv_wind(map_w2(1)+k)) ) &
                                 * adv_wind(map_w2(1)+k) / detj(stencil_map(1,2) + k)    &
                               + (0.5_r_tran - sign(0.5_r_tran, adv_wind(map_w2(1)+k)) ) &
                                 * adv_wind(map_w2(1)+k) / detj(stencil_map(1,1) + k)

        dep_wind(map_w2h(2)+k) = (0.5_r_tran - sign(0.5_r_tran, adv_wind(map_w2(2)+k)) ) &
                                 * adv_wind(map_w2(2)+k) / detj(stencil_map(1,3) + k)    &
                               + (0.5_r_tran + sign(0.5_r_tran, adv_wind(map_w2(2)+k)) ) &
                                 * adv_wind(map_w2(2)+k) / detj(stencil_map(1,1) + k)

        dep_wind(map_w2h(3)+k) = (0.5_r_tran + sign(0.5_r_tran, adv_wind(map_w2(3)+k)) ) &
                                 * adv_wind(map_w2(3)+k) / detj(stencil_map(1,1) + k)    &
                               + (0.5_r_tran - sign(0.5_r_tran, adv_wind(map_w2(3)+k)) ) &
                                 * adv_wind(map_w2(3)+k) / detj(stencil_map(1,4) + k)

        dep_wind(map_w2h(4)+k) = (0.5_r_tran - sign(0.5_r_tran, adv_wind(map_w2(4)+k)) ) &
                                 * adv_wind(map_w2(4)+k) / detj(stencil_map(1,1) + k)    &
                               + (0.5_r_tran + sign(0.5_r_tran, adv_wind(map_w2(4)+k)) ) &
                                 * adv_wind(map_w2(4)+k) / detj(stencil_map(1,5) + k)

      end do

    end if

  end subroutine calc_departure_wind_code

end module calc_departure_wind_kernel_mod

