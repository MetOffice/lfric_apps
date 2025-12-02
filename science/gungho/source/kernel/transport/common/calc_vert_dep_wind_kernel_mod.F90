!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the vertical departure wind.
!> @details Computes the wind used for the vertical departure points by dividing
!!          the advecting wind by the Det(J) in the upwind cell.

module calc_vert_dep_wind_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 CELL_COLUMN, GH_WRITE, &
                                 GH_READ
  use constants_mod,      only : i_def, r_tran
  use fs_continuity_mod,  only : W3, W2, W2v
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: calc_vert_dep_wind_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                  &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! dep_wind
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),  & ! adv_wind
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3)   & ! detj
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_vert_dep_wind_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: calc_vert_dep_wind_code

contains

  !> @brief Compute the horizontal departure wind.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] dep_wind       Departure wind on W2v DOFs
  !> @param[in]     adv_wind       Advecting wind on W2 DOFs
  !> @param[in]     detj           Volume factor at W3
  !> @param[in]     ndf_w2v        Number of degrees of freedom for W2v per cell
  !> @param[in]     undf_w2v       Number of unique degrees of freedom for W2v
  !> @param[in]     map_w2v        Map for W2v
  !> @param[in]     ndf_w2         Number of degrees of freedom for W2 per cell
  !> @param[in]     undf_w2        Number of unique degrees of freedom for W2
  !> @param[in]     map_w2         Map for W2
  !> @param[in]     ndf_w3         Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3        Number of unique degrees of freedom for W3
  !> @param[in]     map_w3         Map for W3

  subroutine calc_vert_dep_wind_code( nlayers,      &
                                      dep_wind,     &
                                      adv_wind,     &
                                      detj,         &
                                      ndf_w2v,      &
                                      undf_w2v,     &
                                      map_w2v,      &
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
    integer(kind=i_def), intent(in) :: undf_w2v
    integer(kind=i_def), intent(in) :: ndf_w2v
    integer(kind=i_def), intent(in) :: undf_w2
    integer(kind=i_def), intent(in) :: ndf_w2

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
    integer(kind=i_def), dimension(ndf_w2v), intent(in) :: map_w2v

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_w2v), intent(inout) :: dep_wind
    real(kind=r_tran), dimension(undf_w2),  intent(in)    :: adv_wind
    real(kind=r_tran), dimension(undf_w3),  intent(in)    :: detj

    ! Variables
    integer(kind=i_def) :: k

    ! W2 has vertical DOFs 5 and 6
    ! W2v has vertical DOFs 1 and 2

    ! Set top and bottom departure winds to zero
    dep_wind(map_w2v(1)) = 0.0_r_tran
    dep_wind(map_w2v(2)+nlayers-1) = 0.0_r_tran

    do k=1, nlayers-1

      ! For the bottom W2v dof of each cell divide by the upwind Det(J)

      dep_wind(map_w2v(1)+k) = (0.5_r_tran + sign(0.5_r_tran, adv_wind(map_w2(5)+k)) )   &
                               * adv_wind(map_w2(5)+k) / detj(map_w3(1) + k - 1)         &
                               + (0.5_r_tran - sign(0.5_r_tran, adv_wind(map_w2(5)+k)) ) &
                               * adv_wind(map_w2(5)+k) / detj(map_w3(1) + k)

    end do

  end subroutine calc_vert_dep_wind_code

end module calc_vert_dep_wind_kernel_mod












