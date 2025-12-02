!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Selects the upwind Det(J) at vertical W2 locations.
!> @details Based on the sign of the wind, set the Det(J) at vertical W2
!>          locations to either be from the cell above (negative wind) or
!>          below (positive wind). The output is the upwind Det(J)
!>          at vertical W2 locations.
!>          This kernel is only designed for the lowest order finite-element spaces.
!>
module calc_upwind_detj_at_w2_kernel_mod

  use argument_mod,          only : arg_type, GH_FIELD,  &
                                    GH_REAL, GH_WRITE,   &
                                    GH_READ, CELL_COLUMN

  use constants_mod,         only : r_tran, r_def, i_def
  use fs_continuity_mod,     only : W2
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : T

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: calc_upwind_detj_at_w2_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/            &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_upwind_detj_at_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: calc_upwind_detj_at_w2_code

contains

!> @brief Selects the upwind Det(J) at vertical W2 locations based on the
!!        sign of the wind.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] detj_w2        Output field containing the upwind Det(J) values at W2 locations
!> @param[in]     detj_w2_above  Input field containing the Det(J) values at W2 from the cell above
!> @param[in]     detj_w2_below  Input field containing the Det(J) values at W2 from the cell below
!> @param[in]     wind           Input field containing the wind at W2 locations
!> @param[in]     ndf_w2         The number of degrees of freedom per cell for the output field
!> @param[in]     undf_w2        The number of unique degrees of freedom for the output field
!> @param[in]     map_w2         Array holding the dofmap for the cell at the base
!!                               of the column for the output field
!>
subroutine calc_upwind_detj_at_w2_code( nlayers,       &
                                        detj_w2,       &
                                        detj_w2_above, &
                                        detj_w2_below, &
                                        wind,          &
                                        ndf_w2,        &
                                        undf_w2,       &
                                        map_w2         &
                                      )

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in)    :: nlayers
  integer(kind=i_def),                    intent(in)    :: ndf_w2
  integer(kind=i_def),                    intent(in)    :: undf_w2
  real(kind=r_tran), dimension(undf_w2),  intent(inout) :: detj_w2
  real(kind=r_def),  dimension(undf_w2),  intent(in)    :: detj_w2_above
  real(kind=r_def),  dimension(undf_w2),  intent(in)    :: detj_w2_below
  real(kind=r_tran), dimension(undf_w2),  intent(in)    :: wind
  integer(kind=i_def), dimension(ndf_w2), intent(in)    :: map_w2

  ! Internal variables
  integer(kind=i_def) :: df, k

  ! Use existing Det(J) values in top and bottom levels and for horizontal dofs

  do k = 1, nlayers-2

    ! Set Det(J) based on the sign of the vertical wind
    ! Note: This assumes that the sign of wind is positive for upwards,
    !       and negative for downwards

    df = T

    detj_w2(map_w2(df)+k) =                                                      &
            0.5_r_tran*(1.0_r_tran - sign( 1.0_r_tran, wind(map_w2(df)+k) ) )*   &
            real(detj_w2_above(map_w2(df)+k), r_tran)                            &
            + 0.5_r_tran*(1.0_r_tran + sign( 1.0_r_tran, wind(map_w2(df)+k) ) )* &
            real(detj_w2_below(map_w2(df)+k), r_tran)

  end do

end subroutine calc_upwind_detj_at_w2_code

end module calc_upwind_detj_at_w2_kernel_mod
