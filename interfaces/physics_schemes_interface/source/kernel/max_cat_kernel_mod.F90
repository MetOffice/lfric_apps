!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Running maximum of the clear air turbulence predictor over
!>        pressure levels, and the pressure at which it occurs.

module max_cat_kernel_mod

  use argument_mod,      only: arg_type,                   &
                               GH_FIELD, GH_SCALAR,        &
                               GH_REAL, GH_INTEGER,        &
                               GH_READ, GH_READWRITE,      &
                               CELL_COLUMN,                &
                               ANY_DISCONTINUOUS_SPACE_1,  &
                               ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: max_cat_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                             &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,              &
                  ANY_DISCONTINUOUS_SPACE_1),                      & ! max_cat
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,              &
                  ANY_DISCONTINUOUS_SPACE_1),                      & ! max_press
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,                   &
                  ANY_DISCONTINUOUS_SPACE_2),                      & ! cat
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                  & ! kp
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)                   & ! press
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: max_cat_code
  end type max_cat_kernel_type

  public :: max_cat_code

contains

  !> @brief   Update the running maximum of the predictor with one level.
  !> @details Where the predictor on the given level exceeds the running
  !>          maximum, the maximum and its pressure are replaced. The test
  !>          is strict, so the first level visited wins any tie, and a
  !>          missing predictor never replaces the running value. This is
  !>          the level loop of the UM PWS driver (pws_diags_driver_mod)
  !>          for the maximum clear air turbulence diagnostics, one level
  !>          per call; the caller initialises both outputs to missing data
  !>          before the first level.
  !> @param[in]     nlayers    Number of layers
  !> @param[in,out] max_cat    Maximum predictor over the levels so far
  !> @param[in,out] max_press  Pressure of that maximum
  !> @param[in]     cat        Clear air turbulence predictor on pressure
  !>                           levels
  !> @param[in]     kp         Slot of the pressure-level field to visit
  !> @param[in]     press_lev  Pressure of that slot (Pa)
  !> @param[in]     ndf_2d     Number of degrees of freedom per cell for the
  !>                           single-level fields
  !> @param[in]     undf_2d    Number of unique degrees of freedom for the
  !>                           single-level fields
  !> @param[in]     map_2d     Dofmap for the cell for the single-level
  !>                           fields
  !> @param[in]     ndf_pl     Number of degrees of freedom per cell for the
  !>                           pressure-level field
  !> @param[in]     undf_pl    Number of unique degrees of freedom for the
  !>                           pressure-level field
  !> @param[in]     map_pl     Dofmap for the cell for the pressure-level
  !>                           field
  subroutine max_cat_code(nlayers,                       &
                          max_cat, max_press,            &
                          cat,                           &
                          kp, press_lev,                 &
                          ndf_2d, undf_2d, map_2d,       &
                          ndf_pl, undf_pl, map_pl)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d
    integer(kind=i_def), intent(in) :: ndf_pl, undf_pl
    integer(kind=i_def), intent(in), dimension(ndf_pl) :: map_pl

    real(kind=r_def), intent(inout), dimension(undf_2d) :: max_cat
    real(kind=r_def), intent(inout), dimension(undf_2d) :: max_press
    real(kind=r_def), intent(in),    dimension(undf_pl) :: cat

    integer(kind=i_def), intent(in) :: kp
    real(kind=r_def),    intent(in) :: press_lev

    ! Internal variables
    integer(kind=i_def) :: slot

    slot = map_pl(1) + kp - 1_i_def

    if ( cat(slot) > max_cat(map_2d(1)) ) then
      max_cat(map_2d(1))   = cat(slot)
      max_press(map_2d(1)) = press_lev
    end if

  end subroutine max_cat_code

end module max_cat_kernel_mod
