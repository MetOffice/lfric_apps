!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Generate the land fraction
!>
module generate_land_frac_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_REAL,         &
                                      GH_READ, GH_WRITE,         &
                                      ANY_DISCONTINUOUS_SPACE_1, &
                                      ANY_DISCONTINUOUS_SPACE_2, &
                                      CELL_COLUMN
  use constants_mod,           only : i_def, r_def
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: generate_land_frac_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! land_fraction
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2)  & ! tile_fraction
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: generate_land_frac_code
  end type

  public :: generate_land_frac_code

contains

  !> @brief Generate the land fraction
  !> @details Get the 2D land fraction from the tile fractions
  !> @param[in]     nlayers                Number of layers
  !> @param[in,out] land_fraction          The fraction of land
  !> @param[in]     tile_fraction          The fraction of each tile
  subroutine generate_land_frac_code(      &
               nlayers, land_fraction,     &
               tile_fraction,              &
               ndf_2d,                     &
               undf_2d,                    &
               map_2d,                     &
               ndf_tile,                   &
               undf_tile,                  &
               map_tile                    &
                          )

    use jules_control_init_mod, only: n_land_tile

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d, ndf_tile, undf_tile
    integer(kind=i_def), dimension(ndf_2d),    intent(in) :: map_2d
    integer(kind=i_def), dimension(ndf_tile),  intent(in) :: map_tile

    real(kind=r_def), intent(inout) :: land_fraction(undf_2d)
    real(kind=r_def), intent(in)    :: tile_fraction(undf_tile)

    ! Local variables
    integer(kind=i_def) :: n

    land_fraction(map_2d(1)) = 0.0_r_def
    do n = 1, n_land_tile
      land_fraction(map_2d(1)) = land_fraction(map_2d(1)) + tile_fraction(map_tile(1)+n-1)
    end do

  end subroutine generate_land_frac_code

end module generate_land_frac_kernel_mod
