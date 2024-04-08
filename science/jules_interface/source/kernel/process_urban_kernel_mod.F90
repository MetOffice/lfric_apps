!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Process Jules soil ancillaries
module process_urban_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_REAL,         &
                           GH_READ, GH_WRITE,         &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           CELL_COLUMN
  use constants_mod, only: r_def, i_def, r_um
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: process_urban_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),  &! wrr
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),  &! hwr
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),  &! hgt
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &! ztm
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &! disp
       /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: process_urban_code
  end type process_urban_kernel_type

  public :: process_urban_code

contains

  !> @param[in]     nlayers                The number of layers
  !> @param[in]     urbwrr                 Urban repeating width ratio
  !> @param[in]     urbhwr                 Urban height-to-width ratio
  !> @param[in]     urbhgt                 Urban building height
  !> @param[in,out] urbztm                 Urban effective roughness length
  !> @param[in,out] urbdisp                Urban displacement height
  !> @param[in]     ndf_2d                 Number of DOFs per cell for 2d fields
  !> @param[in]     undf_2d                Number of total DOFs for 2d fields
  !> @param[in]     map_2d                 Dofmap for cell for surface 2d fields
  subroutine process_urban_code(nlayers,                       &
                                ! IN
                                urbwrr, urbhwr, urbhgt,        &
                                ! OUT
                                urbztm, urbdisp,               &
                                ndf_2d, undf_2d, map_2d)

  use calc_urban_aero_fields_mod, only: calc_urban_aero_fields

  implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

    real(kind=r_def), intent(in)    :: urbwrr(undf_2d)
    real(kind=r_def), intent(in)    :: urbhwr(undf_2d)
    real(kind=r_def), intent(in)    :: urbhgt(undf_2d)
    real(kind=r_def), intent(inout) :: urbztm(undf_2d)
    real(kind=r_def), intent(inout) :: urbdisp(undf_2d)

    real(r_um), dimension(1) :: ztm, disp

    call calc_urban_aero_fields(1,                                             &
                                real(urbwrr(map_2d(1):map_2d(1)), r_um),       &
                                real(urbhwr(map_2d(1):map_2d(1)), r_um),       &
                                real(urbhgt(map_2d(1):map_2d(1)), r_um),       &
                                ztm, disp )

    urbztm(map_2d(1))  = real(ztm(1),  r_def)
    urbdisp(map_2d(1)) = real(disp(1), r_def)

  end subroutine process_urban_code

end module process_urban_kernel_mod
