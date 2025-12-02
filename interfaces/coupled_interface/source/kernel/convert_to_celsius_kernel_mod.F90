!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Kernel which converts from Kelvin to Celsius
!> @brief but only for temperatures > 1K

module convert_to_celsius_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,            &
                                    GH_REAL, GH_INTEGER, &
                                    GH_FIELD, GH_SCALAR, &
                                    GH_WRITE, GH_READ,   &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def

use sci_constants_mod,       only : zero_C_in_K

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: convert_to_celsius_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)          &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: convert_to_celsius_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: convert_to_celsius_code
contains

!> @brief Kernel which converts from Kelvin to Celsius
!> @brief but only for temperatures > 1K
!! @param[in] nlayers Number of layers
!! @param[in,out] output_field Output field to write to
!! @param[in]     input_field  Input field to read from
!! @param[in] ndata The size of the multidata dimension
!! @param[in] ndf Number of degrees of freedom per cell for the function space
!! @param[in] undf Number of unique degrees of freedom  for the function space
!! @param[in] map Dofmap for the cell at the base of the column for the function space

subroutine convert_to_celsius_code(nlayers,                     &
                      output_field, input_field,                &
                      ndata,                                    &
                      ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndata
  integer(kind=i_def), intent(in) :: ndf
  integer(kind=i_def), intent(in) :: undf
  integer(kind=i_def), dimension(ndf),  intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout) :: output_field
  real(kind=r_def), dimension(undf), intent(in)    :: input_field

  ! Internal variables
  integer(kind=i_def) :: i, df

  do i = 0, ndata - 1
    do df = 1,ndf
      if (input_field(map(df) + i) > 1.0_r_def) then
        output_field(map(df) + i) = input_field(map(df) + i) - zero_C_in_K
      else
        output_field(map(df) + i) = 0.0_r_def
      end if
    end do
  end do

end subroutine convert_to_celsius_code

end module convert_to_celsius_kernel_mod
