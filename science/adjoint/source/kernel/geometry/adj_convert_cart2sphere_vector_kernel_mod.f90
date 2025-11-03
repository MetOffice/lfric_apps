!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Adjoint of the sci_psykal_builtin_light_mod routines
module adj_convert_cart2sphere_vector_kernel_mod

use argument_mod,            only : arg_type, GH_FIELD,     &
                                    GH_READWRITE, GH_READ,  &
                                    GH_REAL, DOF, ANY_SPACE_1
use constants_mod,           only : r_def
use kernel_mod,              only : kernel_type

implicit none

!---------------------------------------------------------------------------
! Public types
!---------------------------------------------------------------------------

type, public, extends(kernel_type) :: adj_convert_cart2sphere_vector_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                             &
        arg_type(GH_FIELD*3, GH_REAL, GH_READWRITE, ANY_SPACE_1), &
        arg_type(GH_FIELD*3, GH_REAL, GH_READ,      ANY_SPACE_1)  &
        /)
  integer :: operates_on = DOF
contains
  procedure, nopass :: adj_convert_cart2sphere_vector_code
end type adj_convert_cart2sphere_vector_kernel_type


!---------------------------------------------------------------------------
! Contained functions/subroutines
!---------------------------------------------------------------------------
public adj_convert_cart2sphere_vector_code
contains

!---------------------------------------------------------------------
subroutine adj_convert_cart2sphere_vector_code(field_vec_1, field_vec_2, &
                                               field_vec_3, coord_vec_1, &
                                               coord_vec_2, coord_vec_3  )
  use adj_coord_transform_mod, only : adj_cart2sphere_scalar
  implicit none

  real(kind=r_def), intent(inout) :: field_vec_1, field_vec_2, &
                                     field_vec_3
  real(kind=r_def), intent(in)    :: coord_vec_1, coord_vec_2, &
                                     coord_vec_3

  call adj_cart2sphere_scalar(coord_vec_1, coord_vec_2, &
                              coord_vec_3, field_vec_1, &
                              field_vec_2, field_vec_3 )

end subroutine adj_convert_cart2sphere_vector_code

end module adj_convert_cart2sphere_vector_kernel_mod
