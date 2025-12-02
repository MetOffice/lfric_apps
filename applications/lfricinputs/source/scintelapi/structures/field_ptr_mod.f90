! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module field_ptr_mod
!
! This module defines and provides access to the field pointer object derived
! type. This object simply contains a pointer to a LFRic field type. Its
! purpose is to allow the ability to get around the (current) restriction that
! one cannot create allocatable pointer arrays in Fortran.
!

use field_mod, only: field_type

implicit none

type field_ptr_object
  type(field_type), pointer :: field_ptr => null()
end type field_ptr_object

end module field_ptr_mod
