! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE field_ptr_mod
!
! This module defines and provides access to the field pointer object derived
! type. This object simply contains a pointer to a LFRic field type. Its
! purpose is to allow the ability to get around the (current) restriction that
! one cannot create allocatable pointer arrays in Fortran.
!

USE field_mod, ONLY: field_type

IMPLICIT NONE

TYPE field_ptr_object
  TYPE(field_type), POINTER :: field_ptr => NULL()
END TYPE field_ptr_object

END MODULE field_ptr_mod
