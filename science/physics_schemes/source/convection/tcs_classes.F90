! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing a collection of tcs classes and overloading
! allocation and deallocation subroutines
!
module tcs_classes

use tcs_class_interface,  only: allocate_interface_input,                      &
                                deallocate_interface_input,                    &
                                assign_interface_input
use tcs_class_similarity, only: allocate_similarity,                           &
                                deallocate_similarity,                         &
                                similarity
use tcs_class_scales,     only: allocate_scales, deallocate_scales
use tcs_class_cloud,      only: allocate_cloud_input,                          &
                                deallocate_cloud_input,                        &
                                assign_cloud_input,                            &
                                cloud_input

implicit none
!
! Description:
! This module defines the tcs warm rain "cloud_input" derived type and
! subroutines for allocating and deallocating an instance of this class.
!
!
! Method:
!  Creates new interfaces to overloads tcs_allocate and tcs_deallocate
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.

interface tcs_allocate
  module procedure allocate_scales
  module procedure allocate_similarity
  module procedure allocate_interface_input
  module procedure allocate_cloud_input
end interface

interface tcs_deallocate
  module procedure deallocate_scales
  module procedure deallocate_similarity
  module procedure deallocate_interface_input
  module procedure deallocate_cloud_input
end interface

end module tcs_classes
