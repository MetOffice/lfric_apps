module conv_type_defs

implicit none

! Description:
!  Defines convection types for tests in code
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection



integer, parameter ::                                                          &
  shallow_conv   = 1     & ! Indicates convection is shallow
 ,congestus_conv = 2     & ! Indicates convection is congestus
 ,deep_conv      = 3     & ! Indicates convection is deep
 ,mid_conv       = 4       ! Indicates convection is mid-level


end module conv_type_defs
