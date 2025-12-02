! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************
!
!  Data module for switches for UM inputs,carbon options namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Carbon

module carbon_options_mod

implicit none

!===========================================================================
! items not set in namelist
!===========================================================================
! Carbon cycle
! Interactive 3D CO2 field for use with carbon cycle model
! This is set by the check_carbon_options routine based on the value of
! i_co2_opt
logical :: l_co2_interactive  = .false.



end module carbon_options_mod
