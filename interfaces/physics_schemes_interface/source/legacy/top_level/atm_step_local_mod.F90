! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Declares variables and allocatable arrays which are local to atm_step
!  (and below)
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level

module atm_step_local

implicit none


integer :: rhc_row_length
integer :: rhc_rows

integer :: dim_cs1=1 ! soil C dimension: 1 for single, 4 for 4-pool

integer :: co2_dim_len ! For dimension 3-D CO2 field to be passed
integer :: co2_dim_row !     to NI_bl_ctl

!  stashworki = stashwork for section i
real, allocatable:: stashwork34(:),stashwork38(:), stashwork50(:)

end module atm_step_local

