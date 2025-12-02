! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Global bypass module for switches concerned with microphysics

module mphys_bypass_mod

use um_types, only: real_umphys

implicit none
save

!  Description:
!   Module containing logical switches used by the microphysics code
!   that are not part of the run_precip or run_cloud namelists

!  Method:
!   Variables declared here are initialised in atm_step and then
!   used in the microphysics scheme.

!   This module is only for switches not in the run_cloud or
!   run_precip namelists. New microphysics or cloud scheme
!   logicals should be put in the run_precip or run_cloud namelist
!   and not here.

!   In theory, this module can be made redundant if all the switches
!   below are put into other modules.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

!-----------------------------------------------------
! 1. Microphysics Switches not part of run_precip
!-----------------------------------------------------

logical :: l_crystals     = .false.
! Controls whether a crystal mode is turned on or not.

logical :: l_ref_diag     = .false.
! Controls whether the subroutine mphys_reflec is called
! or not (called if .true.). If no reflectivity diagnostics
! are required, this should save time by skipping this bit
! of code.

!-----------------------------------------------------
! 1. Top of model;
! required in microphysics
!-----------------------------------------------------

real(kind=real_umphys) :: mphys_mod_top

!-----------------------------------------------------
! 2. Dimensions of qcf2; required for electric scheme
!-----------------------------------------------------

integer :: qcf2_idims_start
integer :: qcf2_idims_end
integer :: qcf2_jdims_start
integer :: qcf2_jdims_end

integer, parameter :: qcf2_kdims_start = 1

integer :: qcf2_kdims_end

end module mphys_bypass_mod
