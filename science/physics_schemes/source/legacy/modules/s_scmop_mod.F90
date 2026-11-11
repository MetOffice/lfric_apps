!
! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
! Single column module for SCM diagnostic output
!

module s_scmop_mod

use scmoptype_defn, only: i64, SCMop_type, inot_written, maxnstreams, incdf,   &
                          maxndomprof

implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   Parameters/functions for the SCM diagnostics output system
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: scm
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!-----------------------------------------------------------------------------

! Integers to represent the different time profiles. All must
! be non-negative and less than "only_radsteps".

integer, parameter ::                                                          &
  t_inst        = 1   &! Give the instantaneous value
, t_avg         = 2   &! Construct the average value
, t_max         = 3   &! " maximum value
, t_min         = 4   &! " minimum value
, t_acc         = 5   &! " accumulated value
, t_div         = 7   &! " average value divided
                       !   by another diagnostic
, t_mult        = 8   &! " average value multiplied
                       !   by another diagnostic
, t_acc_div     = 9   &! " accumulated value divided
                       !   by another diagnostic
, t_acc_mult    = 10  &! " accumulated value multiplied
                       !   by another diagnostic
, t_const       = 11  &! The value is constant.
, only_radsteps = 100  ! When added to one of the above parameters,
                       ! flags that the diagnostic is only available
                       ! on radiation timesteps


! Integers to represent the different domain profiles
integer, parameter ::                                                          &
  d_sl      = 1                                                                &
, d_soilt   = 2                                                                &
, d_bl      = 3                                                                &
, d_wet     = 4                                                                &
, d_all     = 5                                                                &
, d_soilm   = 6                                                                &
, d_tile    = 7                                                                &
, d_vis     = 9                                                                &
, d_point   = 13                                                               &
, d_allxtra = 14                                                               &
, d_land    = 15                                                               &
, d_cloud   = 16                                                               &
, d_mlsnow  = 17

integer, parameter ::                                                          &
  default_streams = 63

! Integers to represent the different diagnostics packages
integer, parameter ::                                                          &
  SCMDiag_gen   = 1   & ! General diagnostics
, SCMDiag_rad   = 2   & ! Radiation
, SCMDiag_bl    = 3   & ! Boundary layer
, SCMDiag_surf  = 4   & ! Surface
, SCMDiag_land  = 5   & ! Land points only
, SCMDiag_sea   = 6   & ! Sea points only
, SCMDiag_lsp   = 7   & ! Large scale precip
, SCMDiag_conv  = 8   & ! Convection
, SCMDiag_lscld = 9   & ! Large scale cloud
, SCMDiag_pc2   = 10  & ! PC2
, SCMDiag_forc  = 11  & ! Forcing
, SCMDiag_incs  = 12  & ! Increments
, SCMDiag_gwd   = 13  & ! Gravity Wave Drag
, SCMDiag_MLSnow= 14  & ! Multilayer Snow
, SCMDiag_convss= 15    ! Extra convection diagnostics for every conv sub-step

! These need to be declared for the statement functions below
integer(i64) :: streamlist

contains

function stream(strm)

implicit none

! Statement function to encode a stream number into an integer
! The default streams for diagnostics to go to will be 1,2,3,4,5 and 6.
! The following should thus be equal to:
!
! Stream(1) [2^0=1] + Stream(2) [2^1=2]  + Stream(3) [2^2=4]
! Stream(4) [2^3=8] + Stream(5) [2^4=16] + Stream(6) [2^5=32]
! Total = 63
!
! where Stream() is the statement function defined above.
! default is 63 (all)

integer :: stream, strm

stream = 2**(strm-1)

end function stream

function stepnmbr(SCMop)

implicit none

! Statement function to translate daycount and stepcount into
! an integer representing the number of timesteps since the
! start of the run

integer :: stepnmbr
type(SCMop_type) :: SCMop

stepnmbr = (SCMop%daycount-1)*SCMop%full_daysteps                              &
         +  SCMop%stepcount

end function stepnmbr

function DoNotWrite(streamlist)

implicit none


! Statement function to modify an encoded stream list
! to flag that the streams should not be written to file

integer      :: DoNotWrite
integer(i64) :: streamlist
DoNotWrite = streamlist+2**(inot_written-1)

end function DoNotWrite


function StreamIsOn(streamlist,strm)

implicit none

! Statement function to test whether a particular stream is
! switched on in an encoded stream list

logical :: StreamIsOn
integer :: strm
integer(i64) :: streamlist

StreamIsOn =                                                                   &
   (streamlist-int(streamlist/Stream(strm+1))*Stream(strm+1))/                 &
   Stream(strm) >= 1

end function StreamIsOn

function NotWritten(streamlist)

implicit none

! Statement function to test whether the "do not write" flag
! has been switched on in an encoded stream list

logical :: NotWritten
integer(i64) :: streamlist

NotWritten = StreamIsOn(streamlist,inot_written)

end function NotWritten

function AnyStreamOn(streamlist)

implicit none

! Statement function to test whether any stream is switched
! on in an encoded stream list
logical :: AnyStreamOn
integer(i64) :: streamlist
AnyStreamOn =                                                                  &
     (.not. NotWritten(streamlist) .and. streamlist /= 0_i64)                  &
     .or.                                                                      &
     (NotWritten(streamlist) .and. streamlist /= DoNotWrite(0_i64))

end function AnyStreamOn

end module s_scmop_mod

