! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Get STASH item number (in section 50) for a given diagnostics
!
! Method:
!   Get the name of the diagnostic as input argument and use a
!   case select statement to return the corresponding STASH
!   item number in section 50.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
module get_emdiag_stash_mod

implicit none

character(len=*), parameter, private :: ModuleName='GET_EMDIAG_STASH_MOD'

contains

function get_emdiag_stash (diag_name)

use ereport_mod,     only: ereport
use errormessagelength_mod, only: errormessagelength
use parkind1,        only: jpim,  jprb     ! DrHook
use yomhook,         only: lhook, dr_hook  ! DrHook

implicit none

! Function arguments
character (len=10), intent(in) :: diag_name ! name of the diagnostic
                                            ! (e.g. emission field)
! Function return value
integer :: get_emdiag_stash

! Local variables
integer                            :: ierr
character (len=errormessagelength) :: cmessage

integer (kind=jpim), parameter :: zhook_in  = 0  ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1  ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle   ! DrHook tracing

character(len=*), parameter :: RoutineName='GET_EMDIAG_STASH'

! End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0
get_emdiag_stash = -99

! Return STASH item in Section 50 for a given UKCA diagnostic.
! The item numbers below should agree with STASHmaster_A.
select case ( diag_name )
case ('NO        ')
  get_emdiag_stash = 156

case ('CH4       ')
  get_emdiag_stash = 157

case ('CO        ')
  get_emdiag_stash = 158

case ('HCHO      ')
  get_emdiag_stash = 159

case ('C2H6      ')
  get_emdiag_stash = 160

case ('C3H8      ')
  get_emdiag_stash = 161

case ('Me2CO     ')
  get_emdiag_stash = 162

case ('MeCHO     ')
  get_emdiag_stash = 163

case ('C5H8      ')
  get_emdiag_stash = 164

case ('C4H10     ')
  get_emdiag_stash = 165

case ('C2H4      ')
  get_emdiag_stash = 166

case ('C3H6      ')
  get_emdiag_stash = 167

case ('TOLUENE   ')
  get_emdiag_stash = 168

case ('oXYLENE   ')
  get_emdiag_stash = 169

case ('CH3OH     ')
  get_emdiag_stash = 170

case ('H2        ')
  get_emdiag_stash = 171

case ('NO_aircrft')
  get_emdiag_stash = 172

case ('Monoterp  ')
  get_emdiag_stash = 211

case ('MeOH      ')
  get_emdiag_stash = 212

case ('NH3       ')
  get_emdiag_stash = 213

case ('DMS       ')
  get_emdiag_stash = 214

case ('SO2_low   ')
  get_emdiag_stash = 215

case ('SO2_high  ')
  get_emdiag_stash = 216

case ('SO2_nat   ')
  get_emdiag_stash = 217

  ! New emission diagnostics for CRI species:
case ('EtOH      ')
  get_emdiag_stash = 304

case ('C2H2      ')
  get_emdiag_stash = 305

case ('TBUT2ENE  ')
  get_emdiag_stash = 306

case ('APINENE   ')
  get_emdiag_stash = 307

case ('BPINENE   ')
  get_emdiag_stash = 308

case ('BENZENE   ')
  get_emdiag_stash = 309

case ('HCOOH     ')
  get_emdiag_stash = 310

case ('MeCO2H    ')
  get_emdiag_stash = 311

case ('EtCHO     ')
  get_emdiag_stash = 312

case ('HOCH2CHO  ')
  get_emdiag_stash = 313

case ('MEK       ')
  get_emdiag_stash = 314

case DEFAULT
  ! Report error
  ierr = 1
end select

! Report error if diagnostics not found
if (ierr == 1) then
  cmessage = 'Unexpected UKCA emiss diagnostic: ' // trim (diag_name)
  call ereport ('GET_EMDIAG_STASH', ierr, cmessage)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

return

end function get_emdiag_stash

end module get_emdiag_stash_mod
