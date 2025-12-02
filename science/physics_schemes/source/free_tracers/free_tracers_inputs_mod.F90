! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!  Data module for free tracers related switches and settings
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Free_Tracers

module free_tracers_inputs_mod

use max_dimensions,   only: ntype_max, snow_layers_max

use missing_data_mod, only: rmdi, imdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength

implicit none

!-----------------------------------------------------
! Items set in module
!-----------------------------------------------------

integer,parameter:: a_tracer_first = 1
! First atmospheric tracer (STASH No)

integer,parameter:: a_tracer_last = 150
! Last atmospheric tracer  (STASH No)

integer, parameter :: a_max_trvars = (a_tracer_last - a_tracer_first) + 1
! Max number of tracers allowed

integer, parameter :: max_wtrac = 20
! Max number of water tracers

logical :: tracer_a (0:a_max_trvars)

integer :: a_tr_index(a_max_trvars)
! Index to relative position.
! a_tr_index(n) gives position in jtracer for tracer number N.
! Set in set_atm_pointers.
! a_tr_index(n) is the position, in the list of tracers
! actually present in D1, that tracer number N (in the list
! of all tracers selectable from the user interface) occupies,
! if it is present.
! If tracer number N is absent then A_TR_INDEX(N) is -1.

integer :: a_tr_stashitem(a_max_trvars)
! a_tr_stashitem is set up in set_atm_pointers

integer :: a_tr_lbc_stashitem(a_max_trvars)
! a_tr_lbc_stashitem is set up in inbounda and is only
! referenced if LBC code is active.

integer :: a_tr_active_lbc_index(a_max_trvars)

! Water tracer variables set in module

integer :: n_wtrac = 1 ! required for indexing and allocation

character(len=7), parameter :: norm_class_wtrac   = 'norm   '
character(len=7), parameter :: noniso_class_wtrac = 'non_iso'
character(len=7), parameter :: h218o_class_wtrac  = 'h218o  '
character(len=7), parameter :: hdo_class_wtrac    = 'hdo    '
character(len=7)            :: class_wtrac(max_wtrac)

integer :: i_wtrac_start_non_jls = 0
! Index of the first water tracer which does not go through JULES
! (set in JULES routine wtrac_setup_jls)

!-----------------------------------------------------
! Items set in namelist
!-----------------------------------------------------

logical  ::  l_free_tracer = .false.
! Include tracers in the atmosphere

logical  ::  l_free_tracer_lbc = .false.
! Turn on free tracer LBCs

integer  ::  i_free_tracer(a_max_trvars) = 0
! Specifies which STASH tracer items to be included,
! these will have STASHmaster entries from section 33
!  0 - Do not include
!  1 - Include from dump

integer  ::  i_free_tracer_lbc(a_max_trvars) = 0
! Specifies which tracers have lateral boundary condition
! data in the LBC input file.
! 0 - No LBC data
! 1 - LBC data present

logical :: l_bl_tracer_mix = .false.
! Boundary layer tracer mixing

! ---- PV tracer variables
logical :: l_pv_tracer = .false.

logical :: l_pv_dyn = .false.

logical :: l_pv_tr_rad = .false.

logical :: l_pv_tr_conv = .false.

logical :: l_pv_tr_bl = .false.

logical :: l_calc_pv_full = .false.

integer :: num_dPV = imdi

! ---- Theta tracer variables

logical :: l_theta_tracer= .false.

logical :: l_theta_tr_rad= .false.

logical :: l_theta_tr_bl= .false.

integer :: num_dtheta = imdi

! ---- Water isotopes
logical :: l_wtrac = .false.

logical :: l_h218o_wtrac = .false.

logical :: l_hdo_wtrac = .false.

integer :: n_norm_wtrac = imdi

integer :: n_noniso_wtrac = imdi

integer :: n_noniso_wtrac_jls = imdi

real :: qlimit_h218o_wtrac = rmdi

real :: qlimit_hdo_wtrac = rmdi

namelist /run_free_tracers/ l_free_tracer, l_free_tracer_lbc,                  &
    i_free_tracer, i_free_tracer_lbc, l_bl_tracer_mix,                         &
    l_pv_tracer, l_pv_dyn, l_pv_tr_rad, l_pv_tr_conv, l_pv_tr_bl,              &
    l_calc_pv_full, l_theta_tracer,l_theta_tr_rad,l_theta_tr_bl,               &
    l_wtrac, n_norm_wtrac, n_noniso_wtrac, n_noniso_wtrac_jls,                 &
    qlimit_h218o_wtrac, qlimit_hdo_wtrac, l_h218o_wtrac, l_hdo_wtrac



!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='FREE_TRACERS_INPUTS_MOD'

contains

subroutine print_nlist_run_free_tracers()
use umPrintMgr, only: umPrint
implicit none

character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_FREE_TRACERS'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_free_tracers',                          &
     src='free_tracers_inputs_mod')

write(lineBuffer,*) ' l_free_tracer = ',l_free_tracer
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,*) ' l_free_tracer_lbc  = ',l_free_tracer_lbc
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,*) ' i_free_tracer = ',i_free_tracer
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,*) ' i_free_tracer_lbc = ',i_free_tracer_lbc
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,*) ' l_bl_tracer_mix = ',l_bl_tracer_mix
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_pv_tracer = ',l_pv_tracer
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_pv_dyn = ',l_pv_dyn
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_pv_tr_rad = ',l_pv_tr_rad
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_pv_tr_bl = ',l_pv_tr_bl
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_pv_tr_conv = ',l_pv_tr_conv
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_calc_pv_full = ',l_calc_pv_full
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_theta_tracer = ',l_theta_tracer
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_theta_tr_rad = ',l_theta_tr_rad
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_theta_tr_bl = ',l_theta_tr_bl
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_wtrac = ',l_wtrac
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_h218o_wtrac = ',l_h218o_wtrac
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_hdo_wtrac = ',l_hdo_wtrac
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,I0)') ' n_norm_wtrac = ',n_norm_wtrac
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,I0)') ' n_noniso_wtrac = ',n_noniso_wtrac
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,I0)') ' n_noniso_wtrac_jls = ',n_noniso_wtrac_jls
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,E13.5)') 'qlimit_h218o_wtrac = ',qlimit_h218o_wtrac
call umPrint(lineBuffer,src='free_tracers_inputs_mod')
write(lineBuffer,'(A,E13.5)') 'qlimit_hdo_wtrac = ',qlimit_hdo_wtrac
call umPrint(lineBuffer,src='free_tracers_inputs_mod')

write(lineBuffer,'(A)') '- - - - - - end of namelist - - - - - -'
call umPrint(lineBuffer,src='free_tracers_inputs_mod')


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_free_tracers

!-------------------------------------------------------------------------


!-------------------------------------------------------------------------

subroutine check_run_free_tracers

! Check the water tracer namelist settings and also calculate the total
! number of water tracers (n_wtrac) being used.

use chk_opts_mod,        only: chk_var
implicit none

character(len=*), parameter :: RoutineName='CHECK_RUN_FREE_TRACERS'

character(len=100) :: ChkStr

integer :: ErrorStatus

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0

if (l_wtrac) then

  write(ChkStr,'(A,I0,A)') '[1:',max_wtrac,']'
  call chk_var(n_norm_wtrac, 'n_norm_wtrac', ChkStr)

  write(ChkStr,'(A,I0,A)') '[0:',max_wtrac,']'
  call chk_var(n_noniso_wtrac, 'n_noniso_wtrac', ChkStr)

  ! Check that n_noniso_wtrac_jls is either 0 or equal to n_noniso_wtrac
  write(ChkStr,'(A,I0,A)') '[0,',n_noniso_wtrac,']'
  call chk_var(n_noniso_wtrac_jls, 'n_noniso_wtrac_jls', ChkStr)

  if (l_h218o_wtrac) then
    write(ChkStr,'(A,I0,A)') 'qlimit_h218o_wtrac'
    call chk_var(qlimit_h218o_wtrac, ChkStr, '[-1.0:1.0]')
  end if

  if (l_hdo_wtrac) then
    write(ChkStr,'(A,I0,A)') 'qlimit_hdo_wtrac'
    call chk_var(qlimit_hdo_wtrac, ChkStr, '[-1.0:1.0]')
  end if

  n_wtrac=n_norm_wtrac+n_noniso_wtrac
  if (l_h218o_wtrac) n_wtrac=n_wtrac+1
  if (l_hdo_wtrac) n_wtrac=n_wtrac+1

  write(ChkStr,'(A,I0,A)') '[1:',max_wtrac,']'
  call chk_var(n_wtrac, 'n_wtrac', ChkStr)

else
  n_wtrac=1              ! Required for indexing if l_wtrac=F
  n_noniso_wtrac = 0     ! For safety
  n_noniso_wtrac_jls = 0 ! For safety
end if

if (n_noniso_wtrac == 0) then
  n_noniso_wtrac_jls = 0 ! For safety
end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine check_run_free_tracers

!--------------------------------------------------------------------------

subroutine wtrac_setup_class

!Assign a water tracer class to each member of the water tracer structure,
!based on numbers from the input namelist.

implicit none

character(len=*), parameter  :: routinename='WTRAC_SETUP_CLASS'
integer :: icode ! error code

!DrHook variables
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

integer :: i, i_wtrac

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0

i_wtrac=1
! "Normal" water tracers
do i=1,n_norm_wtrac
  class_wtrac(i_wtrac)= norm_class_wtrac
  i_wtrac=i_wtrac+1
end do

! H2 18O water tracer
if (l_h218o_wtrac) then
  class_wtrac(i_wtrac)= h218o_class_wtrac
  i_wtrac=i_wtrac+1
end if

! HDO tracer
if (l_hdo_wtrac) then
  class_wtrac(i_wtrac)= hdo_class_wtrac
  i_wtrac=i_wtrac+1
end if

! Non-isotropic water tracers
if (n_noniso_wtrac > 0) then
  do i=1,n_noniso_wtrac
    class_wtrac(i_wtrac)= noniso_class_wtrac
    i_wtrac=i_wtrac+1
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine wtrac_setup_class

end module free_tracers_inputs_mod
