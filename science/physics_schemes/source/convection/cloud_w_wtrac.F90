! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module cloud_w_wtrac_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Updates water tracers for any precipitation produced in lifting parcel
!   from layer k to k+1 (as calculated in cloud_w)
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'CLOUD_W_WTRAC_MOD'

contains

! Subroutine interface:
subroutine cloud_w_wtrac(k, npnts, n_wtrac, ni, idx, flxkp1, qclpkp1, qcfpkp1, &
                         prekp1, wtrac_conv_old, wtrac_p)

use planet_constants_mod,    only: g
use wtrac_conv_mod,          only: conv_p_wtrac_type
use wtrac_conv_store_mod,    only: conv_old_wtrac_type,                        &
                                   wtrac_dealloc_conv_store1
use wtrac_all_phase_chg_mod, only: wtrac_all_phase_chg

use parkind1,       only: jprb, jpim
use yomhook,        only: lhook, dr_hook

implicit none

! Subroutine arguments:
integer,intent(in) :: k              ! Present model layer number
integer,intent(in) :: npnts          ! Vector length
integer,intent(in) :: n_wtrac        ! No. of water tracers
integer,intent(in) :: ni             ! No. of working points

integer, intent(in) :: idx(npnts)    ! Working points index

real(kind=real_umphys),intent(in) :: flxkp1(npnts)
                                     ! Parcel mass flux in layer k+1 (Pa/s)
real(kind=real_umphys),intent(in) :: qclpkp1(npnts)
                                     ! Parcel liquid condensate mixing ratio
                                     ! in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qcfpkp1(npnts)
                                     ! Parcel frozen condensate mixing ratio
                                     ! in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: prekp1(npnts)
                                     ! precipitation from parcel as it rises
                                     ! from layer k to k+1 (kg/m**2/s)

type(conv_old_wtrac_type), intent(in out) :: wtrac_conv_old
                                     ! Water values before phase change

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)
                                     ! Structure containing parcel
                                     ! water tracer fields
! Local variables
integer :: i, m, i_wt       ! Loop counters

real(kind=real_umphys) :: prekp1_old(ni)
                                     ! Store input prekp1 (kg/m**2/s)
real(kind=real_umphys) :: prekp1_new(ni)
                                     ! Compressed prekp1 after phase change
real(kind=real_umphys) :: qclkp1_new(ni)
                                     ! Compressed qclkp1 after phase change
real(kind=real_umphys) :: qcfkp1_new(ni)
                                     ! Compressed qcfkp1 after phase change
real(kind=real_umphys) :: qcl_change(ni)
                                     ! Amount of liquid changing phase
real(kind=real_umphys) :: qcf_change(ni)
                                     ! Amount of ice changing phase

! Work arrays for water tracer
real(kind=real_umphys) :: qclpkp1_wtrac(ni,n_wtrac)
                                     ! Parcel water tracer
                                     ! liquid condensate (kg/kg) at layer k+1
real(kind=real_umphys) :: qcfpkp1_wtrac(ni,n_wtrac)
                                     ! Parcel water tracer
                                     ! ice condensate (kg/kg) at layer k+1
real(kind=real_umphys) :: prekp1_wtrac(ni,n_wtrac)
                                     ! Water tracer precip from parcel as it
                                     ! rises from layer k to k+1 (kg/m**2/s)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CLOUD_W_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up compressed working arrays for water tracer k+1 values and initialise
! water tracer precip to zero
do i_wt = 1, n_wtrac
  do m = 1, ni
    i = idx(m)
    qclpkp1_wtrac(m,i_wt) = wtrac_p(i_wt)%qcl(i,k+1)
    qcfpkp1_wtrac(m,i_wt) = wtrac_p(i_wt)%qcf(i,k+1)
    prekp1_wtrac(m,i_wt)  = 0.0
  end do
end do

! Set up compressed arrays and initialise precip to zero and calculate
! amount of water changing phase
do m = 1, ni
  i = idx(m)
  prekp1_old(m)  = 0.0
  qcl_change(m)  = max(wtrac_conv_old%qclpkp1(m) -  qclpkp1(i), 0.0)
  qcf_change(m)  = max(wtrac_conv_old%qcfpkp1(m) -  qcfpkp1(i), 0.0)
  qclkp1_new(m)  = qclpkp1(i)
  qcfkp1_new(m)  = qcfpkp1(i)
  prekp1_new(m)  = prekp1(i)
end do

! -----------------------------------------------------------------
! Update water tracers for phase change
! -----------------------------------------------------------------

! Rain formation (liquid -> precip)
call wtrac_all_phase_chg(ni, n_wtrac, wtrac_conv_old%qclpkp1, prekp1_old,      &
                         qcl_change, qclkp1_new, prekp1_new,                   &
                         'liq', 'rai', 'one_way',                              &
                         qclpkp1_wtrac, prekp1_wtrac)

! Snow formation (ice -> precip)
call wtrac_all_phase_chg(ni, n_wtrac, wtrac_conv_old%qcfpkp1, prekp1_old,      &
                         qcf_change, qcfkp1_new, prekp1_new,                   &
                         'ice', 'sno', 'one_way',                              &
                         qcfpkp1_wtrac, prekp1_wtrac)

! Conversion of water tracer precip into rate
do i_wt = 1, n_wtrac
  do m = 1, ni
    i = idx(m)
    prekp1_wtrac(m,i_wt) =  prekp1_wtrac(m,i_wt) * flxkp1(i) /g
  end do
end do

! Update k+1 values in main structure
do i_wt = 1, n_wtrac
  do m = 1, ni
    i = idx(m)
    wtrac_p(i_wt)%qcl(i,k+1)    = qclpkp1_wtrac(m,i_wt)
    wtrac_p(i_wt)%qcf(i,k+1)    = qcfpkp1_wtrac(m,i_wt)
    wtrac_p(i_wt)%precip(i,k+1) = prekp1_wtrac(m,i_wt)
  end do
end do

! Deallocate temporary arrays
call wtrac_dealloc_conv_store1(wtrac_conv_old)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine cloud_w_wtrac

end module cloud_w_wtrac_mod
