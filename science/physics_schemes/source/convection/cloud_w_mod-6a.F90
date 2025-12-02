! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Convective cloud microphysics routine
!
module cloud_w_6a_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Calculates preciptation produced in lifting parcel from layer k to k+1
!
!   Calls CON_RAD to calculate parameters for radiation calculation
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


character(len=*), parameter, private :: ModuleName='CLOUD_W_6A_MOD'

contains

! Subroutine Interface:

subroutine cloud_w_6a (k, npnts, start_lev,                                    &
                       flxkp1, qclpk, qcfpk,                                   &
                       qsekp1,                                                 &
                       delpkp1,                                                &
                       bwkp1, bland, bterm, l_q_interact,                      &
                       lcbase, lctop,                                          &
                       qclpkp1, qcfpkp1,                                       &
                       tcw, depth, cclwp, lcca,                                &
                       iccb, icct, prekp1, cca, ccwkp1, idx, ni)

use planet_constants_mod, only: g

use cv_run_mod, only:                                                          &
    mparwtr, qlmin, fac_qsat, ccw_for_precip_opt

use parkind1,       only: jprb, jpim
use yomhook,        only: lhook, dr_hook
use con_rad_6a_mod, only: con_rad_6a

implicit none

! Subroutine arguments

!----------------------------------------------------------------------
! Variables that are input
!----------------------------------------------------------------------
integer,intent(in) :: k                   ! Present model layer number
integer,intent(in) :: npnts               ! Vector length

integer,intent(in) :: start_lev(npnts)    ! Level at which convection initiated

real(kind=real_umphys),intent(in) :: flxkp1(npnts)
                                    ! Parcel mass flux in layer k+1 (Pa/s)
real(kind=real_umphys),intent(in) :: qclpk(npnts)
                                    ! Par. qcl in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfpk(npnts)
                                    ! Par. qcf in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qsekp1(npnts)
                                    ! Env. saturated specific humidity in
                                    ! in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: delpkp1(npnts)
                                    ! pressure difference across layer k+1 (Pa)

logical,intent(in) :: bwkp1(npnts)  ! mask for parcels which have liquid
                                    ! condensate in layer k+1
logical,intent(in) :: bland(npnts)  ! Land/sea mask
logical,intent(in) :: bterm(npnts)  ! Mask for parcels which terminate
                                    ! in layer k+1
logical,intent(in) :: l_q_interact  ! True if PC2 is switched on

!----------------------------------------------------------------------
! Variables that are input and output
!----------------------------------------------------------------------
integer,intent(in out) :: lcbase(npnts)! Lowest conv. cloud base level
integer,intent(in out) :: lctop(npnts) ! Lowest conv. cloud top level

real(kind=real_umphys),intent(in out) :: qclpkp1(npnts)
                                     ! Parcel liquid condensate mixing ratio
                                    ! in layer k+1 (kg/kg)
                                    ! in:  before precipitation
                                    ! out: after precipitation
real(kind=real_umphys),intent(in out) :: qcfpkp1(npnts)
                                     ! Parcel liquid condensate mixing ratio
                                    ! in layer k+1 (kg/kg)
                                    ! in:  before precipitation
                                    ! out: after precipitation
real(kind=real_umphys),intent(in out) :: tcw(npnts)
                                     ! Total condensed water (kg/m**2/s)
                                    ! in summed to layer k
                                    ! out summed to layer k+1
real(kind=real_umphys),intent(in out) :: depth(npnts)
                                     ! Depth of convective cloud (m)
                                    ! in summed to layer k
                                    ! out summed to layer k+1
real(kind=real_umphys),intent(in out) :: cclwp(npnts)
                                     ! Condensed water path (kg/m**2)
                                    ! in summed to layer k
                                    ! out summed to layer k+1
real(kind=real_umphys),intent(in out) :: lcca(npnts)
                                     ! Lowest conv. cloud amount (%)

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
integer,intent(in out) :: iccb(npnts)  ! convective cloud base_level
integer,intent(in out) :: icct(npnts)  ! convective cloud top level
! Need to be inout as need to be remembered between subsequent calls to convec2

real(kind=real_umphys),intent(in out) :: prekp1(npnts)
                                       ! precipitation from parcel as it rises
                                       ! from layer k to k+1 (kg/m**2/s)
real(kind=real_umphys),intent(in out) :: cca(npnts)
                                       ! convective cloud amount (%)
real(kind=real_umphys),intent(in out) :: ccwkp1(npnts)
                                       ! Total condensate in level k+1 (kg/kg)


! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni


!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

integer :: i,m              ! Loop counter

real(kind=real_umphys) :: mparmult
                            ! Factor used to multiply mparwtr, value between
                            ! 1. and ~3. being 1.0 for deeper clouds.
real(kind=real_umphys) :: tmp_ccwkp1(npnts)
                            ! Total condensate in level k+1 before precipitation
                            ! (kg/kg)
real(kind=real_umphys) :: ccwk(npnts)
                            ! Total condensate in level k (kg/kg)

real(kind=real_umphys) :: xmin(npnts)
                            ! Maxmimum amount of cloud water retained by the
                            ! parcel before the parcel precipitates (kg/kg)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CLOUD_W_6A'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! Store convective cloud liquid water before Precipitation
!-------------------------------------------------------------------------------
!
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  tmp_ccwkp1(i) = qclpkp1(i) + qcfpkp1(i)  ! Total condensate in level k+1
                                           ! before precipitation
  ccwk(i)       = qclpk(i)   + qcfpk(i)    ! Total condensate in level k
  prekp1(i)     = 0.0                      ! initialise precipitation to zero
end do

!-------------------------------------------------------------------------------
! Calculate Precipitation from layer k+1 and adjust cloud water
!-------------------------------------------------------------------------------

select case (ccw_for_precip_opt)

case (4)            ! xmin profile based on qsat with user defined
                    ! minimum, maximum and qsat scaling.
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    xmin(i) = max(min(mparwtr, fac_qsat*qsekp1(i)), qlmin)

  end do

case (3)            ! Manoj's function for congestus
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    mparmult = 1.5 + 0.5*tanh((3000.0 -depth(i))/800.0)
    xmin(i)  = min(mparwtr*mparmult, fac_qsat*qsekp1(i))

    ! If a land point and liquid water in the layer
    ! increase the minimum cloud water for Precipitation
    ! The reasons for this are an attempt to take some account of
    ! more aerosols over land leading to more small cloud drops and
    ! therefore more cloud water before Precipitation.

    if (bwkp1(i) .and. bland(i)) xmin(i) =xmin(i)*2.0

    ! limit max value to 0.003 kg/kg
    xmin(i) = min(xmin(i),0.003)

    if (l_q_interact) then   ! PC2
      ! Limit xmin(i)
      xmin(i) = max(xmin(i), qlmin)
    end if

  end do

case (2)            ! Manoj's function dependent on depth of cloud
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    mparmult = 2.0 + 1.0*tanh((1500.0 -depth(i))/1000.0)
    xmin(i)  = min(mparwtr*mparmult, fac_qsat*qsekp1(i))

    ! If a land point and liquid water in the layer
    ! increase the minimum cloud water for Precipitation
    ! The reasons for this are an attempt to take some account of
    ! more aerosols over land leading to more small cloud drops and
    ! therefore more cloud water before Precipitation.

    ! If a land point and liquid water in the layer
    if (bwkp1(i) .and. bland(i)) xmin(i) = xmin(i)*2.0

    if (l_q_interact) then
      ! Limit xmin
      xmin(i) = max(xmin(i), qlmin)
    end if

  end do

case (1)            ! No test on a critical depth for Precipitation
                    ! Option in use for HadGEM2
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    xmin(i) = min(mparwtr, fac_qsat*qsekp1(i))

    ! If a land point and liquid water in the layer
    if (bwkp1(i) .and. bland(i)) xmin(i) = xmin(i)*2.0

    if (l_q_interact) then
      ! Limit xmin
      xmin(i) = max(xmin(i), qlmin)
    end if

  end do

end select      ! test on value of ccw_for_precip_opt

!DIR$ IVDEP
do m=1, ni
  i = idx(m)

  ! Precipitate if cloud water in the layer > xmin

  if ( bwkp1(i) ) then  !If condensate is liquid
    if ( qclpkp1(i) > xmin(i) ) then
      prekp1(i)  = (qclpkp1(i) - xmin(i)) * flxkp1(i) /g
      qclpkp1(i) = xmin(i)
    end if
  else !Condensate is frozen
    if ( qcfpkp1(i) > xmin(i) ) then
      prekp1(i)  = (qcfpkp1(i) - xmin(i)) * flxkp1(i) /g
      qcfpkp1(i) = xmin(i)
    end if
  end if

  !Update the convective cloud water. Do not permit negative CCW.
  ccwkp1(i) = max(qclpkp1(i) + qcfpkp1(i), 0.0)

end do

!-------------------------------------------------------------------------------
! Calculate convective cloud base, convective cloud top, total condensed
! water/ice and convective cloud amount
!
! subroutine CON_RAD - now called after rainout has occurred
! UM Documentation paper 27
! Section (9)
!-------------------------------------------------------------------------------

call con_rad_6a (k, npnts, start_lev,                                          &
                 ccwk, ccwkp1, tmp_ccwkp1, flxkp1, delpkp1,                    &
                 l_q_interact, bterm,                                          &
                 tcw, cclwp, lcca, lcbase, lctop,                              &
                 iccb, icct, cca, idx, ni)

!-------------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine cloud_w_6a
end module cloud_w_6a_mod
