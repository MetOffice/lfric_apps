! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates layer dependent constants for layer k, downdraught code
!
module layer_dd_6a_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'LAYER_DD_6A_MOD'
contains

subroutine layer_dd_6a(npnts,k,kct,the_k,the_km1,flx_strt,                     &
                    p_layer_centres,                                           &
                    p_layer_boundaries,                                        &
                    exner_layer_centres,                                       &
                    exner_km12,                                                &
                    pstar,pk,pkm1,delpk,                                       &
                    delpkm1,exk,exkm1,amdetk,ekm14,ekm34,kmin,                 &
                    bddi, recip_pstar )

use cv_param_mod, only:                                                        &
    ae2, ddcoef1, det_lyr

use water_constants_mod, only: tm

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! ------------------------------------------------------------------------------
! Description:
! Calculates layer dependent constants for layer k, downdraught code
!
!  See UM Documentation paper No 27
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts                & ! Number of points
 ,k                    & ! present model layer
 ,kct                    ! convective cloud top layer


logical, intent(in) ::                                                         &
  bddi(npnts)             ! Mask for points where downdraught may initiate

real(kind=real_umphys), intent(in) ::                                          &
  p_layer_centres(npnts,0:kct+1)         & ! Pressure (Pa)
 ,p_layer_boundaries(npnts,0:kct)      & ! Pressure at half level above
                                           ! p_layer_centres (Pa)
 ,exner_layer_centres(npnts,0:kct)       ! exner pressure

real(kind=real_umphys), intent(in) ::                                          &
  pstar(npnts)       & ! Surface pressure (Pa)
 ,exner_km12(npnts)  & ! Exner function at layer k-1/2
 ,flx_strt(npnts)    & ! Updraught mas flux at level where downdraught starts
                       !  (Pa/s)
 ,the_k(npnts)       & ! Pontential temperature of environment
                       ! in layer k (K)
 ,the_km1(npnts)     & ! Pontential temperature of environment
                       ! in layer k-1 (K)
 ,recip_pstar(npnts)   ! 1/pstar (Pa)


!---------------------------------------------------------------------
! Variables which are output
!---------------------------------------------------------------------

integer , intent(out) ::                                                       &
  kmin(npnts)              ! freezing level

real(kind=real_umphys), intent(out) ::                                         &
  pk(npnts)            & ! Pressure of layer k (Pa)
 ,amdetk(npnts)        & ! Mixing detrainment at level k multiplied by
                         ! appropriate layer thickness
 ,ekm14(npnts)         & ! exner ratio at layer k-1/4
 ,ekm34(npnts)           ! exner ratio at layer k-3/4

real(kind=real_umphys), intent(in out) ::                                      &
    pkm1(npnts)          & ! Pressure of layer k-1 (Pa)
   ,delpk(npnts)         & ! Pressure difference across layer k   (Pa)
   ,delpkm1(npnts)       & ! Pressure difference across layer k-1 (Pa)
   ,exkm1(npnts)         & ! Exner ratio at for layer k-1
   ,exk(npnts)             ! Exner ratio at for layer k

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

integer ::                                                                     &
  i                 ! loop counter

real(kind=real_umphys)  ::                                                     &
  ttk        & ! Temperature store at layer k
 ,ttkm1      & ! Temperature store at layer k-1
 ,thkm12     & ! Potential temperature store at layer k-1/2
 ,ttkm12     & ! Temperature store at layer k-1/2
 ,incr_fac   & ! Increment factor for entrainment rates at freezing level
 ,ddcoef2a     ! coefficient used in calculation of downdraught
               ! entrainment rates

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LAYER_DD_6A'

!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Set kmin to initial value
!  Calculate PK, DELPK and EXNER function - If k = KCT then
!  values for previous pass through routine at (k-1)+1 are taken.
!----------------------------------------------------------------------

if (k == kct) then
  do i=1,npnts
    kmin(i) = kct+1   !kmin set to impossible value as unset
    pk(i) = p_layer_centres(i,k)
    delpk(i) =  -(p_layer_boundaries(i,k) - p_layer_boundaries(i,k-1))
    exk(i) = exner_layer_centres(i,k)
  end do
else
  do i=1,npnts
    pk(i) = pkm1(i)
    delpk(i) = delpkm1(i)
    exk(i) = exkm1(i)
  end do
end if

! ---------------------------------------------------------------------
!  Calculate PKM1, DELPKM1
!  Calculate EXNER functions at mid-layer k and k-1, and
!  difference of exner function across layer k
! ---------------------------------------------------------------------

do i=1,npnts
  pkm1(i) = p_layer_centres(i,k-1)
  delpkm1(i) =  -(p_layer_boundaries(i,k-1)                                    &
                - p_layer_boundaries(i,k-2))
  exkm1(i) = exner_layer_centres(i,k-1)
end do

!---------------------------------------------------------------------
! Set DDCOEF2A depending upon which revision of the DD code is used.
!---------------------------------------------------------------------

ddcoef2a=2.0

!
! ---------------------------------------------------------------------
!  Calculate freezing level : Check if freezing level in this layer
! ---------------------------------------------------------------------
!
do i=1,npnts
  if (kmin(i) == kct+1) then  !If kmin not set
    ttk = the_k(i)*exk(i)
    ttkm1 = the_km1(i)*exkm1(i)
    thkm12 = (the_km1(i)+the_k(i))*0.5
    ttkm12 = thkm12*exner_km12(i)
    if (ttkm12  >=  tm .and. ttk  <   tm) then
      kmin(i) = k
    else if (ttkm1  >=  tm .and. ttkm12  <   tm) then
      kmin(i) = k-1
    end if
  end if


  ! ---------------------------------------------------------------------
  !  Calculate entrainment coefficients multiplied by
  !  appropriate layer thickness
  !
  !  Calculate mixing detrainment coefficient multiplied by
  !  appropriate layer thickness
  !
  !  UM DOCUMENTATION PAPER 27
  !  Section (2C), Equation(14)
  ! ---------------------------------------------------------------------

  if (pk(i) <  pstar(i)-det_lyr) then
    ekm14(i) = ae2 * (p_layer_boundaries(i,k-1)-pk(i)) * recip_pstar(i)

    ekm34(i) = ae2 * (pkm1(i)-p_layer_boundaries(i,k-1)) * recip_pstar(i)

    amdetk(i) = (ekm14(i)+ekm34(i)) * (1.0-1.0/ae2)
  else
    ekm14(i) = 0.0
    ekm34(i) = 0.0
    amdetk(i) = delpk(i) /(pstar(i)- p_layer_boundaries(i,k))
  end if

  if (bddi(i) .and. pk(i) <  pstar(i)-det_lyr) then

    if (k == kmin(i)) then
      incr_fac = flx_strt(i)*ddcoef1*recip_pstar(i)
      if (incr_fac >  6.0) incr_fac=6.0
      ekm14(i)  = ekm14(i)*incr_fac
      ekm34(i)  = ekm34(i)*incr_fac
    else
      ekm14(i)  = ekm14(i)*ddcoef2a
      ekm34(i)  = ekm34(i)*ddcoef2a
      amdetk(i) = amdetk(i)*ddcoef2a
    end if

  end if     ! bddi
end do       ! i

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine layer_dd_6a
end module layer_dd_6a_mod
