! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Completes lifting of the parcel from layer k to k+1

module parcel_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Completes lifting of the parcel from layer k to k+1
!   Calls detrain, term_con and cloud_w
!   Detrain  - carries out the forced detrainment.
!   term_con - tests for any convection which is terminating in layer k+1
!   cloud_w  - carries out the cloud microphysics calculation.
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


character(len=*), parameter, private :: ModuleName='PARCEL_6A_MOD'

contains

subroutine parcel_6a (k, npnts, np_full, nlev, ntra, trlev, n_wtrac,           &
                      ad_on, new_termc,                                        &
                      start_lev,                                               &
                      pk, pkp1, delpkp1, exk, exkp1,                           &
                      thek, thekp1, qek, qekp1,                                &
                      qclek, qclekp1, qcfek, qcfekp1,                          &
                      qsekp1,                                                  &
                      thpk, qpk, qclpk, qcfpk,                                 &
                      rbuoyk, rbuoykp1, rbuoyukp1,                             &
                      watldek, watldekp1, watldpk, watldpkp1,                  &
                      Qlkp1, Qfkp1, Frezkp1,                                   &
                      ekp14, ekp34, amdetk, flxk, flx_init,                    &
                      uek, uekp1, vek, vekp1,                                  &
                      upk, vpk,                                                &
                      trae, wtrac_e,                                           &
                      l_q_interact, l_mom_gk, l_tracer,                        &
                      bgmk, bgmkp1, bwk, bwkp1, bland,                         &
                      ! In/out
                      lcbase, lctop,                                           &
                      thpkp1, qpkp1, qclpkp1, qcfpkp1,                         &
                      tcw, depth, cclwp, lcca,                                 &
                      bterm, blatent, xsbmin,                                  &
                      trap, wtrac_p,                                           &
                      ! Out
                      iccb, icct, prekp1,                                      &
                      thrk, qrk, qrk_wtrac, deltak, flxkp1, flxkp12,           &
                      cca, ccwkp1,                                             &
                      upkp1, vpkp1,                                            &
                      !Indirect indexing
                      idx,ni)


use cv_run_mod,               only: r_det, cpress_term
use water_constants_mod,      only: lc, lf
use cv_derived_constants_mod, only: ls
use planet_constants_mod,     only: cp
use cv_param_mod,             only: term_undil
use yomhook,                  only: lhook, dr_hook
use parkind1,                 only: jprb, jpim
use cloud_w_6a_mod,           only: cloud_w_6a
use thetar_6a_mod,            only: thetar_6a
use thp_det_6a_mod,           only: thp_det_6a
use det_rate_6a_mod,          only: det_rate_6a
use wtrac_conv_mod,           only: l_wtrac_conv, conv_e_wtrac_type,           &
                                    conv_p_wtrac_type
use wtrac_conv_store_mod,     only: conv_old_wtrac_type,                       &
                                    wtrac_alloc_conv_store1
use parcel_wtrac_mod,         only: parcel_wtrac
use cloud_w_wtrac_mod,        only: cloud_w_wtrac

implicit none

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

integer,intent(in) :: k             ! present model layer
integer,intent(in) :: npnts         ! Number of points
integer,intent(in) :: np_full       ! Full vector length
integer,intent(in) :: nlev          ! Number of model levels for calculations
integer,intent(in) :: ntra          ! Number of tracer variables
integer,intent(in) :: trlev         ! No. of model levels on which
                                    ! tracers are included
integer,intent(in) :: n_wtrac       ! Number of water tracers

integer,intent(in) :: ad_on         ! Flag for adaptive detrainment
integer,intent(in) :: new_termc     ! Method of terminating the convective
                                    ! parcel ascent

integer,intent(in) :: start_lev(npnts) ! Level at which convection is initiated

real(kind=real_umphys),intent(in) :: pk(npnts)
                                    ! pressure at mid-point of layer k (Pa)
real(kind=real_umphys),intent(in) :: pkp1(npnts)
                                    ! pressure at mid-point of layer k+1 (Pa)
real(kind=real_umphys),intent(in) :: delpkp1(npnts)
                                    ! pressure difference across layer k+1 (Pa)
real(kind=real_umphys),intent(in) :: exk(npnts)
                                    ! Exner ratio at mid-point of layer k
real(kind=real_umphys),intent(in) :: exkp1(npnts)
                                    ! Exner ratio at mid-point of layer k+1
real(kind=real_umphys),intent(in) :: thek(npnts)
                                    ! Env. pot. temperature in layer k (K)
real(kind=real_umphys),intent(in) :: thekp1(npnts)
                                    ! Env. pot. temperature in layer k+1 (K)
real(kind=real_umphys),intent(in) :: qek(npnts)
                                    ! Env. specific humidity in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qekp1(npnts)
                                    ! Env. spec. humidity in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qclek(npnts)
                                    ! Env. qcl in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qclekp1(npnts)
                                    ! Env. qcl in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qcfek(npnts)
                                    ! Env. qcf in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfekp1(npnts)
                                    ! Env. qcf in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qsekp1(npnts)
                                    ! Env. saturated specific humidity in
                                    ! in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: thpk(npnts)
                                    ! Par. pot. temperature in layer k (K)
real(kind=real_umphys),intent(in) :: qpk(npnts)
                                    ! Par. specific humidity in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qclpk(npnts)
                                    ! Par. qcl in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfpk(npnts)
                                    ! Par. qcf in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: rbuoyk(npnts)
                                    ! Par. buoyancy in layer k (K)
real(kind=real_umphys),intent(in) :: rbuoykp1(npnts)
                                    ! Par. buoyancy in layer k+1 (K)
real(kind=real_umphys),intent(in) :: rbuoyukp1(npnts)
                                    ! Undilute Par. buoyancy in layer k+1 (K)
real(kind=real_umphys),intent(in) :: watldek(npnts)
                                    ! Env. water loading in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: watldekp1(npnts)
                                    ! Env. water loading in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: watldpk(npnts)
                                    ! Par. water loading in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: watldpkp1(npnts)
                                    ! Par. water loading in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: Qlkp1(npnts)
                                    ! Amount of condensation to liquid water
                                    ! in the parcel (kg/kg)
real(kind=real_umphys),intent(in) :: Qfkp1(npnts)
                                    ! Amount of deposition to ice water
                                    ! in the parcel (kg/kg)
real(kind=real_umphys),intent(in) :: Frezkp1(npnts)
                                    ! Amount of freezing from liquid
                                    ! to ice water in the parcel (kg/kg)
real(kind=real_umphys),intent(in) :: ekp14(npnts)
                                    ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: ekp34(npnts)
                                    ! Entrainment coefficient at level k+3/4
                                    ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: amdetk(npnts)
                                    ! Mixing detrainment coefficient at level k
                                    ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: flxk(npnts)
                                    ! Parcel massflux in layer k (Pa/s)
real(kind=real_umphys),intent(in) :: flx_init(npnts)
                                    ! Initial par. massflux at cloud base (Pa/s)
real(kind=real_umphys),intent(in) :: uek(npnts)       ! Env. U in layer k (m/s)
real(kind=real_umphys),intent(in) :: uekp1(npnts)
                                    ! Env. U in layer k+1 (m/s)
real(kind=real_umphys),intent(in) :: vek(npnts)       ! Env. V in layer k (m/s)
real(kind=real_umphys),intent(in) :: vekp1(npnts)
                                    ! Env. V in layer k+1 (m/s)
real(kind=real_umphys),intent(in) :: upk(npnts)       ! Par. U in layer k (m/s)
real(kind=real_umphys),intent(in) :: vpk(npnts)       ! Par. V in layer k (m/s)
real(kind=real_umphys),intent(in) :: trae(np_full,trlev,ntra)
                                            ! Env. tracer content
                                            ! (kg/kg)
type(conv_e_wtrac_type), intent(in) :: wtrac_e(n_wtrac)
                                    ! Structure containing environment
                                    ! water tracer fields

logical,intent(in) :: l_q_interact  ! True if PC2 is switched on
logical,intent(in) :: l_mom_gk      ! Switch for inclusion of Gregory-Kershaw
                                    ! CMT
logical,intent(in) :: l_tracer      ! Switch for tracers

logical,intent(in) :: bgmk(npnts)   ! mask for parcels which are saturated
                                    ! in layer k
logical,intent(in) :: bgmkp1(npnts) ! Mask for parcels which are saturated
                                    ! in layer k+1
logical,intent(in) :: bwk(npnts)    ! mask for parcels which have liquid
                                    ! condensate in layer k
logical,intent(in) :: bwkp1(npnts)  ! mask for parcels which have liquid
                                    ! condensate in layer k+1
logical,intent(in) :: bland(npnts)  ! Land/sea mask

!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------
integer,intent(in out) :: lcbase(npnts)! Lowest conv. cloud base level
integer,intent(in out) :: lctop(npnts) ! Lowest conv. cloud top level

real(kind=real_umphys),intent(in out) :: thpkp1(npnts)
                                     ! Par. pot. temperature in layer k+1 (K)
                                    ! in after entrainment and latent heating
                                    ! out after forced detrainment
real(kind=real_umphys),intent(in out) :: qpkp1(npnts)
                                     ! Par. spec. humidity in layer k+1 (kg/kg)
                                    ! in after entrainment and latent heating
                                    ! out after forced detrainment
real(kind=real_umphys),intent(in out) :: qclpkp1(npnts)
                                     ! Par. qcl in layer k+1 (kg/kg)
                                    ! in after entrainment and latent heating
                                    ! out after forced detrainment
real(kind=real_umphys),intent(in out) :: qcfpkp1(npnts)
                                     ! Par. qcf in layer k+1 (kg/kg)
                                    ! in after entrainment and latent heating
                                    ! out after forced detrainment
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

logical,intent(in out) :: bterm(npnts)   ! Mask for parcels which terminate
                                        ! in layer k+1
logical,intent(in out) :: blatent(npnts) ! Mask for points where latent heat has
                                        ! been released
real(kind=real_umphys),intent(in out) :: xsbmin(npnts)
                                         ! Threshold buoyancy for forced
                                        ! detrainment (K)
real(kind=real_umphys),intent(in out) :: trap(np_full,nlev,ntra)
                                               ! Par. tracer content
                                               ! (kg/kg)

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)
                                    ! Structure containing parcel
                                    ! water tracer fields

!---------------------------------------------------------------------
! Variables which are output
!---------------------------------------------------------------------
integer,intent(in out) :: iccb(npnts)  ! convective cloud base_level
integer,intent(in out) :: icct(npnts)  ! convective cloud top level
! Need to be inout as need to be remembered between subsequent calls to convec2

real(kind=real_umphys),intent(in out) :: prekp1(npnts)
                                       ! precipitation from parcel as it rises
                                       ! from layer k to k+1 (kg/m**2/s)
real(kind=real_umphys),intent(in out) :: thrk(npnts)
                                       ! pot. temperature of forced detrained
                                       ! parcel in layer k (K)
real(kind=real_umphys),intent(in out) :: qrk(npnts)
                                       ! Specific humidity of forced detrained
                                       ! parcel in layer k (kg/kg)
real(kind=real_umphys),intent(in out) :: qrk_wtrac(npnts,n_wtrac)
                                       ! Water tracer vapour content of forced
                                       ! detrained parcel in layer k (kg/kg)
real(kind=real_umphys),intent(in out) :: deltak(npnts)
                                       ! Parcel forced detrainment rate in
                                       ! layer k multiplied by layer thickness
real(kind=real_umphys),intent(in out) :: flxkp1(npnts)
                                       ! parcel massflux in layer k+1 (Pa/s)
real(kind=real_umphys),intent(in out) :: flxkp12(npnts)
                                       ! parcel massflux in layer k+1/2 (Pa/s)
real(kind=real_umphys),intent(in out) :: cca(npnts)
                                       ! convective cloud amount (%)
real(kind=real_umphys),intent(in out) :: ccwkp1(npnts)
                                       ! Total condensate in level k+1 (kg/kg)
real(kind=real_umphys),intent(in out) :: upkp1(npnts)
                                       ! Par. U in layer k+1 (m/s)
real(kind=real_umphys),intent(in out) :: vpkp1(npnts)
                                       ! Par. V in layer k+1 (m/s)

! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

integer :: i, m            ! loop counter

integer :: ndet            ! Compress vector length for the detrainment
                           ! calculation
integer :: ktra, i_wt      ! Loop counter for tracers

integer :: index1(npnts)   ! Index for compress and expand

real(kind=real_umphys) :: Factor1            ! Factor used in update calculation
real(kind=real_umphys) :: Factor2            ! Factor used in update calculation

real(kind=real_umphys), parameter :: tiny_value_condensate = 1.0e-18
                           ! Condensate is removed if below this value

logical :: bdetk(npnts)    ! Mask for points under going forced detrainment

type(conv_old_wtrac_type) :: wtrac_conv_old
                           ! Store water values prior to phase change

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='PARCEL_6A'
!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------

!DIR$ IVDEP
do m=1, ni
  i = idx(m)

  ! ---------------------------------------------------------------------
  !  Calculate mask for those points under going forced detrainment
  !
  !  UM Documentation paper 27
  !  Section (6), equation (23)
  ! ---------------------------------------------------------------------

  !  Calculate XSBMIN (threshold buoyancy for forced detrainment)
  !  Use adaptive detrainment if ADAPTIVE ON flag is set to 1
  if (ad_on  ==  1) then
    ! Overwrite xsbmin with a new adaptive value
    !  (this is the threshold minimum buoyancy for forced detrainment, and is
    !   also the target parcel buoyancy at k+1 after forced detrainment).
    xsbmin(i) = r_det * rbuoyk(i) + (1-r_det)*rbuoykp1(i)
    ! Make sure that adaptive min. buoyancy doesn't fall below 0.0
    xsbmin(i) = max(xsbmin(i),0.0)
  end if
  ! If not using adaptive detrainment, xsbmin is left set to the value passed
  ! down from a higher level (set in the main convection routines).

  !  Set mask of points for forced detrainment
  bdetk(i) = rbuoykp1(i)  <   xsbmin(i)

end do

!Initialise properties of the forced detrained parcel
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  deltak(i) = 0.0
  thrk(i)   = 0.0
  qrk(i)    = 0.0
end do

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
!DIR$ IVDEP
    do m=1, ni
      i = idx(m)
      qrk_wtrac(i,i_wt) = 0.0
    end do
  end do
end if

!----------------------------------------------------------------------
! Compress all input arrays for the forced detrainment calculations
! NB the compression is in the same order as the argument list
!----------------------------------------------------------------------

ndet = 0
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  if (bdetk(i)) then
    ndet = ndet + 1
    index1(ndet) = i
  end if
end do

if (ndet > 0) then

  ! -------------------------------------------------------------------
  !  detrainment calculation
  !
  !  UM Documentation paper 27
  !  Section (6)
  ! -------------------------------------------------------------------

  ! Initialise deltak
!DIR$ IVDEP
  do m=1, ndet
    i = index1(m)
    deltak(i) = 0.0
    thrk(i)   = 0.0
    qrk(i)    = 0.0
  end do

  ! ----------------------------------------------------------------------
  !   Calculate the potential temperature and humidity of the parcel
  !   undergoing forced detrainment
  ! ----------------------------------------------------------------------

  call thetar_6a (npnts, exk, exkp1, pk, thek, thekp1, qek, qekp1,             &
                  qpk, watldek, watldekp1, watldpk, watldpkp1,                 &
                  bwk, bgmk, thrk, qrk, index1, ndet)

  ! ----------------------------------------------------------------------
  !   Calculate the potential temperature and humidity in layer k+1
  !   at the points where detrainment is taking place
  ! ----------------------------------------------------------------------

  call thp_det_6a (npnts, exkp1, pkp1, thekp1, qekp1, watldekp1, watldpkp1,    &
                   xsbmin, bwkp1, bgmkp1, qpkp1, thpkp1, index1, ndet)

  !----------------------------------------------------------------------
  !  calculate forced detrainment rate
  !----------------------------------------------------------------------

  call det_rate_6a (npnts, exkp1,                                              &
                    thek, thekp1, thpk, thpkp1, thrk,                          &
                    qek, qekp1, qpk, qpkp1, qrk,                               &
                    Qlkp1, Qfkp1, Frezkp1,                                     &
                    ekp14, ekp34,                                              &
                    deltak, index1, ndet)

end if     !  ndet > 0

! ----------------------------------------------------------------------
! Update parcel properties at level k+1. NB precipitation is calculated
! later.
! ----------------------------------------------------------------------

!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  flxkp1(i)  = (1.0-amdetk(i))*(1.0-deltak(i))                                 &
               *(1.0+ekp14(i))*(1.0+ekp34(i))*flxk(i)

  ! ---------------------------------------------------------------------
  !  Test for points at which convection terminates in layer k+1
  !  Convection will terminate when the mass flux falls below a small
  !  fraction of the cloud base mass flux or if the forced detrainment
  !  rate is large or if approaching the top of the model.
  !  It will also terminate if bterm has been set to true at a higher
  !  level.
  !----------------------------------------------------------------------

  bterm(i) =  bterm(i)                         .or.                            &
             (flxkp1(i)  <   0.05*flx_init(i)) .or.                            &
             (deltak(i)  >=  0.95)             .or.                            &
             (k+1        ==  nlev)
  if (new_termc == term_undil) then
    bterm(i) =  bterm(i)                       .or.                            &
             (rbuoyukp1(i)<  0.0)
  end if
end do  ! i

!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  if (bterm(i)) then
    ! ---------------------------------------------------------------------
    !   The parcel is terminating so calculate the
    !   the properties at level k+1 for a terminating parcel
    ! ---------------------------------------------------------------------
    bterm(i)   = .true.
    flxkp1(i)  = 0.0
    flxkp12(i) = 0.0
    thpkp1(i)  = 0.0
    qpkp1(i)   = 0.0
    qclpkp1(i) = 0.0
    qcfpkp1(i) = 0.0
    thrk(i)    = thpk(i)
    qrk(i)     = qpk(i)
    deltak(i)  = 1.0
    if (l_tracer) then
      do ktra = 1,ntra
        trap(i,k+1,ktra) = 0.0
      end do
    end if
    if (l_mom_gk) then
      upkp1(i)   = 0.0
      vpkp1(i)   = 0.0
    end if       ! l_mom_gk test

  else
    ! ---------------------------------------------------------------------
    !   The parcel does not terminate if the forced detrainment rate is
    !   not too large and the parcel mass flux is not to small. The parcel
    !   properties at level k+1 are therefore calculated.
    ! ---------------------------------------------------------------------
    !    bterm(i)   = .false.
    Factor1     = 1.0/((1.0+ekp14(i))*(1.0+ekp34(i)))
    Factor2     = Factor1/(1.0-deltak(i))
    !   flxkp1 is calculated earlier
    flxkp12(i)  = (1.0-amdetk(i))*(1.0-deltak(i))*(1.0+ekp14(i))*flxk(i)
    thpkp1(i)   = ( thpk(i) - deltak(i)*thrk(i)                                &
                + (1.0-deltak(i))*ekp14(i)*thek(i)                             &
                + (1.0-deltak(i))*(1.0+ekp14(i))*ekp34(i)*thekp1(i) )          &
                * Factor2                                                      &
                + (lc*Qlkp1(i) + ls*Qfkp1(i) + lf*Frezkp1(i))/(cp*exkp1(i))

    qpkp1(i)    = ( qpk(i) - deltak(i)*qrk(i)                                  &
                + (1.0-deltak(i))*ekp14(i)*qek(i)                              &
                + (1.0-deltak(i))*(1.0+ekp14(i))*ekp34(i)*qekp1(i) )           &
                * Factor2                                                      &
                - Qlkp1(i) - Qfkp1(i)

    if (l_q_interact) then
      !PC2 so entrainment from the environment
      qclpkp1(i) = ( qclpk(i)                                                  &
                 + ekp14(i)*qclek(i)                                           &
                 + (1.0+ekp14(i))*ekp34(i)*qclekp1(i) )                        &
                 * Factor1                                                     &
                 + Qlkp1(i) - Frezkp1(i)
      qcfpkp1(i) = ( qcfpk(i)                                                  &
                 + ekp14(i)*qcfek(i)                                           &
                 + (1.0+ekp14(i))*ekp34(i)*qcfekp1(i) )                        &
                 * Factor1                                                     &
                 + Qfkp1(i) + Frezkp1(i)
    else
      !Not PC2 so no entrainment from the environment
      qclpkp1(i) = qclpk(i)                                                    &
                 * Factor1                                                     &
                 + Qlkp1(i) - Frezkp1(i)
      qcfpkp1(i) = qcfpk(i)                                                    &
                 * Factor1                                                     &
                 + Qfkp1(i) + Frezkp1(i)
    end if ! l_q_interact

    !Due to floating point arithmetic, qclpkp1 or qcfpkp1 could
    !be tiny and/or negative when they should be 0.0. This sometimes
    !happens at phase changes. Therefore reset very small values to zero.
    if (abs(qclpkp1(i))  <  tiny_value_condensate) qclpkp1(i) = 0.0
    if (abs(qcfpkp1(i))  <  tiny_value_condensate) qcfpkp1(i) = 0.0

    blatent(i) = (blatent(i) .or. Qlkp1(i) > 0.0 .or. Qfkp1(i) > 0.0)

    if (l_mom_gk) then   ! l_mom_gk set according to type of CMT required
                         ! This code does Gregory-Kershaw CMT
      upkp1(i)   = ( upk(i)                                                    &
                 + ekp14(i)*uek(i)                                             &
                 + (1.0+ekp14(i))*ekp34(i)*uekp1(i) )                          &
                 * Factor1                                                     &
                 - cpress_term*(uek(i)-uekp1(i))/(1.0+ekp34(i))

      vpkp1(i)   = ( vpk(i)                                                    &
                 + ekp14(i)*vek(i)                                             &
                 + (1.0+ekp14(i))*ekp34(i)*vekp1(i) )                          &
                 * Factor1                                                     &
                 - cpress_term*(vek(i)-vekp1(i))/(1.0+ekp34(i))

    end if       ! l_mom_gk test

    if (l_tracer) then
      do ktra = 1,ntra
        trap(i,k+1,ktra) = ( trap(i,k,ktra)                                    &
                        + ekp14(i)*trae(i,k,ktra)                              &
                        + (1.0+ekp14(i))*ekp34(i)*trae(i,k+1,ktra) )           &
                        * Factor1

      end do
    end if

  end if ! b_term
end do ! m

if (l_wtrac_conv) then

  call parcel_wtrac(k, npnts, ni, n_wtrac, idx, tiny_value_condensate,         &
                    deltak, ekp14, ekp34, qpk, qrk, qclpkp1, qcfpkp1,          &
                    wtrac_e, bterm, wtrac_p, qrk_wtrac)
end if

! ----------------------------------------------------------------------
!  Cloud microphysics calculation
!
!  subroutine cloud_w
!
!  UM Documentation paper 27
!  Section (8), (9)
! ----------------------------------------------------------------------

if (l_wtrac_conv) then
  ! Store water values prior to phase change in cloud_w
  call wtrac_alloc_conv_store1(npnts, ni, idx, qpkp1, qclpkp1, qcfpkp1,        &
                               wtrac_conv_old)
end if

call cloud_w_6a (k, npnts, start_lev,                                          &
                 flxkp1, qclpk, qcfpk,                                         &
                 qsekp1,                                                       &
                 delpkp1,                                                      &
                 bwkp1, bland, bterm, l_q_interact,                            &
                 lcbase, lctop,                                                &
                 qclpkp1, qcfpkp1,                                             &
                 tcw, depth, cclwp, lcca,                                      &
                 iccb, icct, prekp1, cca, ccwkp1, idx, ni)

if (l_wtrac_conv) then
  ! Update water tracers for phase change in cloud_w
  call cloud_w_wtrac (k, npnts, n_wtrac, ni, idx, flxkp1, qclpkp1, qcfpkp1,    &
                      prekp1, wtrac_conv_old, wtrac_p)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine parcel_6a

end module parcel_6a_mod
