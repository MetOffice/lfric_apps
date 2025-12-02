! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the effect of convection upon the large-scale atmosphere
!
module environ_6a_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Calculate the effect of convection upon the large-scale atmosphere
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

character(len=*), parameter, private :: ModuleName='ENVIRON_6A_MOD'

contains

! Subroutine Interface:
subroutine environ_6a (k, npnts, np_full, nlev, ntra, trlev, n_wtrac,          &
                    timestep, delpk, delpkp1, delp_uv_k, delp_uv_kp1,          &
                    exk, exkp1,                                                &
                    thek, thekp1, qek, qekp1,                                  &
                    qclek, qclekp1, qcfek, qcfekp1,                            &
                    thpk, qpk,                                                 &
                    qclpk, qcfpk,                                              &
                    thrk, qrk, qrk_wtrac,                                      &
                    ekp14, amdetk, deltak, flxk,                               &
                    uekp1, vekp1,                                              &
                    upk, upkp1, vpk, vpkp1,                                    &
                    trae, trap,                                                &
                    l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,         &
                    blowst, bterm, wtrac_p,                                    &
                    ! In/out
                    dthek, dqek, dqclek, dqcfek,                               &
                    eflux_u_ud, eflux_v_ud,                                    &
                    duek, dvek,                                                &
                    dtrae, wtrac_e,                                            &
                    ! Out
                    dthekp1, dqekp1, dqclekp1, dqcfekp1,                       &
                    duekp1, dvekp1,tnuc_nlcl,                                  &
                    !Indirect indexing
                    idx,ni)

use water_constants_mod,      only: lf, tm
use cv_derived_constants_mod, only: lcrcp, lfrcp, lsrcp
use planet_constants_mod,     only: cp

use cloud_inputs_mod,         only: i_pc2_conv_coupling,                       &
                                    starticeTKelvin, alliceTdegC
use science_fixes_mod,        only: l_fix_pc2_cnv_mix_phase
use mphys_inputs_mod,         only: l_progn_tnuc
use lsprec_mod,               only: zerodegc
use wtrac_conv_mod,           only: l_wtrac_conv, conv_e_wtrac_type,           &
                                    conv_p_wtrac_type
use environ_wtrac_mod,        only: environ_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer,intent(in) :: k             ! present model layer
integer,intent(in) :: npnts         ! Number of points
integer,intent(in) :: np_full       ! Full vector length
integer,intent(in) :: nlev          ! Number of model levels for calculations
integer,intent(in) :: ntra          ! Number of tracer variables
integer,intent(in) :: trlev         ! No. of model levels on which
                                    ! tracers are included
integer,intent(in) :: n_wtrac       ! Number of water tracer variables

real(kind=real_umphys),intent(in) :: timestep        ! Timestep
real(kind=real_umphys),intent(in) :: delpk(npnts)
                                    ! pressure difference across layer k (Pa)
real(kind=real_umphys),intent(in) :: delpkp1(npnts)
                                    ! pressure difference across layer k+1 (Pa)
real(kind=real_umphys),intent(in) :: delp_uv_k(npnts)
                                    ! pressure difference across UV
                                    ! layer k (Pa)
real(kind=real_umphys),intent(in) :: delp_uv_kp1(npnts)
                                      ! pressure difference across UV
                                    ! layer k+1 (Pa)
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
real(kind=real_umphys),intent(in) :: thpk(npnts)
                                    ! Par. pot. temperature in layer k (K)
real(kind=real_umphys),intent(in) :: qpk(npnts)
                                    ! Par. specific humidity in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qclpk(npnts)
                                    ! Par. qcl in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfpk(npnts)
                                    ! Par. qcf in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: thrk(npnts)
                                    ! pot. temperature of forced detrained
                                    ! parcel in layer k (K)
real(kind=real_umphys),intent(in) :: qrk(npnts)
                                    ! Specific humidity of forced detrained
                                    ! parcel in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qrk_wtrac(npnts,n_wtrac)
                                    ! Water tracer specific humidity of forced
                                    ! detrained parcel in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: ekp14(npnts)
                                    ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: amdetk(npnts)
                                    ! Mixing detrainment coefficient at level k
                                    ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: deltak(npnts)
                                    ! Parcel forced detrainment rate in
                                    ! layer k multiplied by layer thickness
real(kind=real_umphys),intent(in) :: flxk(npnts)
                                    ! Parcel massflux in layer k (Pa/s)
real(kind=real_umphys),intent(in) :: uekp1(npnts)
                                    ! Env. U in layer k+1 (m/s)
real(kind=real_umphys),intent(in) :: vekp1(npnts)
                                    ! Env. V in layer k+1 (m/s)
real(kind=real_umphys),intent(in) :: upk(npnts)       ! Par. U in layer k (m/s)
real(kind=real_umphys),intent(in) :: upkp1(npnts)
                                    ! Par. U in layer k+1 (m/s)
real(kind=real_umphys),intent(in) :: vpk(npnts)       ! Par. V in layer k (m/s)
real(kind=real_umphys),intent(in) :: vpkp1(npnts)
                                    ! Par. V in layer k+1 (m/s)
real(kind=real_umphys),intent(in) :: trae(np_full,trlev,ntra)
                                            ! Env. tracer content
                                            ! (kg/kg)
real(kind=real_umphys),intent(in) :: trap(np_full,nlev,ntra)
                                          ! Par. tracer content
                                          ! (kg/kg)
real(kind=real_umphys),intent(in) :: tnuc_nlcl(npnts)
                                    !nucleation temperature as function of dust
                                    !indexed using nlcl(deg cel)
real(kind=real_umphys) :: n_starticetK
                      !new detrainment temperature (startice) to be used
                      !with dust-tnuc switch(Kelvin)
real(kind=real_umphys) :: n_allicedegC
                      !new detrainment temperature (allice) to be used
                      !with dust-tnuc switch(deg cel)
logical,intent(in) :: l_q_interact  ! True if PC2 is switched on
logical,intent(in) :: l_mom_gk      ! Switch for inclusion of Gregory-Kershaw
                                    ! CMT
logical,intent(in) :: l_mom_gk_stable ! Switch for stabilized Gregory-Kershaw
                                      ! CMT
logical,intent(in) :: l_tracer      ! Switch for tracers

logical,intent(in) :: blowst(npnts) ! mask for those points at which stability
                                    ! is low enough for convection to occur
logical,intent(in) :: bterm(npnts)  ! Mask for parcels which terminate
                                    ! in layer k+1

type(conv_p_wtrac_type), intent(in) :: wtrac_p(n_wtrac)
                                        ! Structure containing parcel
                                        ! water tracer fields

!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------

! Convection increments to model fields at level k
! in before processes at level k. This may be non-zero because of smoothed
! forced detrainment and the initial perturbation
! out after processes at level k
real(kind=real_umphys),intent(in out) :: dthek(npnts)
                                     ! Increment to p. temperature
                                    ! in layer k (K/s)
real(kind=real_umphys),intent(in out) :: dqek(npnts)
                                     ! Increment to spec. humidity
                                    ! in layer k (kg/kg/s)
real(kind=real_umphys),intent(in out) :: dqclek(npnts)
                                     ! Increment to qcl in layer k (kg/kg/s)
real(kind=real_umphys),intent(in out) :: dqcfek(npnts)
                                     ! Increment to qcf in layer k (kg/kg/s)

real(kind=real_umphys),intent(in out) :: eflux_u_ud(npnts)
                                         ! Eddy flux of momentum to UD
                                    ! in at bottom of layer
                                    ! out at top of layer
real(kind=real_umphys),intent(in out) :: eflux_v_ud(npnts)
                                         ! Eddy flux of momentum to UD
                                    ! in at bottom of layer
                                    ! out at top of layer
! Convection increments to winds and tracers at level k
! in before processes at level k. This may be non-zero because of smoothed
! forced detrainment and the initial perturbation
! out after processes at level k
real(kind=real_umphys),intent(in out) :: duek(npnts)
                                     ! Increment to U in layer k (m/s**2)
real(kind=real_umphys),intent(in out) :: dvek(npnts)
                                     ! Increment to V in layer k (m/s**2)
real(kind=real_umphys),intent(in out) :: dtrae(np_full,nlev, ntra)
                                                  ! Increment to tracer
                                                  ! (kg/kg/s)

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                        ! Structure containing environment
                                        ! water tracer fields

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
! Convection increments to model fields at level k+1
real(kind=real_umphys),intent(in out) :: dthekp1(npnts)
                                       ! Increment to p. temperature
                                       ! in layer k+1 (K/s)
real(kind=real_umphys),intent(in out) :: dqekp1(npnts)
                                       ! Increment to spec. humidity
                                       ! in layer k+1 (kg/kg/s)
real(kind=real_umphys),intent(in out) :: dqclekp1(npnts)
                                       ! Increment to qcl in layer k+1 (kg/kg/s)
real(kind=real_umphys),intent(in out) :: dqcfekp1(npnts)
                                       ! Increment to qcf in layer k+1 (kg/kg/s)
real(kind=real_umphys),intent(in out) :: duekp1(npnts)
                                       ! Increment to U in layer k+1 (m/s**2)
real(kind=real_umphys),intent(in out) :: dvekp1(npnts)
                                       ! Increment to V in layer k+1 (m/s**2)

! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni


!-----------------------------------------------------------------------
! Variables that are defined locally
!-----------------------------------------------------------------------

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='ENVIRON_6A'

integer :: i,m            ! loop counter
integer :: ktra           ! Loop counter for tracers

! Smoothed forced detrainment variables
real(kind=real_umphys) :: tmp_fd_dthek
                          ! F. detrainment P.temp inc across levels k and k+1
real(kind=real_umphys) :: tmp_fd_dqek
                          ! F. detrainment humidity inc across levels k and k+1
real(kind=real_umphys) :: tmp_fd_dqclek
                          ! F. detrainment qcl inc across levels k and k+1
real(kind=real_umphys) :: tmp_fd_dqcfek
                          ! F. detrainment qcf across levels k and k+1
real(kind=real_umphys) :: tmp_fd_dtraek
                          ! F. detrainment tracer across levels k and k+1

real(kind=real_umphys) :: flxbydpk(npnts)
                          ! mass flux divided by layer thickness (1/s)
real(kind=real_umphys) :: flx_u_kp0p5
                          ! Flux of zonal momentum in cloud at top of current
                          ! uv layer
real(kind=real_umphys) :: flx_v_kp0p5
                          ! Flux of meridional momentum in cloud at top of
                          ! current uv layer
real(kind=real_umphys) :: tmp_dqclek
                          ! Storage space for liquid condensate rate.
real(kind=real_umphys) :: tmp_dqcfek
                          ! Storage space for frozen condensate rate.
real(kind=real_umphys) :: dth_detcond(npnts)
                          ! Used to calculate the contribution to the theta
                          ! increment from the evaporation of detrained
                          ! 0.0 if PC2 is on.
real(kind=real_umphys) :: dq_detcond(npnts)
                          ! Used to calculate the contribution to the humidity
                          ! increment from the evaporation of detrained
                          ! 0.0 if PC2 is on.
real(kind=real_umphys) :: frac_icek         ! Fraction of ice water at level k
real(kind=real_umphys) :: frac_liqk
                          ! Fraction of liquid water at level k
real(kind=real_umphys) :: frac_icekp1       ! Fraction of ice water at level k+1
real(kind=real_umphys) :: frac_liqkp1
                          ! Fraction of liquid water at level k+1
real(kind=real_umphys) :: tmp_dqcfekp1
                          ! Storage space for ice condensate rate at level k+1
real(kind=real_umphys) :: tmp_dqclekp1
                          ! Storage space for liquid condensate rate, level k+1

real(kind=real_umphys),allocatable :: up_ref(:)
                              ! Par. U to be referred for momentum flux
                              ! at k+1/2
real(kind=real_umphys),allocatable :: vp_ref(:)
                              ! Par. V to be referred for momentum flux
                              ! at k+1/2

! Phase change amounts needed for water tracer
real(kind=real_umphys),allocatable :: dqclek_frz(:)
                              ! Freezing rate of liq condensate at k
real(kind=real_umphys),allocatable :: dqcfek_mlt(:)
                              ! Melting rate of ice condensate at k
real(kind=real_umphys),allocatable :: dqclekp1_frz(:)
                              ! Freezing rate of liq condensate at k+1
real(kind=real_umphys),allocatable :: dqcfekp1_mlt(:)
                              ! Melting rate of ice condensate at k+1
real(kind=real_umphys),allocatable :: qclek_temp(:)
                              ! Temporary updated qclek field at k
real(kind=real_umphys),allocatable :: qcfek_temp(:)
                              ! Temporary updated qcfek field at k

! Parameters

real(kind=real_umphys),parameter :: a_smth    = 0.5
                                      ! Parameter determining the weighting
                                      ! between the forced detrainment
                                      ! increments at k and k-1
real(kind=real_umphys),parameter :: bridge_temp = 10.0
                                      ! Parameter determining the
                                      ! bridging value between phase
                                      ! detrainment temperatures
                                      ! while using dust-tnuc scheme

! Several options are available:
integer,parameter :: pc2_conv_original = 1
! As described in Wilson et al (2008). Condensate increments from
! detrainment and subsidence advection (combined) are turned into cloud
! fraction increments using the inhomogeneous forcing method. Note that
! the abrupt change in phase of condensate in the convective plume is
! replicated in the condensate increments due to convection.
integer,parameter :: pc2_conv_maxincldqc = 2
! Unpublished. As above + Protects (in ni_conv_ctl) against
! generation of inconsistently low cloud fraction implying
! very high in-cloud condensate amounts.
!integer,parameter :: pc2_conv_smooth_liqice = 3
! Unpublished. As above + The phase of the detrained condensate varies
! smoothly according to the ambient temperature rather than having
! an abrupt change.

!---------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_q_interact) then
  ! PC2 is on and therefore the detrained condensate is not evaporated
  ! and does not affect the theta and q increments
  ! Note, water tracer code below assumes that dq_detcond=0 if using PC2
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    dth_detcond(i)  = 0.0
    dq_detcond(i)   = 0.0
  end do
else
  ! PC2 is off and therefore the detrained condensate is evaporated
  ! and will affect the theta and q increments
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    dth_detcond(i)  = (lcrcp*qclpk(i) + lsrcp*qcfpk(i))/exk(i)
    dq_detcond(i)   =        qclpk(i) +       qcfpk(i)
  end do
end if

if (l_wtrac_conv) then
  ! Allocate and initialise arrays for storing melt/freezing amounts
  allocate(dqclek_frz(npnts))
  allocate(dqcfek_mlt(npnts))
  allocate(dqclekp1_frz(npnts))
  allocate(dqcfekp1_mlt(npnts))
  allocate(qclek_temp(npnts))
  allocate(qcfek_temp(npnts))

  do i = 1, npnts
    dqclek_frz(i)    = 0.0
    dqcfek_mlt(i)    = 0.0
    dqclekp1_frz(i)  = 0.0
    dqcfekp1_mlt(i)  = 0.0
    qclek_temp(i)    = qclek(i)
    qcfek_temp(i)    = qcfek(i)
  end do

end if


!DIR$ IVDEP
do m=1, ni
  i = idx(m)

  !----------------------------------------------------------------------
  ! Calculate parcel mass flux divided by the thickness of layer k.
  ! This value is used in several places in the subroutine.
  !----------------------------------------------------------------------
  flxbydpk(i)   = flxk(i)/delpk(i)

  !-----------------------------------------------------------------------
  !         Potential temperature increment
  !-----------------------------------------------------------------------
  !                Smoothed forced detrainment term
  tmp_fd_dthek  = deltak(i) * (1.0-amdetk(i))                                  &
                * ( thrk(i) - thek(i) - dth_detcond(i) )

  !           Theta increment at level k
  dthek(i)      = dthek(i) + flxbydpk(i)                                       &
  !                Compensating subsidence
                  * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))               &
                  * ( thekp1(i)-thek(i) )                                      &
  !                Smoothed forced detrainment
                  + a_smth*tmp_fd_dthek                                        &
  !                Mixing detrainment
                  + amdetk(i)                                                  &
                  * ( thpk(i) - thek(i) - dth_detcond(i) ) )

  !           Theta increment at level k+1 due
  !           to smoothed  forced detrainment
  dthekp1(i)    = flxbydpk(i) * (1.0-a_smth)                                   &
                * exk(i)/exkp1(i)                                              &
                * delpk(i)/delpkp1(i) * tmp_fd_dthek

  !-----------------------------------------------------------------------
  !         Humidity increment
  !-----------------------------------------------------------------------
  !                Smoothed forced detrainment term
  tmp_fd_dqek   = deltak(i) * (1.0-amdetk(i))                                  &
                * ( qrk(i) - qek(i) + dq_detcond(i) )

  !           q increment at level k
  dqek(i)       = dqek(i) + flxbydpk(i)                                        &
  !                Compensating subsidence
                  * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))               &
                  * ( qekp1(i)-qek(i) )                                        &
  !                Smoothed forced detrainment
                  + a_smth*tmp_fd_dqek                                         &
  !                Mixing detrainment
                  + amdetk(i)                                                  &
                  * ( qpk(i) - qek(i) + dq_detcond(i) ) )

  !           q increment at level k+1 due
  !           to smoothed  forced detrainment
  dqekp1(i)     = flxbydpk(i) * (1.0-a_smth)                                   &
                * delpk(i)/delpkp1(i)*tmp_fd_dqek
end do



if (l_q_interact) then
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    !-----------------------------------------------------------------------
    !         PC2 is on. Set the qcl and qcf increments appropriately
    !-----------------------------------------------------------------------
    if (l_progn_tnuc) then   !prognostic tnuc
      n_starticetK = (tnuc_nlcl(i) + zerodegc)
      n_allicedegC = (tnuc_nlcl(i) - bridge_temp)  !bridge temp detrainment
    else
      n_starticetK = starticeTKelvin
      n_allicedegC = alliceTdegC
    end if
    frac_icekp1   = (thekp1(i)*exkp1(i)-n_starticetK) /                        &
                        (n_allicedegC-(n_starticetK-tm))
    frac_icekp1   = min(max(0.0,frac_icekp1),1.0)
    frac_liqkp1   = 1.0-frac_icekp1

    !-----------------------------------------------------------------------
    !         Liquid condensate increment
    !         A temporary variable is used because PC2 modifies the
    !         increment
    !-----------------------------------------------------------------------
    !                Smoothed forced detrainment term
    tmp_fd_dqclek = deltak(i) * (1.0-amdetk(i))                                &
                  * ( qclpk(i)-qclek(i) )

    !           qcl increment at level k
    !           previous contribution at level k not included here
    !           but is added later because of PC2 messing about with
    !           this increment.
    tmp_dqclek    = flxbydpk(i)                                                &
    !                  Compensating subsidence
                      * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))           &
                      * ( qclekp1(i)-qclek(i) )                                &
    !                  Smoothed forced detrainment
                      + a_smth*tmp_fd_dqclek                                   &
    !                  Mixing detrainment
                      + amdetk(i)                                              &
                      * ( qclpk(i) - qclek(i) ) )

    if (l_wtrac_conv) qclek_temp(i) = qclek(i) + tmp_dqclek*timestep
    !-----------------------------------------------------------------------
    !         Frozen condensate increment
    !         A temporary variable is used because PC2 modifies the
    !         increment
    !-----------------------------------------------------------------------
    !                Smoothed forced detrainment term
    tmp_fd_dqcfek = deltak(i) * (1.0-amdetk(i))                                &
                  * ( qcfpk(i)-qcfek(i) )

    !           qcf increment at level k
    !           previous contribution at level k not included here
    !           but is added later because of PC2 messing about with
    !           this increment.
    tmp_dqcfek    = flxbydpk(i)                                                &
    !                  Compensating subsidence
                      * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))           &
                      * ( qcfekp1(i)-qcfek(i) )                                &
    !                  Smoothed forced detrainment
                      + a_smth*tmp_fd_dqcfek                                   &
    !                  Mixing detrainment
                      + amdetk(i)                                              &
                      * ( qcfpk(i) - qcfek(i) ) )

    if (l_wtrac_conv) qcfek_temp(i) = qcfek(i) + tmp_dqcfek*timestep

    ! Temporary storage of condensate increments on level k+1.
    ! These will be used for doing smooth phase transition
    ! (and associated correction to temperatures).
    tmp_dqclekp1      = flxbydpk(i) * (1.0-a_smth)                             &
                      * delpk(i)/delpkp1(i)*tmp_fd_dqclek

    tmp_dqcfekp1      = flxbydpk(i) * (1.0-a_smth)                             &
                      * delpk(i)/delpkp1(i)*tmp_fd_dqcfek

    if (tmp_dqclekp1 > 0.0) then
      ! Add some of the liquid increment to the ice.
      dqcfekp1(i) = frac_icekp1 * tmp_dqclekp1
      ! Adjust for latent heating from this freezing.
      dthekp1(i) = dthekp1(i) + frac_icekp1 * tmp_dqclekp1 * lf/(exkp1(i)*cp)
      ! Add the rest of the liquid increment to the liquid.
      dqclekp1(i) = frac_liqkp1 * tmp_dqclekp1
      ! Store freezing for water tracer use
      if (l_wtrac_conv) dqclekp1_frz(i) = frac_icekp1*tmp_dqclekp1

    else
      dqclekp1(i) = tmp_dqclekp1
      dqcfekp1(i) = 0.0
    end if

    if (l_fix_pc2_cnv_mix_phase) then
      !Corrected version that conserves water
      if (tmp_dqcfekp1 > 0.0) then
        ! Add some of the ice increment to the liquid.
        dqclekp1(i) = dqclekp1(i) + frac_liqkp1 * tmp_dqcfekp1
        ! Adjust for latent cooling from this melting.
        dthekp1(i) = dthekp1(i) - frac_liqkp1 * tmp_dqcfekp1 * lf/(exkp1(i)*cp)
        ! Add the rest of the ice increment to the ice.
        dqcfekp1(i) = dqcfekp1(i) + frac_icekp1 * tmp_dqcfekp1
        ! Store melt for water tracer use
        if (l_wtrac_conv) dqcfekp1_mlt(i) = frac_liqkp1*tmp_dqcfekp1
      else
        dqcfekp1(i) = dqcfekp1(i) + tmp_dqcfekp1
      end if
    else
      !Original version that does not always conserve water
      if (tmp_dqcfekp1 > 0.0) then
        ! Add some of the ice increment to the liquid.
        dqclekp1(i) = frac_liqkp1 * tmp_dqcfekp1
        ! Adjust for latent cooling from this melting.
        dthekp1(i) = dthekp1(i) - frac_liqkp1 * tmp_dqcfekp1 * lf/(exkp1(i)*cp)
        ! Add the rest of the ice increment to the ice.
        dqcfekp1(i) = frac_icekp1 * tmp_dqcfekp1
        if (l_wtrac_conv) dqcfekp1_mlt(i) = frac_liqkp1*tmp_dqcfekp1
      else
        dqcfekp1(i) = tmp_dqcfekp1
      end if
    end if

    ! ----------------------------------------------------------------------
    !       Adjust temperature increment and condensate increments to
    !       take account of a smoother transistion between water and ice
    ! ----------------------------------------------------------------------

    if (i_pc2_conv_coupling == pc2_conv_original    .or.                       &
        i_pc2_conv_coupling == pc2_conv_maxincldqc) then
      ! Original method
      dqclek(i) = dqclek(i) + tmp_dqclek
      dqcfek(i) = dqcfek(i) + tmp_dqcfek
    else
      ! Convective plume is either liq or ice with an abrupt change.
      ! Partition the condensate increments detrained from the
      ! convective plume onto the large-scale more smoothly over
      ! a range of temperature. Adjust for latent heating.
      frac_icek   = (thek(i)*exk(i)-n_starticetK)/                             &
             (n_allicedegC-(n_starticetK-tm))
      frac_icek   = min(max(0.0,frac_icek),1.0)
      frac_liqk   = 1.0-frac_icek

      if (tmp_dqclek > 0.0) then
        dqcfek(i) = dqcfek(i) + frac_icek*tmp_dqclek
        ! Extra freezing heats
        dthek(i)  = dthek(i)  + frac_icek*tmp_dqclek*lfrcp/exk(i)
        dqclek(i) = dqclek(i) + frac_liqk*tmp_dqclek
        ! Store freezing for water tracer use
        if (l_wtrac_conv) dqclek_frz(i) = frac_icek*tmp_dqclek
      else
        dqclek(i) = dqclek(i) + tmp_dqclek
      end if

      if (tmp_dqcfek > 0.0) then
        dqclek(i) = dqclek(i) + frac_liqk*tmp_dqcfek
        ! Melting cools
        dthek(i)  = dthek(i)  - frac_liqk*tmp_dqcfek*lfrcp/exk(i)
        dqcfek(i) = dqcfek(i) + frac_icek*tmp_dqcfek
        ! Store melt for water tracer use
        if (l_wtrac_conv) dqcfek_mlt(i) = frac_liqk*tmp_dqcfek
      else
        dqcfek(i) = dqcfek(i) + tmp_dqcfek
      end if

    end if

  end do
else
  !-----------------------------------------------------------------------
  !         PC2 is off. Set the qcl and qcf increments to zero.
  !-----------------------------------------------------------------------
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    dqclek(i)     = 0.0
    dqclekp1(i)   = 0.0

    dqcfek(i)     = 0.0
    dqcfekp1(i)   = 0.0
  end do
end if  !l_q_interact

! ---------------------------------------------------------------------
!  Calculate effect of convection upon momentum of layer k and
!  do terminal detrainment of momentum.
!
!  Rate of change of wind field by convection is estimated using a
!  divergence of vertical eddy momentum flux across the layer.
! --------------------------------------------------------------------
!
! All convective momentum transport calculations for the cumulus convection
! (deep and shallow) done when convection terminates.

if (l_mom_gk) then      ! Gregory-Kershaw CMT

  !-------------------------------------
  ! Set up and vp levels to be referred
  !------------------------------------
  allocate(up_ref(npnts))
  allocate(vp_ref(npnts))

  if (l_mom_gk_stable) then ! more stable version
!DIR$ IVDEP
    do m=1, ni
      i = idx(m)
      up_ref(i) = upkp1(i)
      vp_ref(i) = vpkp1(i)
    end do
  else                     ! original discretisation
!DIR$ IVDEP
    do m=1, ni
      i = idx(m)
      up_ref(i) = upk(i)
      vp_ref(i) = vpk(i)
    end do
  end if

!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    !----------------------------------------------------------------------
    ! Estimate eddy flux at top of current UV layer due to convection
    !----------------------------------------------------------------------

    flx_u_kp0p5 = flxk(i) * (1.0-amdetk(i)) * (1.0-deltak(i)) *                &
                            (1.0+ekp14(i)) * (up_ref(i)-uekp1(i))
    flx_v_kp0p5 = flxk(i) * (1.0-amdetk(i)) * (1.0-deltak(i)) *                &
                            (1.0+ekp14(i)) * (vp_ref(i)-vekp1(i))

    if (blowst(i)) then
      !----------------------------------------------------------------------
      ! Initial convecting layer - no flux at base of layer
      !----------------------------------------------------------------------
      duek(i) = duek(i) - flx_u_kp0p5 / delp_uv_k(i)
      dvek(i) = dvek(i) - flx_v_kp0p5 / delp_uv_k(i)
      !----------------------------------------------------------------------
      ! Store eddy flux at top of current UV layer ready for calculation
      ! of next layer.
      !----------------------------------------------------------------------
      eflux_u_ud(i) = flx_u_kp0p5
      eflux_v_ud(i) = flx_v_kp0p5

    else
      !----------------------------------------------------------------------
      ! Convecting layer - take eddy flux divergence across the layer
      !----------------------------------------------------------------------
      duek(i) = duek(i) - ( (flx_u_kp0p5 - eflux_u_ud(i)) /delp_uv_k(i) )

      dvek(i) = dvek(i) - ( (flx_v_kp0p5 - eflux_v_ud(i)) /delp_uv_k(i) )

      !----------------------------------------------------------------------
      ! Store eddy flux at top of curent UV layer ready for calculation of
      ! next layer
      !----------------------------------------------------------------------
      eflux_u_ud(i) = flx_u_kp0p5
      eflux_v_ud(i) = flx_v_kp0p5

    end if

    if (bterm(i)) then
      !----------------------------------------------------------------------
      ! Convection terminates - calculate increment due to convection in top
      ! layer - no flux out of top layer.
      !----------------------------------------------------------------------
      duekp1(i)  = eflux_u_ud(i) / delp_uv_kp1(i)
      dvekp1(i)  = eflux_v_ud(i) / delp_uv_kp1(i)
      !----------------------------------------------------------------------
      ! Zero eddy flux out of top layer.
      !----------------------------------------------------------------------
      eflux_u_ud(i) = 0.0
      eflux_v_ud(i) = 0.0
    else
      duekp1(i)  = 0.0
      dvekp1(i)  = 0.0
    end if

  end do

  deallocate(vp_ref)
  deallocate(up_ref)

end if   ! l_mom_gk

!----------------------------------------------------------------------
!  Effect of convection on tracer content of layer k.
!  (looping over number of tracer variables)
!  and do terminal detrainment of tracer.
!----------------------------------------------------------------------

if (l_tracer) then
  do ktra = 1,ntra
!DIR$ IVDEP
    do m=1, ni
      i = idx(m)

      !                   Smoothed forced detrainment term
      tmp_fd_dtraek   =  deltak(i) * (1.0-amdetk(i))                           &
                      * ( trap(i,k,ktra)-trae(i,k,ktra) )

      !                   tracer increment at level k
      dtrae(i,k,ktra)  = dtrae(i,k,ktra) + flxbydpk(i)                         &
      !                   Compensating subsidence
                            * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))     &
                            * (trae(i,k+1,ktra)-trae(i,k,ktra))                &
      !                   Smoothed forced detrainment
                            + a_smth*tmp_fd_dtraek                             &
      !                   Mixing detrainment
                            + amdetk(i)                                        &
                            * ( trap(i,k,ktra)-trae(i,k,ktra) ) )

      !                   tracer increment at level k+1 due
      !                   to smoothed  forced detrainment
      dtrae(i,k+1,ktra)= flxbydpk(i) *(1.0-a_smth)                             &
                      * delpk(i)/delpkp1(i)*tmp_fd_dtraek

    end do  ! i
  end do    ! ktra
end if      ! l_tracer

! Water tracers (assuming use of PC2 here)
if (l_wtrac_conv) then    ! Water tracers

  call environ_wtrac(k, npnts, ni, n_wtrac,                                    &
                     idx, timestep, a_smth, deltak,                            &
                     amdetk, flxbydpk, ekp14, delpk, delpkp1,                  &
                     qclek_temp, qcfek_temp,                                   &
                     dqclek_frz, dqcfek_mlt, dqclekp1_frz, dqcfekp1_mlt,       &
                     qrk_wtrac, wtrac_p, wtrac_e)

  deallocate(qclek_temp)
  deallocate(qcfek_temp)
  deallocate(dqcfekp1_mlt)
  deallocate(dqclekp1_frz)
  deallocate(dqcfek_mlt)
  deallocate(dqclek_frz)

end if ! l_wtrac_conv

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine environ_6a
end module environ_6a_mod
