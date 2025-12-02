! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Completes lifting of the parcel from layer k to k+1

module convec2_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Completes lifting of the parcel from layer k to k+1.
!
!   Calls subroutines parcel and environ
!
!   Subroutine parcel carries out the forced
!   detrainment calculation, tests to see if convection is termintating
!   and calculates the precipitation rate from layer k+1.
!
!   Subroutine environ calculates the effect of convection upon the
!   large-scale atmosphere.
!
!   The parcel properties at level k+1 are then fully updated. CAPE and
!   dCAPEbydt are calculated for use in the CAPE closure.
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


character(len=*), parameter, private :: ModuleName='CONVEC2_6A_MOD'

contains

subroutine convec2_6a(k, npnts, np_full, nlev, ntra, n_wtrac, trlev,           &
                      ad_on, new_termc, start_lev, timestep,                   &
                      pk, pkp1, delpk,                                         &
                      delpkp1, delp_uv_k, delp_uv_kp1,                         &
                      exk, exkp1,                                              &
                      thek, thekp1, qek, qekp1,                                &
                      qclek, qclekp1, qcfek, qcfekp1,                          &
                      qsek, qsekp1,                                            &
                      cflek, cflekp1,  cffek,  cffekp1,                        &
                      thpk, qpk, qclpk, qcfpk,                                 &
                      rbuoyk, rbuoykp1, rbuoyukp1,                             &
                      watldek, watldekp1, watldpk, watldpkp1,                  &
                      Qlkp1, Qfkp1, Frezkp1,                                   &
                      ekp14, ekp34, amdetk, flxk, flx_init,                    &
                      uek, uekp1, vek, vekp1,                                  &
                      upk, vpk,                                                &
                      trae,                                                    &
                      zk, zkp12, zkp1,                                         &
                      l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,       &
                      bgmk, bgmkp1, bwk,                                       &
                      bwkp1, blowst, bland,                                    &

                      ! In/out
                      lcbase, lctop,                                           &
                      thpkp1, qpkp1, qclpkp1, qcfpkp1,                         &
                      dthek, dqek, dqclek, dqcfek,                             &
                      tcw, depth, cclwp, lcca,                                 &
                      cape, fcape, dcpbydt,                                    &
                      relh, dptot, deltaktot, max_cfl,                         &
                      eflux_u_ud, eflux_v_ud,                                  &
                      duek, dvek,                                              &
                      dtrae, trap, wtrac_e, wtrac_p,                           &
                      w2p_k, bterm, blatent, xsbmin,                           &

                      ! Out
                      iccb, icct,                                              &
                      dcflek, dcffek, dbcfek,                                  &
                      dthekp1, dqekp1, dqclekp1, dqcfekp1,                     &
                      dcflekp1, dcffekp1, dbcfekp1,                            &
                      prekp1, thrk, qrk, deltak,                               &
                      flxkp12, flxkp1,                                         &
                      cca, ccwkp1,                                             &
                      upkp1, vpkp1,                                            &
                      duekp1, dvekp1,                                          &
                      w2p_kp1,tnuc_nlcl,                                       &

                      !Indirect indexing
                      idx,ni )

use cv_run_mod,           only: cnv_wat_load_opt
use cv_stash_flg_mod,     only: flg_w_eqn
use planet_constants_mod, only: c_virtual, r

use ukca_option_mod,     only: l_ukca, l_ukca_plume_scav
use ukca_scavenging_mod, only: ukca_plume_scav

use yomhook,         only: lhook, dr_hook
use parkind1,        only: jprb, jpim
use parcel_6a_mod,   only: parcel_6a
use calc_w_eqn_mod,  only: calc_w_eqn
use environ_6a_mod,  only: environ_6a
use pc2_environ_mod, only: pc2_environ
use wtrac_conv_mod,  only: l_wtrac_conv, conv_e_wtrac_type, conv_p_wtrac_type

implicit none

! Subroutine arguments

integer,intent(in) :: k             ! present model layer
integer,intent(in) :: npnts         ! Number of points
integer,intent(in) :: np_full       ! Full vector length
integer,intent(in) :: nlev          ! Number of model levels for calculations
integer,intent(in) :: ntra          ! Number of tracer variables
integer,intent(in) :: n_wtrac       ! Number of water tracer variables
integer,intent(in) :: trlev         ! No. of model levels on which
                                    ! tracers are included
integer,intent(in) :: ad_on         ! Flag for adaptive detrainment
integer,intent(in) :: new_termc     ! Method of terminating the convective
                                    ! parcel ascent

integer,intent(in) :: start_lev(npnts) ! Level at which convection is initiated

real(kind=real_umphys),intent(in) :: timestep      ! Timestep
real(kind=real_umphys),intent(in) :: pk(npnts)
                                    ! pressure at mid-point of layer k (Pa)
real(kind=real_umphys),intent(in) :: pkp1(npnts)
                                    ! pressure at mid-point of layer k+1 (Pa)
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
real(kind=real_umphys),intent(in) :: qsek(npnts)
                                    ! Env. saturated specific humidity in
                                    ! in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qsekp1(npnts)
                                    ! Env. saturated specific humidity in
                                    ! in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: cflek(npnts)
                                    ! Env. liquid cloud volume fraction
                                    ! in layer k
real(kind=real_umphys),intent(in) :: cflekp1(npnts)
                                    ! Env. liquid cloud volume fraction
                                    ! in layer k+1
real(kind=real_umphys),intent(in) :: cffek(npnts)
                                    ! Env. frozen cloud volume fraction
                                    ! in layer k
real(kind=real_umphys),intent(in) :: cffekp1(npnts)
                                    ! Env. frozen cloud volume fraction
                                    ! in layer k+1
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
real(kind=real_umphys),intent(in) :: watldpk(npnts)
                                    ! Par. water loading in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: watldekp1(npnts)
                                    ! Env. water loading in layer k+1 (kg/kg)
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

! Height above surface of model level ...
real(kind=real_umphys),intent(in) :: zk(npnts)    ! ...k     [m]
real(kind=real_umphys),intent(in) :: zkp12(npnts) ! ...k+1/2 [m]
real(kind=real_umphys),intent(in) :: zkp1(npnts)  ! ...k+1   [m]

logical,intent(in) :: l_q_interact  ! True if PC2 is switched on
logical,intent(in) :: l_mom_gk      ! Switch for inclusion of Gregory-Kershaw
                                    ! CMT
logical,intent(in) :: l_mom_gk_stable ! Switch for stabilised Gregory-Kershaw
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
logical,intent(in) :: blowst(npnts) ! mask for those points at which stability
                                    ! is low enough for convection to occur
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
real(kind=real_umphys),intent(in out) :: dcflek(npnts)
                                     ! Increment to liquid cloud volume fraction
                                    ! in layer k
real(kind=real_umphys),intent(in out) :: dcffek(npnts)
                                     ! Increment to frozen cloud volume fraction
                                    ! in layer k
real(kind=real_umphys),intent(in out) :: dbcfek(npnts)
                                     ! Increment to total cloud volume fraction
                                    ! in layer k
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
real(kind=real_umphys),intent(in out) :: cape(npnts)
                                     ! Convective available potential energy
                                    ! in up to layer k-1
                                    ! out up to layer k
real(kind=real_umphys),intent(in out) :: fcape(npnts)
                                     ! CAPE weighted by f.det profile
                                    ! in up to layer k-1
                                    ! out up to layer k
real(kind=real_umphys),intent(in out) :: dcpbydt(npnts)! Rate of change of CAPE
                                    ! in up to layer k-1
                                    ! out up to layer k
real(kind=real_umphys),intent(in out) :: relh(npnts)
                                     ! Relative humidity integral (averaged
                                    ! when convection terminates)
real(kind=real_umphys),intent(in out) :: dptot(npnts)  ! Delta P integral
real(kind=real_umphys),intent(in out) :: deltaktot(npnts)
                                         ! Integrated forced detrainment
real(kind=real_umphys),intent(in out) :: max_cfl(npnts)! CFL ratio
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
real(kind=real_umphys),intent(in out) :: dtrae(np_full,nlev,ntra)
                                                ! Increment to tracer
                                                ! (kg/kg/s)

real(kind=real_umphys),intent(in out) :: trap(np_full,nlev,ntra)
                                                ! Par. tracer content
                                                ! (kg/kg)

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                                ! Structure containing env
                                                ! water tracer fields

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)
                                                ! Structure containing parcel
                                                ! water tracer fields

real(kind=real_umphys), intent(in out) :: w2p_k(npnts)
                                        ! (Parcel vertical velocity)^2
                                        ! on level k, [(m/s)^2]

logical,intent(in out) :: bterm(npnts)  ! Mask for parcels which terminate
                                        ! in layer k+1
logical,intent(in out) :: blatent(npnts)! Mask for points where latent heat has
                                        ! been released
real(kind=real_umphys),intent(in out) :: xsbmin(npnts)
                                        ! Threshold buoyancy for forced
                                        ! detrainment (K)

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
integer,intent(in out) :: iccb(npnts)  ! convective cloud base_level
integer,intent(in out) :: icct(npnts)  ! convective cloud top level
! Need to be inout as need to be remembered between subsequent calls to convec2

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
real(kind=real_umphys),intent(in out) :: dcflekp1(npnts)
                                       ! Increment to liquid cloud volume
                                       ! fraction in layer k+1
real(kind=real_umphys),intent(in out) :: dcffekp1(npnts)
                                       ! Increment to frozen cloud volume
                                       ! fraction in layer k+1
real(kind=real_umphys),intent(in out) :: dbcfekp1(npnts)
                                       ! Increment to total cloud volume
                                       ! fraction in layer k+1
real(kind=real_umphys),intent(in out) :: prekp1(npnts)
                                       ! precipitation from parcel as it rises
                                       ! from layer k to k+1 (kg/m**2/s)
real(kind=real_umphys),intent(in out) :: thrk(npnts)
                                       ! pot. temperature of forced detrained
                                       ! parcel in layer k (K)
real(kind=real_umphys),intent(in out) :: qrk(npnts)
                                       ! Specific humidity of forced detrained
                                       ! parcel in layer k (kg/kg)
real(kind=real_umphys),intent(in out) :: deltak(npnts)
                                       ! Parcel forced detrainment rate in
                                       ! layer k multiplied by layer thickness
real(kind=real_umphys),intent(in out) :: flxkp12(npnts)
                                       ! parcel massflux in layer k+1/2 (Pa/s)
real(kind=real_umphys),intent(in out) :: flxkp1(npnts)
                                       ! parcel massflux in layer k+1 (Pa/s)
real(kind=real_umphys),intent(in out) :: cca(npnts)
                                       ! convective cloud amount (%)
real(kind=real_umphys),intent(in out) :: ccwkp1(npnts)
                                       ! Total condensate in level k+1 (kg/kg)
real(kind=real_umphys),intent(in out) :: upkp1(npnts)
                                       ! Par. U in layer k+1 (m/s)
real(kind=real_umphys),intent(in out) :: vpkp1(npnts)
                                       ! Par. V in layer k+1 (m/s)
real(kind=real_umphys),intent(in out) :: duekp1(npnts)
                                       ! Increment to U in layer k+1 (m/s**2)
real(kind=real_umphys),intent(in out) :: dvekp1(npnts)
                                       ! Increment to V in layer k+1 (m/s**2)
real(kind=real_umphys),intent(in out) :: w2p_kp1(npnts)
                                       ! (Parcel vertical velocity)^2
                                       ! on level k+1, [(m/s)^2]
real(kind=real_umphys),intent(in) :: tnuc_nlcl(npnts)
                                       ! nucleation temperature as function of
                                       ! dust indexed using nlcl(deg cel)
! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni


!-------------------------------------------------------------------------------
! Local variables

integer :: i, m, ktra     ! loop counter

real(kind=real_umphys) :: thvp              ! Virtual tempature of parcel
real(kind=real_umphys) :: thve              ! Virtual tempature of environment
real(kind=real_umphys) :: rho
                          ! Density required in CAPE calculations
real(kind=real_umphys) :: tmp_dcpbydt       ! Temporary dcpbydt

real(kind=real_umphys), allocatable :: qrk_wtrac(:,:)
                                    ! Water tracer vapour content of forced
                                    ! detrained parcel in layer k (kg/kg)

! Compressed variable
real(kind=real_umphys) :: trapkp1_c(ni,ntra)
real(kind=real_umphys) :: qclpkp1_c(ni)
real(kind=real_umphys) :: prekp1_c(ni)
real(kind=real_umphys) :: flxkp1_c(ni)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CONVEC2_6A'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  Complete lifting parcel to layer k+1.
!
!  subroutine parcel
!
!  UM Documentation paper 27
!  Sections (5),(6),(7),(8),(9)
! ----------------------------------------------------------------------

! Allocate water tracer array to hold equivalent of qrk
if (l_wtrac_conv) then
  allocate(qrk_wtrac(npnts,n_wtrac))
else
  allocate(qrk_wtrac(1,1))
end if

call parcel_6a (k, npnts, np_full, nlev, ntra, trlev, n_wtrac,                 &
                ad_on, new_termc,                                              &
                start_lev,                                                     &
                pk, pkp1, delpkp1, exk, exkp1,                                 &
                thek, thekp1, qek, qekp1,                                      &
                qclek, qclekp1, qcfek, qcfekp1,                                &
                qsekp1,                                                        &
                thpk, qpk, qclpk, qcfpk,                                       &
                rbuoyk, rbuoykp1, rbuoyukp1,                                   &
                watldek, watldekp1, watldpk, watldpkp1,                        &
                Qlkp1, Qfkp1, Frezkp1,                                         &
                ekp14, ekp34, amdetk, flxk, flx_init,                          &
                uek, uekp1, vek, vekp1,                                        &
                upk, vpk,                                                      &
                trae, wtrac_e,                                                 &
                l_q_interact, l_mom_gk, l_tracer,                              &
                bgmk, bgmkp1, bwk, bwkp1, bland,                               &
                ! In/out
                lcbase, lctop,                                                 &
                thpkp1, qpkp1, qclpkp1, qcfpkp1,                               &
                tcw, depth, cclwp, lcca,                                       &
                bterm, blatent, xsbmin,                                        &
                trap, wtrac_p,                                                 &
                ! Out
                iccb, icct, prekp1,                                            &
                thrk, qrk, qrk_wtrac, deltak, flxkp1, flxkp12,                 &
                cca, ccwkp1,                                                   &
                upkp1, vpkp1,                                                  &
                ! Indirect indexing
                idx, ni)

! ----------------------------------------------------------------------
!  Calculate the change in parcel tracer content due to scavenging
!  by precipitation for UKCA tracers.
!  UM Documentation paper 84
! ----------------------------------------------------------------------

if (l_tracer .and. l_ukca .and. l_ukca_plume_scav .and. npnts > 0              &
    .and. any(prekp1 > 1e-30)) then

  ! Compression

!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    qclpkp1_c(m) = qclpkp1(i)
    prekp1_c(m)  = prekp1(i)
    flxkp1_c(m)  = flxkp1(i)

    do ktra = 1,ntra
      trapkp1_c(m,ktra) = trap(i,k+1,ktra)
    end do

  end do

  call ukca_plume_scav(ntra, ni, trapkp1_c,                                    &
                       qclpkp1_c, prekp1_c, flxkp1_c)

  ! Decompression

!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    do ktra = 1,ntra
      trap(i,k+1,ktra) = trapkp1_c(m,ktra)
    end do
  end do

end if

! ----------------------------------------------------------------------
!  Calculate the increments to environmental theta, q, qcl, qcf,
!  winds and tracers.
! ----------------------------------------------------------------------

call environ_6a (k, npnts, np_full, nlev, ntra, trlev, n_wtrac,                &
              timestep, delpk, delpkp1, delp_uv_k, delp_uv_kp1,                &
              exk, exkp1,                                                      &
              thek, thekp1, qek, qekp1,                                        &
              qclek, qclekp1, qcfek, qcfekp1,                                  &
              thpk,  qpk,                                                      &
              qclpk, qcfpk,                                                    &
              thrk, qrk, qrk_wtrac,                                            &
              ekp14, amdetk, deltak, flxk,                                     &
              uekp1, vekp1,                                                    &
              upk, upkp1, vpk, vpkp1,                                          &
              trae,                                                            &
              trap,                                                            &
              l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,               &
              blowst, bterm, wtrac_p,                                          &
              ! In/out
              dthek, dqek, dqclek, dqcfek,                                     &
              eflux_u_ud, eflux_v_ud,                                          &
              duek, dvek,                                                      &
              dtrae, wtrac_e,                                                  &
              ! Out
              dthekp1, dqekp1, dqclekp1, dqcfekp1,                             &
              duekp1, dvekp1,tnuc_nlcl,                                        &
              ! Indirect indexing
              idx, ni)

deallocate(qrk_wtrac)

! ----------------------------------------------------------------------
!  Calculate the increments to environmental cloud fractions
! ----------------------------------------------------------------------
call pc2_environ (   npnts,                                                    &
              qclek, qclekp1, qcfek, qcfekp1,                                  &
              cflek, cflekp1, cffek, cffekp1,                                  &
              qclpk, qcfpk,                                                    &
              dqclek, dqcfek,                                                  &
              dqclekp1, dqcfekp1,                                              &
              l_q_interact,                                                    &
              bterm,                                                           &
              ! Out
              dcflek, dcffek, dbcfek,                                          &
              dcflekp1, dcffekp1, dbcfekp1,                                    &
              ! Indirect indexing
              idx, ni)


! ---------------------------------------------------------------------
!  Calculate contribution to CAPE, the rate of change of CAPE  due to
!  the updraught and Courant number
! ---------------------------------------------------------------------

!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  thvp    = thpk(i)*(1.0+c_virtual*qpk(i)-watldpk(i))
  thve    = thek(i)*(1.0+c_virtual*qek(i)-watldek(i))
  rho     = pk(i)/(r*thve*exk(i))

  cape(i) = cape(i)+(thvp-thve)*delpk(i)/(rho*thve)
  relh(i) = relh(i)+(qek(i)/qsek(i))*delpk(i)
  dptot(i)= dptot(i)+delpk(i)

  fcape(i)= fcape(i)+flxk(i)*deltak(i)*cape(i)
  deltaktot(i)= deltaktot(i)+flxk(i)*deltak(i)

  if (cnv_wat_load_opt == 0) then
    !   Water loading is not included and therefore do not include
    !   the contribution from the condensate increments
    tmp_dcpbydt = ( dthek(i)*(1.0+c_virtual*qek(i))                            &
                  + thek(i)*c_virtual*dqek(i) )                                &
                  * (delpk(i)/(rho*thve))
  else
    !   Water loading is on and therefore include the contribution
    !   from the condensate increments. If PC2 is off the condensate
    !   increments will be zero.
    tmp_dcpbydt = ( dthek(i)*(1.0+c_virtual*qek(i)-watldek(i))                 &
                  + thek(i)*(c_virtual*dqek(i)                                 &
                  - dqclek(i) - dqcfek(i)) )                                   &
                  * (delpk(i)/(rho*thve))
  end if

  if (tmp_dcpbydt  >   0.0) then
    dcpbydt(i) = dcpbydt(i) + tmp_dcpbydt
  end if

  ! Calculate courant number using mass flux at half level
  max_cfl(i)=max(max_cfl(i),                                                   &
             flxk(i)/delpk(i)*(1+ekp14(i))*(1.0-deltak(i))*(1.0-amdetk(i)))

end do


if (flg_w_eqn) then
  call calc_w_eqn ( npnts, bterm, blowst                                       &
                  , ekp14, ekp34, zk, zkp12, zkp1                              &
                  , thek,  thekp1,  thpk,  thpkp1                              &
                  , qek,   qekp1,   qpk,   qpkp1                               &
                  , qclek, qclekp1, qclpk, qclpkp1                             &
                  , qcfek, qcfekp1, qcfpk, qcfpkp1                             &
                  , w2p_k, w2p_kp1, idx, ni )
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine convec2_6a
end module convec2_6a_mod
