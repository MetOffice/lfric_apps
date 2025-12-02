! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Initialise many arrays for glue_conv_6a

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

module init_conv6a_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'INIT_CONV6A_MOD'
contains

subroutine init_6a_rain(npnts, nlev, n_sh, n_dp, n_cg, n_md, rain_3d_sh,       &
                        snow_3d_sh, rain_3d_dp, snow_3d_dp, rain_3d_cg,        &
                        snow_3d_cg, rain_3d_md, snow_3d_md, ind_cape_reduced,  &
                        cape_ts_used, ind_deep, ind_shall, kterm_deep,         &
                        kterm_shall, kterm_congest )

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts,               & ! Full number of points
  nlev,                & ! Number of levels
  n_sh,                & ! Number of points with shallow convection
  n_dp,                & ! Number of points with deep convection
  n_cg,                & ! Number of points with congestus convection
  n_md                   ! Number of points with mid-level convection


real(kind=real_umphys), intent(out) ::                                         &
  rain_3d_sh(n_sh,nlev),  & ! Rain from shallow convection
  snow_3d_sh(n_sh,nlev),  & ! Snow from shallow convection
  rain_3d_dp(n_dp,nlev),  & ! Rain from shallow convection
  snow_3d_dp(n_dp,nlev),  & ! Snow from shallow convection
  rain_3d_cg(n_cg,nlev),  & ! Rain from congestus convection
  snow_3d_cg(n_cg,nlev),  & ! Snow from congestus convection
  rain_3d_md(n_md,nlev),  & ! Rain from mid-level convection
  snow_3d_md(n_md,nlev),  & ! Snow from mid-level convection
  ind_cape_reduced(npnts),& ! indicates cape timescale reduced
  cape_ts_used(npnts),    & ! cape timescale for deep (s)
  ind_deep(npnts),        & ! indicator of real deep
  ind_shall(npnts)          ! indicator of real shallow

integer, intent(out) :: kterm_deep(npnts) ! index deep conv
integer, intent(out) :: kterm_shall(npnts) ! level for shallow termination
integer, intent(out) :: kterm_congest(npnts) ! termination level for congestus

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i, k             ! loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='INIT_6A_RAIN'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise 3d rain variables
!
do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_sh
    rain_3d_sh(i,k) = 0.0
    snow_3d_sh(i,k) = 0.0
  end do
end do

do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_dp
    rain_3d_dp(i,k) = 0.0
    snow_3d_dp(i,k) = 0.0
  end do
end do

do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_cg
    rain_3d_cg(i,k) = 0.0
    snow_3d_cg(i,k) = 0.0
  end do
end do

do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_md
    rain_3d_md(i,k) = 0.0
    snow_3d_md(i,k) = 0.0
  end do
end do

do i = 1,npnts
  ind_cape_reduced(i) = 0.0
  cape_ts_used(i)     = 0.0
  ind_deep(i)         = 0.0
  ind_shall(i)        = 0.0
  kterm_deep(i)       = 0
  kterm_shall(i)      = 0
  kterm_congest(i)    = 0
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine init_6a_rain

subroutine init_6a_w2p(nlev, n_sh, n_dp, n_cg, n_md, w2p_sh, w2p_cg, w2p_md,   &
                     w2p_dp )

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! Subroutine arguments

integer, intent(in) ::                                                         &
  nlev,                & ! Number of levels
  n_sh,                & ! Number of points with shallow convection
  n_dp,                & ! Number of points with deep convection
  n_cg,                & ! Number of points with congestus convection
  n_md                   ! Number of points with mid-level convection

! Arrays for (Parcel vertical velocity)^2 [(m/s)^2] on ...
real(kind=real_umphys), intent(out) :: w2p_sh(n_sh,nlev)
                                       ! ...shallow   convection points
real(kind=real_umphys), intent(out) :: w2p_cg(n_cg,nlev)
                                       ! ...congestus convection points
real(kind=real_umphys), intent(out) :: w2p_md(n_md,nlev)
                                       ! ...mid-level convection points
real(kind=real_umphys), intent(out) :: w2p_dp(n_dp,nlev)
                                       ! ...deep      convection points

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i, k             ! loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='INIT_6A_W2P'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise local arrays for w-eqn
!
do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_sh
    w2p_sh(i,k) = 0.0
  end do
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_md
    w2p_md(i,k) = 0.0
  end do
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_dp
    w2p_dp(i,k) = 0.0
  end do
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_cg
    w2p_cg(i,k) = 0.0
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine init_6a_w2p

subroutine init_6a_npnts(npnts, nlev, n_cca_lev, np_field, n_wtrac,            &
                         dqbydt, dthbydt, dqclbydt, dqcfbydt, dcflbydt,        &
                         dcffbydt, dbcfbydt, dqbydt_wtrac, dqclbydt_wtrac,     &
                         dqcfbydt_wtrac, ccw, cca,                             &
                         dt_dd, dq_dd, du_dd, dv_dd, area_ud, area_dd,         &
                         up_flux,up_flux_half,dwn_flux,                        &
                         entrain_up,detrain_up,   entrain_dwn, detrain_dwn,    &
                         mf_deep, mf_congest,                                  &
                         mf_shall, mf_midlev, dt_deep, dt_congest, dt_shall,   &
                         dt_midlev, dq_deep, dq_congest, dq_shall, dq_midlev,  &
                         du_deep, du_congest, du_shall, du_midlev, dv_deep,    &
                         dv_congest, dv_shall, dv_midlev, wqt_flux_sh,         &
                         wql_flux_sh, wthetal_flux_sh, wthetav_flux_sh,        &
                         uw_deep, vw_deep, uw_shall, vw_shall, uw_mid, vw_mid )

use cv_stash_flg_mod, only:                                                    &
   flg_up_flx, flg_up_flx_half, flg_dwn_flx,                                   &
   flg_entr_up, flg_entr_dwn, flg_detr_up, flg_detr_dwn,                       &
   flg_uw_dp, flg_vw_dp, flg_uw_shall, flg_vw_shall, flg_uw_mid, flg_vw_mid,   &
   flg_wqt_flux, flg_wthetal_flux, flg_wthetav_flux, flg_wql_flux,             &
   flg_mf_deep, flg_mf_congest, flg_mf_shall, flg_mf_midlev,                   &
   flg_dt_deep, flg_dt_congest, flg_dt_shall, flg_dt_midlev,                   &
   flg_dq_deep, flg_dq_congest, flg_dq_shall, flg_dq_midlev,                   &
   flg_du_deep, flg_du_congest, flg_du_shall, flg_du_midlev,                   &
   flg_dv_deep, flg_dv_congest, flg_dv_shall, flg_dv_midlev

use wtrac_conv_mod, only: l_wtrac_conv

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! Subroutine arguments

integer, intent(in) ::                                                         &
  nlev,                & ! Number of levels
  n_cca_lev,           & ! Number of conv cloud amount levels
  npnts,               & ! Number of points
  np_field,            & ! Full size of fields
  n_wtrac                ! Number of water tracers

!Increments due to convection
real(kind=real_umphys), intent(out) :: dqbydt(np_field,nlev) ! Increments to q
real(kind=real_umphys), intent(out) :: dthbydt(np_field,nlev)
                                            ! Increments to potl temp
real(kind=real_umphys), intent(out) :: dqclbydt(np_field,nlev)
                                             ! Increments to liq condensate
real(kind=real_umphys), intent(out) :: dqcfbydt(np_field,nlev)
                                             ! Increments to ice condensate
real(kind=real_umphys), intent(out) :: dcflbydt(np_field,nlev)
                                             ! Increments to liq cloud volume
real(kind=real_umphys), intent(out) :: dcffbydt(np_field,nlev)
                                             ! Increments to ice cloud volume
real(kind=real_umphys), intent(out) :: dbcfbydt(np_field,nlev)
                                             ! Increments to total cld volume
real(kind=real_umphys), intent(out) :: dqbydt_wtrac(np_field,nlev,n_wtrac)
                                             ! Increments to
                                             ! water tracer vapour
real(kind=real_umphys), intent(out) :: dqclbydt_wtrac(np_field,nlev,n_wtrac)
                                             ! Increments to
                                             ! water tracer liq condensate
real(kind=real_umphys), intent(out) :: dqcfbydt_wtrac(np_field,nlev,n_wtrac)
                                             ! Increments to
                                             ! water tracer ice condensate
real(kind=real_umphys), intent(out) :: cca(np_field,n_cca_lev)
                                             ! Conv. cld amount (0-1) ??? nlev
real(kind=real_umphys), intent(out) :: ccw(np_field,nlev)
                                             ! Cnv. in-cld liquid water (kg/kg)
real(kind=real_umphys), intent(out) :: dt_dd(np_field,nlev)
                                             ! dT/dt from DD & evap below base
real(kind=real_umphys), intent(out) :: dq_dd(np_field,nlev)
                                             ! dq/dt from DD & evap below base
real(kind=real_umphys), intent(out) :: du_dd(np_field,nlev)
                                             ! du/dt from DD (m/s2)
real(kind=real_umphys), intent(out) :: dv_dd(np_field,nlev)
                                             ! dv/dt from DD (m/s2)
real(kind=real_umphys), intent(out) :: area_ud(np_field,nlev)
                                             ! updraught fractional area
real(kind=real_umphys), intent(out) :: area_dd(np_field,nlev)
                                             ! downdraught fractional area

! Other diagnostics if requested
real(kind=real_umphys), intent(out) :: up_flux(np_field,nlev) ! mass flux
real(kind=real_umphys), intent(out) :: up_flux_half(np_field,nlev)
                                                 !mass flux on rho levels
real(kind=real_umphys), intent(out) :: dwn_flux(np_field,nlev)
                                             ! Downdraught mass flux (Pa/s)
real(kind=real_umphys), intent(out) :: entrain_up(np_field,nlev)
                                               ! Fractional entrainment rate
real(kind=real_umphys), intent(out) :: detrain_up(np_field,nlev)
                                               ! Fractional detrainment rate
real(kind=real_umphys), intent(out) :: entrain_dwn(np_field,nlev)
                                               ! Fractional entrainment rate
real(kind=real_umphys), intent(out) :: detrain_dwn(np_field,nlev)
                                               ! Fractional detrainment rate
real(kind=real_umphys), intent(out) :: mf_deep(np_field,nlev)
                                               ! mass flux deep
real(kind=real_umphys), intent(out) :: mf_congest(np_field,nlev)
                                               ! mass flux congestus
real(kind=real_umphys), intent(out) :: mf_shall(np_field,nlev)
                                               ! mass flux shallow
real(kind=real_umphys), intent(out) :: mf_midlev(np_field,nlev)
                                               ! mass flux mid-lev
real(kind=real_umphys), intent(out) :: dt_deep(np_field,nlev)
                                               ! dt increment deep (K/s)
real(kind=real_umphys), intent(out) :: dt_congest(np_field,nlev)
                                               ! dt increment congestus (K/s)
real(kind=real_umphys), intent(out) :: dt_shall(np_field,nlev)
                                               ! dt increment shallow (K/s)
real(kind=real_umphys), intent(out) :: dt_midlev(np_field,nlev)
                                               ! dt increment mid-level (K/s)
real(kind=real_umphys), intent(out) :: dq_deep(np_field,nlev)
                                               ! dq increment deep (kg/kg/s)
real(kind=real_umphys), intent(out) :: dq_congest(np_field,nlev)
                                               ! dq increment congestus
real(kind=real_umphys), intent(out) :: dq_shall(np_field,nlev)
                                               ! dq increment shallow (kg/kg/s)
real(kind=real_umphys), intent(out) :: dq_midlev(np_field,nlev)
                                               ! dq increment mid-level
real(kind=real_umphys), intent(out) :: du_deep(np_field,nlev+1)
                                               ! du increment deep (m/s)
real(kind=real_umphys), intent(out) :: du_congest(np_field,nlev+1)
                                                ! du increment congestus (m/s)
real(kind=real_umphys), intent(out) :: du_shall(np_field,nlev+1)
                                               ! du increment shallow (m/s)
real(kind=real_umphys), intent(out) :: du_midlev(np_field,nlev+1)
                                               ! du increment mid-level (m/s)
real(kind=real_umphys), intent(out) :: dv_deep(np_field,nlev+1)
                                               ! dv increment deep (m/s)
real(kind=real_umphys), intent(out) :: dv_congest(np_field,nlev+1)
                                                ! dv increment congestus (m/s)
real(kind=real_umphys), intent(out) :: dv_shall(np_field,nlev+1)
                                               ! dv increment shallow (m/s)
real(kind=real_umphys), intent(out) :: dv_midlev(np_field,nlev+1)
                                               ! dv increment mid-level (m/s)
real(kind=real_umphys), intent(out) :: wqt_flux_sh(np_field,nlev)
                                               ! w'qt' flux (m/s kg/kg)
real(kind=real_umphys), intent(out) :: wthetal_flux_sh(np_field,nlev)
                                                    ! w'thetal' flux  (m/s K)
real(kind=real_umphys), intent(out) :: wthetav_flux_sh(np_field,nlev)
                                                    ! w'thetav' flux  (m/s K)
real(kind=real_umphys), intent(out) :: wql_flux_sh(np_field,nlev)
                                               ! w'ql' flux  (m/s kg/kg)
! Extra variable initialisation for safety
real(kind=real_umphys), intent(out) :: uw_deep(np_field,nlev)
                                               ! X-comp. of deep conv stress
real(kind=real_umphys), intent(out) :: vw_deep(np_field,nlev)
                                               ! Y-comp. of deep conv stress
real(kind=real_umphys), intent(out) :: uw_shall(np_field,nlev)
                                               ! X-comp. of shallow conv stress
real(kind=real_umphys), intent(out) :: vw_shall(np_field,nlev)
                                               ! Y-comp. of shallow conv stress
real(kind=real_umphys), intent(out) :: uw_mid(np_field,nlev)
                                               ! U comp of mid-level stress
real(kind=real_umphys), intent(out) :: vw_mid(np_field,nlev)
                                               ! V comp of mid-level stress

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i, k, i_wt         ! loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='INIT_6A_NPNTS'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise various convection variables, outputs and diagnostics to zero
! at all points. Points not in the compression lists for the respective calls
! to DPCONV, SHCONV, MDCONV or CGCONV will remain zero.

do k = 1,nlev
  do i = 1,npnts
    dqbydt(i,k)              = 0.0
    dthbydt(i,k)             = 0.0
    dqclbydt(i,k)            = 0.0
    dqcfbydt(i,k)            = 0.0
    dcflbydt(i,k)            = 0.0
    dcffbydt(i,k)            = 0.0
    dbcfbydt(i,k)            = 0.0
    ccw(i,k)                 = 0.0
    dt_dd(i,k)               = 0.0
    dq_dd(i,k)               = 0.0
    du_dd(i,k)               = 0.0
    dv_dd(i,k)               = 0.0
    area_ud(i,k)             = 0.0
    area_dd(i,k)             = 0.0
  end do
end do
do k = 1,n_cca_lev
  do i = 1,npnts
    cca(i,k)                 = 0.0
  end do
end do

! Water tracers

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do k = 1,nlev
      do i = 1,npnts
        dqbydt_wtrac(i,k,i_wt)              = 0.0
        dqclbydt_wtrac(i,k,i_wt)            = 0.0
        dqcfbydt_wtrac(i,k,i_wt)            = 0.0
      end do
    end do
  end do
end if

! Initialise diagnostics if requested

if (flg_up_flx) then
  do k = 1,nlev
    do i = 1,npnts
      up_flux(i,k)             = 0.0
    end do
  end do
end if
if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,npnts
      up_flux_half(i,k)        = 0.0
    end do
  end do
end if
if (flg_dwn_flx) then
  do k = 1,nlev
    do i = 1,npnts
      dwn_flux(i,k)            = 0.0
    end do
  end do
end if
if (flg_entr_up) then
  do k = 1,nlev
    do i = 1,npnts
      entrain_up(i,k)          = 0.0
    end do
  end do
end if
if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,npnts
      detrain_up(i,k)          = 0.0
    end do
  end do
end if
if (flg_entr_dwn) then
  do k = 1,nlev
    do i = 1,npnts
      entrain_dwn(i,k)         = 0.0
    end do
  end do
end if
if (flg_detr_dwn) then
  do k = 1,nlev
    do i = 1,npnts
      detrain_dwn(i,k)         = 0.0
    end do
  end do
end if
if (flg_mf_deep) then
  do k = 1,nlev
    do i = 1,npnts
      mf_deep(i,k)         = 0.0
    end do
  end do
end if
if (flg_mf_congest) then
  do k = 1,nlev
    do i = 1,npnts
      mf_congest(i,k)         = 0.0
    end do
  end do
end if
if (flg_mf_shall) then
  do k = 1,nlev
    do i = 1,npnts
      mf_shall(i,k)         = 0.0
    end do
  end do
end if
if (flg_mf_midlev) then
  do k = 1,nlev
    do i = 1,npnts
      mf_midlev(i,k)         = 0.0
    end do
  end do
end if
if (flg_dt_deep) then
  do k = 1,nlev
    do i = 1,npnts
      dt_deep(i,k)         = 0.0
    end do
  end do
end if
if (flg_dt_congest) then
  do k = 1,nlev
    do i = 1,npnts
      dt_congest(i,k)         = 0.0
    end do
  end do
end if
if (flg_dt_shall) then
  do k = 1,nlev
    do i = 1,npnts
      dt_shall(i,k)         = 0.0
    end do
  end do
end if
if (flg_dt_midlev) then
  do k = 1,nlev
    do i = 1,npnts
      dt_midlev(i,k)         = 0.0
    end do
  end do
end if
if (flg_dq_deep) then
  do k = 1,nlev
    do i = 1,npnts
      dq_deep(i,k)         = 0.0
    end do
  end do
end if
if (flg_dq_congest) then
  do k = 1,nlev
    do i = 1,npnts
      dq_congest(i,k)         = 0.0
    end do
  end do
end if
if (flg_dq_shall) then
  do k = 1,nlev
    do i = 1,npnts
      dq_shall(i,k)         = 0.0
    end do
  end do
end if
if (flg_dq_midlev) then
  do k = 1,nlev
    do i = 1,npnts
      dq_midlev(i,k)         = 0.0
    end do
  end do
end if
if (flg_du_deep) then
  do k = 1,nlev+1
    do i = 1,npnts
      du_deep(i,k)         = 0.0
    end do
  end do
end if
if (flg_du_congest) then
  do k = 1,nlev+1
    do i = 1,npnts
      du_congest(i,k)         = 0.0
    end do
  end do
end if
if (flg_du_shall) then
  do k = 1,nlev+1
    do i = 1,npnts
      du_shall(i,k)         = 0.0
    end do
  end do
end if
if (flg_du_midlev) then
  do k = 1,nlev+1
    do i = 1,npnts
      du_midlev(i,k)         = 0.0
    end do
  end do
end if
if (flg_dv_deep) then
  do k = 1,nlev+1
    do i = 1,npnts
      dv_deep(i,k)         = 0.0
    end do
  end do
end if
if (flg_dv_congest) then
  do k = 1,nlev+1
    do i = 1,npnts
      dv_congest(i,k)         = 0.0
    end do
  end do
end if
if (flg_dv_shall) then
  do k = 1,nlev+1
    do i = 1,npnts
      dv_shall(i,k)         = 0.0
    end do
  end do
end if
if (flg_dv_midlev) then
  do k = 1,nlev+1
    do i = 1,npnts
      dv_midlev(i,k)         = 0.0
    end do
  end do
end if
if (flg_wqt_flux) then
  do k = 1,nlev
    do i = 1,npnts
      wqt_flux_sh(i,k)         = 0.0
    end do
  end do
end if
if (flg_wql_flux) then
  do k = 1,nlev
    do i = 1,npnts
      wql_flux_sh(i,k)         = 0.0
    end do
  end do
end if
if (flg_wthetal_flux) then
  do k = 1,nlev
    do i = 1,npnts
      wthetal_flux_sh(i,k)         = 0.0
    end do
  end do
end if
if (flg_wthetav_flux) then
  do k = 1,nlev
    do i = 1,npnts
      wthetav_flux_sh(i,k)         = 0.0
    end do
  end do
end if
!
! Extra variable initialisation for safety
!
if (flg_uw_dp) then
  do k = 1,nlev
    do i = 1,npnts
      uw_deep(i,k)  = 0.0
    end do
  end do
end if
if (flg_vw_dp) then
  do k = 1,nlev
    do i = 1,npnts
      vw_deep(i,k)  = 0.0
    end do
  end do
end if
if (flg_uw_shall) then
  do k = 1,nlev
    do i = 1,npnts
      uw_shall(i,k) = 0.0
    end do
  end do
end if
if (flg_vw_shall) then
  do k = 1,nlev
    do i = 1,npnts
      vw_shall(i,k) = 0.0
    end do
  end do
end if
if (flg_uw_mid) then
  do k = 1,nlev
    do i = 1,npnts
      uw_mid(i,k) = 0.0
    end do
  end do
end if
if (flg_vw_mid) then
  do k = 1,nlev
    do i = 1,npnts
      vw_mid(i,k) = 0.0
    end do
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine init_6a_npnts


subroutine init_6a_moist(npnts, nlev, np_field, ntra, rho_theta,               &
                         rho_dry_theta, mt, dubydt, dvbydt, dtrabydt, l_tracer)

use cv_run_mod,  only: l_mom

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments

integer, intent(in) ::                                                         &
  nlev,                & ! Number of levels
  npnts,               & ! Number of points
  np_field,            & ! Size of fields
  ntra                   ! Number of tracer fields

logical, intent(in) ::                                                         &
  l_tracer               ! Switch for inclusion of tracers

! NOTE - All moist variables passed down to this routine are :
! specific humidities if l_mr_physics = .false.
! mixing ratios       if l_mr_physics = .true.

real(kind=real_umphys), intent(in) ::                                          &
  rho_theta(np_field,nlev),    & ! wet air density on theta levels (kg/m3)
  rho_dry_theta(np_field,nlev)   ! dry air density on theta levels (kg/m3)

real(kind=real_umphys), intent(out) :: mt(npnts,nlev)
                          ! total water mixing ratio (kg/kg)
                          ! rhowet/rhodry = 1 + mt
real(kind=real_umphys), intent(out) :: dubydt(np_field,nlev+1)
                                             ! Increments to U due to CMT (m/s2)
real(kind=real_umphys), intent(out) :: dvbydt(np_field,nlev+1)
                                             ! Increments to V due to CMT (m/s2)
real(kind=real_umphys), intent(out) :: dtrabydt(npnts,nlev,ntra)
                                               ! Increment to tracer due to
                                               ! convection (kg/kg/s)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i, k, ktra       ! loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='INIT_6A_MOIST'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialisation of mt used for mixing/specific conversions
do k=1,nlev
  do i = 1,npnts
    mt(i,k) = (rho_theta(i,k)/rho_dry_theta(i,k)) - 1.0
  end do
end do

if (l_mom) then
  do k=1,nlev+1
    do i = 1,npnts
      dubydt(i,k)              = 0.0
      dvbydt(i,k)              = 0.0
    end do
  end do
end if
if (l_tracer) then
  do ktra = 1,ntra
    do k = 1,nlev
      do i = 1,npnts
        dtrabydt(i,k,ktra)  = 0.0   ! tracer increments
      end do
    end do
  end do
end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine init_6a_moist


end module init_conv6a_mod
