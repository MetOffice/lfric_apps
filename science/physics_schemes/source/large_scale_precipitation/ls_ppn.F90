! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     SUBROUTINES LS_PPN and LS_PPNC------------------------------------
! Purpose:
!    LS_PPN and LS_PPNC:
!       Calculate large-scale (dynamical) precipitation.
!       LS_PPNC is the gather/scatter routine which then
!       calls LSP_ICE.
! Note: in all cases, level counters (incl subscripts) run from
!       1 (lowest model layer) to tdims%k_end
!       (topmost "wet" model layer)
!
! Programming standard: Unified Model Documentation Paper No 3
!
! Documentation: UM Documentation Paper 26.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation
module ls_ppn_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='LS_PPN_MOD'

contains

subroutine ls_ppn(                                                             &
  p_theta_levels, bland,  deltaz, r_theta_levels, r_rho_levels,                &
!-----------------------
! primary fields and
! cloud fractions
!-----------------------
  cf, cfl, cff, precfrac,                                                      &
  rhcrit, f_arr, cos_theta_latitude,                                           &
  lspice_dim1,lspice_dim2,lspice_dim3,                                         &
  rho_r2, rho_dry, q, qcf, qcl, t,                                             &
  qcf2, qrain, qgraup,                                                         &
!------------------------------------
! Wind field for lateral displacement
! of falling ice by shear
!------------------------------------
  u_on_p, v_on_p,                                                              &
!-----------------------
! aerosol variables
!-----------------------
  sea_salt_film, sea_salt_jet,                                                 &
  salt_dim1, salt_dim2, salt_dim3,                                             &
  ukca_cdnc,                                                                   &
  cdnc_dim1, cdnc_dim2, cdnc_dim3,                                             &
  easyaerosol_cdnc,                                                            &
  biogenic,                                                                    &
  snow_depth, land_fract,                                                      &
  so4_acc,                                                                     &
  so4_dis,                                                                     &
  bmass_agd,                                                                   &
  bmass_cld,                                                                   &
  ocff_agd,                                                                    &
  ocff_cld,                                                                    &
  nitr_acc,                                                                    &
  nitr_diss,                                                                   &
  aerosol,                                                                     &
  n_arcl_compnts, i_arcl_compnts, arcl,                                        &
!---------------------------
! Other variables for mphys
!---------------------------
  lsrain,lssnow,lsgraup_mean,                                                  &
  lsrain3d, lssnow3d, lsgraup3d, rainfrac3d,                                   &
  n_drop_tpr, n_drop_3d,                                                       &
  rhc_row_length, rhc_rows,                                                    &
  rhodz_dry, rhodz_moist,                                                      &
! variables for Jules rainfrac code
  ls_rainfrac, land_points, land_index,                                        &
!---------------------------
! COSP variables
!---------------------------
  l_cosp_lsp,                                                                  &
!------------------------------------------------------
! For subgrid orographic rain enhancement
!------------------------------------------------------
  hmteff, zb,                                                                  &
! tnuc_dust field
  tnuc_new,                                                                    &
!---------------------------
! Water tracers
  wtrac_as, wtrac_mp )

     ! General modules
use arcl_mod,              only: npd_arcl_compnts

   ! Stochastic Physics Module
use stochastic_physics_run_mod, only: x1r_rp, l_rp2, i_rp_scheme,              &
                                      i_rp2b, m_ci_rp, rp_idx

   ! Microphysics modules
use mphys_inputs_mod,      only: niters_mp, l_warm_new, l_mcr_qrain,           &
                                 l_mcr_qcf2, l_mcr_qgraup, x1r, l_mcr_precfrac

use mphys_bypass_mod,      only: mphys_mod_top

use mphys_constants_mod,   only: m_ci_sav, l_calc_mp_const, iter_z

use mphys_diags_mod,       only: l_point_diag, mphys_pts

use atm_fields_bounds_mod, only: array_dims, pdims_s, tdims, pdims, tdims_l

use jules_hydrology_mod, only: l_var_rainfrac

   ! Automatic segment size tuning
use tuning_segments_mod,   only: l_autotune_segments, precip_segment_size


   ! Dr Hook modules
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

   ! Large scale precipitation modules
use ls_ppnc_mod,           only: ls_ppnc
use lsp_taper_ndrop_mod,   only: lsp_taper_ndrop
use lspcon_mod,            only: lspcon
use mphys_air_density_mod, only: mphys_air_density
use lsprec_mod,            only: lsprec_set_reals

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use wtrac_mphys_mod,         only: mp_wtrac_type
use wtrac_atm_step_mod,      only: atm_step_wtrac_type

use def_easyaerosol, only: t_easyaerosol_cdnc
implicit none

integer, intent(in) ::                                                         &
  rhc_row_length, rhc_rows,                                                    &
  lspice_dim1,lspice_dim2,lspice_dim3,                                         &
! Dimensions for 3D arrays
    salt_dim1,                                                                 &
                     ! Array dimensions for sea-salt arrays (equal
    salt_dim2,                                                                 &
                     ! either to tdims%i_start:tdims%i_end,
                     ! tdims%j_start:tdims%j_end, and
                     ! 1:tdims%k_end, or
    salt_dim3,                                                                 &
                     ! else 1,1,1, depending on L_SEASALT_CCN).
    cdnc_dim1, cdnc_dim2, cdnc_dim3,                                           &
                     ! UKCA cloud drop number concentration dimensions
    n_arcl_compnts,                                                            &
                      ! Corresponding number of requested components
    i_arcl_compnts(npd_arcl_compnts)
                      ! Array index of each aerosol clim component

real(kind=real_umphys) ::                                                      &
  cf( tdims%i_start:tdims%i_end,                                               &
      tdims%j_start:tdims%j_end,                                               &
      1:tdims%k_end ),                                                         &
                                  ! in Cloud fraction.
  p_theta_levels( tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  1:tdims%k_end ),                                             &
  r_theta_levels( tdims_l%i_start:tdims_l%i_end,                               &
                  tdims_l%j_start:tdims_l%j_end,                               &
                  tdims_l%k_start:tdims_l%k_end ),                             &
                               ! in height of theta levels from centre of Earth
  r_rho_levels( tdims_l%i_start:tdims_l%i_end,                                 &
                tdims_l%j_start:tdims_l%j_end,                                 &
                            1:tdims_l%k_end ),                                 &
                               ! in height of rho levels from centre of earth
    rhcrit( rhc_row_length, rhc_rows,                                          &
            1:tdims%k_end ),                                                   &
                                                 ! in Critical humidity
                                                 ! for cloud formation.
   f_arr(3, tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, tdims%k_end),&
                                        ! in parameters used in fractional
                                        ! standard deviation calculation
   cos_theta_latitude(pdims_s%i_start:pdims_s%i_end,                           &
                      pdims_s%j_start:pdims_s%j_end),                          &
                                        ! in cosine of the latitude
   cfl( tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end,                                             &
        1:tdims%k_end ),                                                       &
                                        !in Cloud liquid fraction.
   cff( tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end,                                             &
        1:tdims%k_end ),                                                       &
                                        !in Cloud ice fraction.
   precfrac( tdims%i_start : tdims%i_end,                                      &
             tdims%j_start : tdims%j_end,                                      &
                         1 : tdims%k_end),                                     &
                                        !INOUT prognostic precipitation fraction
   rho_r2(pdims_s%i_start:pdims_s%i_end,                                       &
          pdims_s%j_start:pdims_s%j_end,                                       &
          pdims_s%k_start:pdims_s%k_end),                                      &
                                ! in Air density * earth radius**2
   rho_dry(pdims_s%i_start:pdims_s%i_end,                                      &
           pdims_s%j_start:pdims_s%j_end,                                      &
           pdims_s%k_start:pdims_s%k_end)
                                ! in Unscaled dry-density / kg m-3

logical :: bland( tdims%i_start : tdims%i_end,                                 &
               tdims%j_start : tdims%j_end )
                                   ! in Land/sea mask

logical, intent(in) ::  l_cosp_lsp ! Switch for COSP LS diagnostics

real(kind=real_umphys), intent(in out) ::                                      &
 tnuc_new( tdims%i_start:tdims%i_end,                                          &
              tdims%j_start:tdims%j_end,                                       &
              1:tdims%k_end ),                                                 &
              ! ice nucln. temp as function of dust (deg cel)

 q( tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end,                                                 &
    1:tdims%k_end ),                                                           &
                                       ! Specific humidity (kg water
 qcf( tdims%i_start:tdims%i_end,                                               &
      tdims%j_start:tdims%j_end,                                               &
      1:tdims%k_end ),                                                         &
                                       ! Cloud ice (kg per kg air).
 qcl( tdims%i_start:tdims%i_end,                                               &
      tdims%j_start:tdims%j_end,                                               &
      1:tdims%k_end ),                                                         &
                                       ! Cloud liquid water (kg per
 qcf2( tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end,                                              &
       1:tdims%k_end ),                                                        &
                                       ! Ice (kg per kg air)
 qrain( tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end,                                             &
        1:tdims%k_end ),                                                       &
                                       ! Rain (kg per kg air)
 qgraup( tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end,                                            &
         1:tdims%k_end ),                                                      &
                                       ! Graupel (kg per kg air)
 t( tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end,                                                 &
    1:tdims%k_end ),                                                           &
                                       ! Temperature (K).
 aerosol( tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end,                                           &
          1:tdims%k_end )
                                       ! 'Murk' tracer aerosol.

real(kind=real_umphys), intent(out) ::                                         &
 rhodz_dry( tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end,                                         &
                        1:tdims%k_end ),                                       &
                                       ! Dry density
 rhodz_moist( tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end,                                       &
                          1:tdims%k_end )
                                       ! Moist density

real(kind=real_umphys),intent(in) ::                                           &
  ! For calculating shear in falling ice cloud fraction calculation.
  u_on_p( pdims%i_start : pdims%i_end,                                         &
          pdims%j_start : pdims%j_end,                                         &
          pdims%k_start : pdims%k_end),                                        &
  v_on_p( pdims%i_start : pdims%i_end,                                         &
          pdims%j_start : pdims%j_end,                                         &
          pdims%k_start : pdims%k_end),                                        &
  ! Aerosol climatology array:
  arcl(   tdims%i_start : tdims%i_end,                                         &
          tdims%j_start : tdims%j_end,                                         &
                      1 : tdims%k_end,                                         &
                      1 : n_arcl_compnts)

real(kind=real_umphys),intent(in) ::                                           &
                                   !Sulphur Cycle tracers (mmr kg/kg)
   so4_acc( tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end,                                         &
            1:tdims%k_end ),                                                   &
   so4_dis( tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end,                                         &
            1:tdims%k_end ),                                                   &

                                   !Biomass smoke tracers
   bmass_agd( tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end,                                       &
              1:tdims%k_end ),                                                 &
   bmass_cld( tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end,                                       &
              1:tdims%k_end ),                                                 &

                                   !Fossil-fuel organic carbon tracers
   ocff_agd( tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end,                                        &
             1:tdims%k_end ),                                                  &
   ocff_cld( tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end,                                        &
             1:tdims%k_end ),                                                  &

                                   !Ammonium nitrate tracers
   nitr_acc( tdims%i_start:tdims%i_end,                                        &
              tdims%j_start:tdims%j_end,                                       &
              1:tdims%k_end ),                                                 &
    nitr_diss( tdims%i_start:tdims%i_end,                                      &
               tdims%j_start:tdims%j_end,                                      &
               1:tdims%k_end )

real(kind=real_umphys) ::                                                      &
    snow_depth( tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end ),                                   &
                                   ! in Snow depth (m)
    land_fract( tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end )
                                   ! in Land fraction

real(kind=real_umphys) ::                                                      &
    sea_salt_film(salt_dim1, salt_dim2, salt_dim3),                            &
                                                       ! (m-3)
    sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)  ! (m-3)

real(kind=real_umphys) ::                                                      &
    biogenic( tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end,                                       &
              1:tdims%k_end )
                                                ! (m.m.r.)

real(kind=real_umphys) :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3)
!         CDNC from UKCA for 2nd indirect effect (m-3)

type (t_easyaerosol_cdnc), intent(in) :: easyaerosol_cdnc
!         CDNC from EasyAerosol for 2nd indirect effect (m-3)

real(kind=real_umphys) ::                                                      &
 lsrain( tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end ),                                          &
                             ! out Surface rainfall rate (kg / sq m /
 lssnow( tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end ),                                          &
                             ! out Surface snowfall rate (kg / sq m /
 lsgraup_mean( tdims%i_start:tdims%i_end,                                      &
               tdims%j_start:tdims%j_end ),                                    &
                             ! out Surface graupel rate (kg / sq m /
 lsgraup( tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end )
                             ! Graupel fall rate (kg/m2/s)

real(kind=real_umphys), intent(out) ::                                         &
                     n_drop_tpr( tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end ),                &
                     n_drop_3d(  tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )
                     ! Tapered droplet number and droplet number
                     ! from autoconversion scheme

real(kind=real_umphys), intent(out) ::                                         &
                      deltaz( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

real(kind=real_umphys) ::                                                      &
    lsrain3d(lspice_dim1,lspice_dim2,lspice_dim3),                             &
                                                      ! out
!                           Rain rate out of each model layer
      lssnow3d(lspice_dim1,lspice_dim2,lspice_dim3),                           &
                                                        ! out
!                           Snow rate out of each model layer
     lsgraup3d(lspice_dim1,lspice_dim2,lspice_dim3),                           &
                                                        ! out
!                           Graupel rate out of each model layer

      rainfrac3d(lspice_dim1,lspice_dim2,lspice_dim3) ! out
!                           rain fraction out of each model layer

! variables for Jules rainfrac code
integer, intent(in) :: land_points
                       ! in No.of land points being processed, can be 0.
integer, intent(in) ::                                                         &
  land_index( land_points )
real(kind=real_umphys), intent(out) ::                                         &
  ls_rainfrac(land_points)
                      ! Rain fraction array on land points to be
                      ! passed to atmos_physics2


! Variable for seeder feeder scheme
real(kind=real_umphys),    intent(in) ::                                       &
                       hmteff(tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end)
real(kind=real_umphys),    intent(in) ::                                       &
                       zb(tdims%i_start : tdims%i_end,                         &
                          tdims%j_start : tdims%j_end)
! Water tracer structures:
!   containing water tracer lsrain and lssnow
type(atm_step_wtrac_type), intent(in out) :: wtrac_as(n_wtrac)
!   containing working arrays of microphysics
type(mp_wtrac_type),       intent(in out) :: wtrac_mp(n_wtrac)

!    Workspace usage ---------------------------------------------------

real(kind=real_umphys) ::                                                      &
         lsrain_mean(tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end )

real(kind=real_umphys) ::                                                      &
         lssnow_mean(tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end )

real(kind=real_umphys) ::                                                      &
         cf_max(tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end )

integer ::                                                                     &
 ix(  tdims%i_len  *                                                           &
      tdims%j_len  , 2 )
                               ! Index for compress/expand.

real(kind=real_umphys) ::                                                      &
        vfall( tdims%i_start:tdims%i_end,                                      &
            tdims%j_start:tdims%j_end )
             ! snow fall velocity (m per s).
real(kind=real_umphys) ::                                                      &
        vfall2( tdims%i_start:tdims%i_end,                                     &
             tdims%j_start:tdims%j_end )
             ! fall velocity for qcf2 (m/s)
real(kind=real_umphys) ::                                                      &
        lssnow2( tdims%i_start:tdims%i_end,                                    &
              tdims%j_start:tdims%j_end )
             ! snowfall rate for qcf2
real(kind=real_umphys) ::                                                      &
        droplet_flux( tdims%i_start:tdims%i_end,                               &
                   tdims%j_start:tdims%j_end )
             ! water drop flux / kg m-2 s-1
real(kind=real_umphys) ::                                                      &
        vfall_rain( tdims%i_start:tdims%i_end,                                 &
                 tdims%j_start:tdims%j_end )
             ! fall velocity for rain (m/s)
real(kind=real_umphys) ::                                                      &
        vfall_graup( tdims%i_start:tdims%i_end,                                &
                  tdims%j_start:tdims%j_end )
             ! fall vel. for graupel (m/s)
real(kind=real_umphys) ::                                                      &
        cttemp( tdims%i_start:tdims%i_end,                                     &
             tdims%j_start:tdims%j_end )
real(kind=real_umphys) ::                                                      &
        rainfrac( tdims%i_start:tdims%i_end,                                   &
               tdims%j_start:tdims%j_end )
real(kind=real_umphys) ::                                                      &
        rainfrac_impr( tdims%i_start:tdims%i_end,                              &
               tdims%j_start:tdims%j_end )
real(kind=real_umphys) ::                                                      &
        frac_ice_above( tdims%i_start:tdims%i_end,                             &
                     tdims%j_start:tdims%j_end )
             ! Cloud ice fraction passed
             ! in layer above

! Precipitation fraction associated with precip falling through
! the model-level interfaces
! (the prognostic precip frac is defined on theta-levels, concurrent
!  with the prognostic precip mass).
real(kind=real_umphys) :: precfrac_fall( tdims%i_start:tdims%i_end,            &
                                         tdims%j_start:tdims%j_end )

real(kind=real_umphys) :: iter_eta    ! eta value at which to
                                       ! start iterative melting

!  Physical constants -------------------------------------------------
real(kind=real_umphys) :: cfmin
parameter (                                                                    &
 cfmin=1.0e-3                                                                  &
                         ! Used for LS_PPNC  compress.
)
!  Define local variables ----------------------------------------------
integer :: i,k,it,i_wt,                                                        &
                      ! Loop counters: I - horizontal field index;
                      ! K - vertical level index.
           kp1,                                                                &
                      ! Index of level above: k=k+1, apart from
                      ! when k=_dims%end when kp1=k.
           n,                                                                  &
                      ! "nval" for WHEN routine.
           i_jul, j_jul
                      ! i and j local indexes for Jules flag section.
                      ! Allows j to be a parameter below.

integer, parameter :: j = 1 ! Arrays passed in by LFRic have a
                            ! size of 1 in the j dimension.

real(kind=real_umphys) :: work
                    ! work variable

logical :: l_3ddiag ! Flag to determine if we want 3d diagnostics

logical :: l_use_gridbox( tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end )
           ! flag for whether grid box is used in microphysics or not

real(kind=real_umphys) :: one_over_niters_mp  ! 1./niters_mp

! Structure for dimensions of the q arrays to pass into mphys_air_density
type(array_dims) :: qdims

!Automatic segment size tuning

!  Variables for Dr Hook:

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LS_PPN'


!-----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Set up automatic segment tuning

! Determine if we want 3d diagnosics
l_3ddiag = (lspice_dim1  ==  tdims%i_len  .and.                                &
   lspice_dim2  ==  tdims%j_len  .and.                                         &
   lspice_dim3  ==  tdims%k_end  )

!----------------------------------------------------------------------

! Define CX and CONSTP values

! Determine whether m_ci_rp or x1r_rp have been udpated

if (l_calc_mp_const .or. m_ci_rp(rp_idx) /= m_ci_sav .or.                      &
    (l_rp2 .and. i_rp_scheme == i_rp2b .and. x1r /= x1r_rp(rp_idx))) then

  ! If the microphysics constants and number of microphysics iterations
  ! (niters_mp) are not set we can compute them
  ! and save to a module. However, if m_ci (random parameters)
  ! or x1r (random parameters version 2b) changes, some constants
  ! will change with it, so we need to recompute the constants in
  ! this case. Continuation runs (e.g. climate) will also need to
  ! recalculate the constants.
  call lspcon()

  !Set the values of scheme constants in the lsprec_mod
  call lsprec_set_reals()

end if

! Set dimensions of the q arrays to pass into mphys_air_density
! (same as tdims but without the zero-level)
qdims = tdims
qdims%k_start = 1

! Calculate air density (use a standard routine for both CASIM
! and Wilson and Ballard microphysics )
call mphys_air_density( r_theta_levels, r_rho_levels, rho_dry, rho_r2, qdims,  &
                        q, qcl, qcf, qcf2, qrain, qgraup,                      &
                        rhodz_dry, rhodz_moist, deltaz )


!-------------------------------------------------------------
! Calculation of cloud droplet number. This is now calculated
! here for all models and not in lsp_autoc as was done in
! older versions of the Unified Model
!-------------------------------------------------------------

call lsp_taper_ndrop(                                                          &
               ! (Full) Aerosol tracers
                        so4_acc, so4_dis, sea_salt_film,                       &
                        biogenic, sea_salt_jet, bmass_agd,                     &
                        bmass_cld, ocff_agd, ocff_cld,                         &
                        nitr_acc, nitr_diss,                                   &
                        n_arcl_compnts,                                        &
                        i_arcl_compnts, arcl,                                  &
               ! Murk aerosol
                        aerosol,                                               &
               ! CDNC from UKCA
                        ukca_cdnc,                                             &
                        cdnc_dim1, cdnc_dim2, cdnc_dim3,                       &
               ! CDNC from EasyAerosol
                        easyaerosol_cdnc,                                      &
               ! Other parameters
                        rhodz_dry,  rhodz_moist,                               &
                        deltaz,                                                &
                        snow_depth, land_fract,                                &
               ! Output parameters
                        n_drop_tpr                                             &
                             )

!-----------------------------------------------------------------------
!  2. Loop round levels from top down (counting bottom level as level 1,
!     as is standard in the Unified model).
!-----------------------------------------------------------------------


!-----------------------------------------------
! Setup for iterative metlting
! Define constants outside of K loop
!-----------------------------------------------

! calculate level independent eta value
! iter_z is theta_height for level 13 in 38 level set.
! hence is now independent of number of levels

iter_eta = iter_z / mphys_mod_top

!-----------------------------------------------------------------------
! Internal structure.
! 2a. Initialise outside of iterative loop.
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED( tdims, lsrain_mean, lssnow_mean, lsgraup_mean, lsrain3d,         &
!$OMP         lssnow3d, lsgraup3d, rainfrac3d, lspice_dim1,                    &
!$OMP         lspice_dim3, l_wtrac, n_wtrac, wtrac_as )                        &
!$OMP private( i, k, i_wt )
!$OMP  do SCHEDULE(STATIC)
do i = tdims%i_start, tdims%i_end
  lsrain_mean(i,j) =0.0
  lssnow_mean(i,j) =0.0
  lsgraup_mean(i,j)=0.0
end do ! Loop over points,i
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do k = 1, lspice_dim3
  do i = 1, lspice_dim1
    lsrain3d(i,j,k)=0.0
    lssnow3d(i,j,k)=0.0
    lsgraup3d(i,j,k)=0.0
    rainfrac3d(i,j,k)=0.0
  end do ! Loop over points,i
end do ! Loop over points,k
!$OMP end do NOWAIT

if (l_wtrac) then
  do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
    do i = tdims%i_start, tdims%i_end
      wtrac_as(i_wt)%ls_rain(i,j) = 0.0
      wtrac_as(i_wt)%ls_snow(i,j) = 0.0
    end do
!$OMP end do NOWAIT
  end do
end if
!$OMP end PARALLEL

ls_rainfrac(:)=0.0
one_over_niters_mp=1.0/niters_mp

do it = 1, niters_mp ! Substep outside of column

  !-----------------------------------------------------------------------
  !  Internal structure - moved inside BS iter loop
  !  2.b Initialise rain and snow to zero.
  !   Initialise scavenged amounts of S Cycle tracers to 0 for full field
  !-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED( tdims, lsrain, lssnow, lssnow2, lsgraup, droplet_flux,           &
!$OMP         cttemp, rainfrac, rainfrac_impr, precfrac_fall, frac_ice_above,  &
!$OMP         vfall, vfall2, vfall_rain, vfall_graup, cf_max, n_drop_3d,       &
!$OMP         l_wtrac, n_wtrac, wtrac_mp)                                      &
!$OMP  private( i, k, i_wt )
!$OMP do SCHEDULE(STATIC)
  do i = tdims%i_start, tdims%i_end

    lsrain(i,j)=0.0
    lssnow(i,j)=0.0
    lssnow2(i,j)=0.0
    lsgraup(i,j)=0.0
    droplet_flux(i,j)=0.0
    cttemp(i,j)=0.0
    rainfrac(i,j)=0.0
    rainfrac_impr(i,j)=0.0
    precfrac_fall(i,j)=0.0
    frac_ice_above(i,j)=0.0
    vfall(i,j)=0.0
    vfall2(i,j)=0.0
    vfall_rain(i,j)=0.0
    vfall_graup(i,j)=0.0
    cf_max(i,j)=-huge(1.0)
  end do ! Loop over points,i
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do i = tdims%i_start, tdims%i_end
      ! Initialise n_drop_3d to zero before passing down code tree
      n_drop_3d(i,j,k)=0.0
    end do ! Loop over points,i
  end do ! Loop over k
!$OMP end do NOWAIT

  if (l_wtrac) then
    do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
      do i = tdims%i_start, tdims%i_end
        wtrac_mp(i_wt)%lsrain(i,j)       = 0.0
        wtrac_mp(i_wt)%lssnow(i,j)       = 0.0
        wtrac_mp(i_wt)%droplet_flux(i,j) = 0.0
      end do ! Loop over points,i
!$OMP end do NOWAIT
    end do ! Loop over wt
  end if
!$OMP end PARALLEL

  do k = tdims%k_end, 1, -1
!$OMP PARALLEL DEFAULT(SHARED) private( i, work )
!$OMP do SCHEDULE(STATIC)
    do i = tdims%i_start, tdims%i_end

      !-----------------------------------------------------------------------
      !  2.5 Form index IX to gather/scatter variables in LS_PPNC
      !-----------------------------------------------------------------------

      !  Set index where cloud fraction > CFMIN or where non-zero pptn
      !  Note: whenimd is functionally equivalent to WHENILE (but autotasks).

      ! Set up if statement to determine whether to call the
      ! microphysics code for this grid box (i.e. if there is
      ! already condensate in the grid box or there is
      ! precipitation about to fall into the grid box)
      work = qcf(i,j,k)

      ! Include extra microphysics variables if in use
      if (l_mcr_qcf2 )  work = work + qcf2(i,j,k) + lssnow2(i,j)
      if (l_mcr_qrain)  work = work + qrain(i,j,k)
      if (l_mcr_qgraup) work = work + qgraup(i,j,k)+lsgraup(i,j)
      work = work + droplet_flux(i,j)   ! droplet settling

      if (cfl(i,j,k) > cfmin .or.                                            &
          (lsrain(i,j)+lssnow(i,j)) > 0.0 .or. work > 0.0) then
            ! include this grid box.
            ! Strictly speaking the CFL > CFMIN clause is too
            ! restrictive since ice nucleation does not require
            ! liquid water, but the code would be very messy.
        l_use_gridbox(i,j) = .true.
            ! Note that mphys is done on this point if diagnostic is
            ! requested
        if (l_point_diag) mphys_pts(i,j,k) = .true.
      else
        l_use_gridbox(i,j) = .false.
      end if

    end do ! tdims%i
!$OMP end do NOWAIT

    if ( .not. l_mcr_precfrac ) then
      ! If remembering the rain fraction from previous timesteps,
      ! the following setup is not needed
!$OMP do SCHEDULE(STATIC)
      do i = tdims%i_start, tdims%i_end
        ! set up improved rain fraction...
        rainfrac_impr(i,j) = rainfrac(i,j)
      end do ! tdims%i
!$OMP end do NOWAIT

      if (l_mcr_qrain) then
!$OMP do SCHEDULE(STATIC)
        do i = tdims%i_start, tdims%i_end
          ! CF_MAX contains the maximum CF values of all the levels
          ! higher than k.
          if (k /= tdims%k_end) then
            cf_max(i,j) = max(cf_max(i,j), cf(i,j,k+1))
          end if

          if ( qrain(i,j,k) > 0.0 ) then
            if (rainfrac_impr(i,j) == 0.0) then
              ! if cloud is present, assume this is a good proxy for the
              ! rain fraction
              rainfrac_impr(i,j) = max( cf(i,j,k), cf_max(i,j))
            end if
            if (rainfrac_impr(i,j) == 0.0) then
              ! otherwise set to 0.5
              rainfrac_impr(i,j) = 0.5
              ! set a lower limit for stability
            else if (rainfrac_impr(i,j) <= 0.01) then
              rainfrac_impr(i,j) = 0.01
            end if
          end if
        end do ! tdims%i
!$OMP end do NOWAIT
      end if

      if (l_warm_new) then
!$OMP do SCHEDULE(STATIC)
        do i = tdims%i_start, tdims%i_end
          ! Always set rain fraction back to improved rain fraction
          rainfrac(i,j) = rainfrac_impr(i,j)
        end do ! tdims%i
!$OMP end do NOWAIT
      end if
    end if  ! ( .not. l_mcr_precfrac )
!$OMP end PARALLEL

    n=0
    do i = tdims%i_start, tdims%i_end
      if (l_use_gridbox(i,j)) then
        n = n + 1
        ix(n,1) = i
        ix(n,2) = j
      end if
    end do ! tdims%i


    if (n > 0) then
      ! kp1 is index of model level above, unless we are the model
      ! top in which case it is set to k. Used to calc vertical wind shear.
      if (k == tdims%k_end) then
        kp1=k
      else
        kp1=k+1
      end if

      call ls_ppnc(k,ix,n,                                                     &
                   lsrain,lssnow,lssnow2,lsgraup,droplet_flux,                 &
                   cf(1,1,k),cfl(1,1,k),cff(1,1,k),                            &
                   qcf(1,1,k),qcl(1,1,k),tnuc_new(1,1,k),t(1,1,k),             &
                   qcf2(1,1,k),qrain(1,1,k),qgraup(1,1,k),                     &
                   n_drop_tpr(1,1,k), n_drop_3d(1,1,k),                        &
                   aerosol(1,1,k),                                             &
                   hmteff, zb,                                                 &
                   q(1,1,k), p_theta_levels(1,1,k), r_theta_levels,            &
                   deltaz(1,1,k), rhodz_dry(1,1,k), rhodz_moist(1,1,k),        &
                   rhc_row_length, rhc_rows,                                   &
                   bland, rhcrit(1,1,k), f_arr(:,1,1,k), cos_theta_latitude,   &
                   vfall, vfall2, vfall_rain, vfall_graup,                     &
                   frac_ice_above,                                             &
                   cttemp, rainfrac, rainfrac_impr,                            &
                   precfrac(1,1,k), precfrac_fall,                             &
                   niters_mp,                                                  &
                   u_on_p(1,1,k),   v_on_p(1,1,k),                             &
                   u_on_p(1,1,kp1), v_on_p(1,1,kp1), k,                        &
                   l_cosp_lsp, wtrac_mp                                        &
                   )
    end if


    ! Copy rainfall and snowfall rates to 3D fields for diagnostic output

    if ( l_3ddiag ) then
      ! Only copy rain and snow to 3D fields if arrays are dimensionalized.
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP SHARED( tdims, lsrain3d, lsrain, droplet_flux, one_over_niters_mp,       &
!$OMP         lssnow3d, lsgraup3d, lssnow, lssnow2, lsgraup, rainfrac3d,       &
!$OMP         rainfrac, k)                                                     &
!$OMP private(i)
      do i = tdims%i_start, tdims%i_end

        lsrain3d(i, j, k)   = lsrain3d(i, j, k) +                            &
            (lsrain(i, j) + droplet_flux(i, j))*one_over_niters_mp

        lssnow3d(i, j, k)   = lssnow3d(i, j, k) +                            &
            (lssnow(i, j) + lssnow2(i, j) + lsgraup(i, j))                   &
            *one_over_niters_mp

        lsgraup3d(i, j, k)  = lsgraup3d(i, j, k) +                           &
                                ( lsgraup(i, j) * one_over_niters_mp )

        rainfrac3d(i, j, k) = rainfrac3d(i, j, k) +                          &
            rainfrac(i, j)*one_over_niters_mp

      end do ! tdims%i
!$OMP end PARALLEL do
    end if

  end do ! Loop over k

  ! save rain fraction on land points to pass to Jules
  ! loop over K is backwards, hence rainfrac now contains the level 1
  ! (surface) value
!$OMP PARALLEL DEFAULT(SHARED) private( k, i_wt, j_jul, i_jul )
  if (l_var_rainfrac) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, land_points
      j_jul = (land_index(k)-1)/ ( tdims%i_len ) + 1
      i_jul = land_index(k) - (j_jul-1)*( tdims%i_len )
      ls_rainfrac(k) = ls_rainfrac(k) +                                      &
          rainfrac(i_jul,j_jul)*one_over_niters_mp
    end do
!$OMP end do NOWAIT
  end if ! Loop over k

  ! If substepping outside of loop over K, then need to accumulate (mean)
  ! precip rates...

!$OMP do SCHEDULE(STATIC)
  do i = tdims%i_start, tdims%i_end

    lsrain_mean(i, j) = lsrain_mean(i, j)                                    &
        + (lsrain(i, j) + droplet_flux(i, j))*one_over_niters_mp

    ! Add together ice crystals, snow aggregates and graupel
    ! for surface snow rate (kg/m2/s)

    lssnow_mean(i, j) = lssnow_mean(i, j)                                    &
                    + (lssnow(i, j)  + lssnow2(i, j)                         &
                    + lsgraup(i, j)) * one_over_niters_mp

    lsgraup_mean(i, j) = lsgraup_mean(i, j)                                  &
                    + lsgraup(i, j) * one_over_niters_mp

  end do ! tdims%i
!$OMP end do NOWAIT

  ! Accumulate rain and snow rates for water tracers
  if (l_wtrac) then
    do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
      do i = tdims%i_start, tdims%i_end
        wtrac_as(i_wt)%ls_rain(i,j) = wtrac_as(i_wt)%ls_rain(i,j)            &
                                      + (wtrac_mp(i_wt)%lsrain(i,j)          &
                                        + wtrac_mp(i_wt)%droplet_flux(i,j))  &
                                      * one_over_niters_mp

        ! No qcf2 or graupel for water tracers
        wtrac_as(i_wt)%ls_snow(i,j) = wtrac_as(i_wt)%ls_snow(i,j)            &
                          + wtrac_mp(i_wt)%lssnow(i,j) * one_over_niters_mp
      end do ! tdims%i
!$OMP end do NOWAIT
    end do ! loop over wt
  end if

!$OMP end PARALLEL
end do ! Outer substepping loop over it

! Recover meaned precip rates
! (No need to do this for water tracers as mean value already stored in
!  lsrain and lssnow water tracer equivalents)
!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED( tdims, lssnow, lssnow_mean, lsrain, lsrain_mean)                &
!$OMP  private( i )
do i = tdims%i_start, tdims%i_end

  lssnow(i,j) = lssnow_mean(i,j)
  lsrain(i,j) = lsrain_mean(i,j)

end do ! tdims%i
!$OMP end PARALLEL do

!If autotuning is active, decide what to do with the
!segment size and report the current status.

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_ppn

end module ls_ppn_mod
