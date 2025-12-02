! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Set flags for use in convection code

module cv_stash_flg_mod

! Description:
!   Module containing stash flags used in top level convection routine.
!
! Method:
!   Declares all flags and contains a subroutine used to set the flags.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:

use science_fixes_mod,  only: l_fix_conv_diags_var

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
save

!=====================================================================
! Logical flags used to control output/setting of various convection
! diagnostics
!=====================================================================

logical ::                                                                     &
  l_apply_diag         ! true if sub-step for calculating diagnostics

logical ::                                                                     &
  flg_up_flx         & ! stash flag for updraught mass flux (theta-levels)
 ,flg_up_flx_half    & ! stash flag for updraught mass flux (rho-levels)
 ,flg_dwn_flx        & ! stash flag for downdraught mass flux (theta-levels)
 ,flg_dwn_flx_half   & ! stash flag for downdraught mass flux (rho-levels)
 ,flg_entr_up        & ! stash flag for updraught entrainment
 ,flg_entr_dwn       & ! stash flag for downdraught entrainmnent
 ,flg_detr_up        & ! stash flag for updraught detrainment
 ,flg_detr_dwn       & ! stash flag for downdraught detrainmnent
 ,flg_conv_rain_3d   & ! stash flag for 3d conv rainfall
 ,flg_conv_snow_3d     ! stash flag for 3d conv snowfall

! CMT diagnostics

logical ::                                                                     &
  flg_uw_dp          & ! stash flag for deep x-comp stress
 ,flg_vw_dp          & ! stash flag for deep y-comp stress
 ,flg_uw_shall       & ! stash flag for shallow x stress
 ,flg_vw_shall       & ! stash flag for shallow y stress
 ,flg_uw_mid         & ! stash flag for mid-level x stress
 ,flg_vw_mid           ! stash flag for mid-level y stress
! Increment diagnostics:

logical ::                                                                     &
  l_qcl_incr_cinh    & ! liquid cloud condensate qCL (inhom)
 ,l_qcf_incr_cinh    & ! frozen cloud condensate qCF (inhom)
 ,l_cfl_incr_cinh    & ! liquid cloud amount cf_liquid (inhom)
 ,l_cff_incr_cinh    & ! frozen cloud amount cf_frozen (inhom)
 ,l_bcf_incr_cinh    & ! total (bulk) cloud amount bulk_cf (inhom)

 ,l_theta_incr_conv  & ! potential temperature increment
 ,l_T_incr_conv      & ! temperature
 ,l_q_incr_conv      & ! humidity
 ,l_qcl_incr_conv    & ! liquid cloud condensate qCL
 ,l_qcf_incr_conv    & ! frozen cloud condensate qCF
 ,l_qcf2_incr_conv   & ! frozen cloud condensate qCF2
 ,l_qgraup_incr_conv & ! frozen condensate qGRAUP
 ,l_qrain_incr_conv  & ! liquid condensate qRAIN
 ,l_cfl_incr_conv    & ! liquid cloud amount cf_liquid
 ,l_cff_incr_conv    & ! frozen cloud amount cf_frozen
 ,l_bcf_incr_conv    & ! total (bulk) cloud amount bulk_cf
 ,l_T_conv_only      & ! temperature from convection only
 ,l_q_conv_only        ! humidity from convection only

! DD increments

logical ::                                                                     &
  flg_dt_dd          & ! stash flag for dT/dt from DD
 ,flg_dq_dd          & ! stash flag for dq/dt from DD
 ,flg_du_dd          & ! stash flag for du/dt from DD
 ,flg_dv_dd            ! stash flag for dv/dt from DD

! Fractional areas
logical ::                                                                     &
  flg_area_ud        & ! stash flag for updraught fractional area
 ,flg_area_dd          ! stash flag for downdraught fractional area


! 5A turbulence based schemes

logical ::                                                                     &
  flg_wqt_flux       & ! stash flag for w'qt'
 ,flg_wql_flux       & ! stash flag for w'ql'
 ,flg_wthetal_flux   & ! stash flag for w'thetal'
 ,flg_wthetav_flux   & ! stash flag for w'thetav'

 ,flg_deep_tops      & ! stash flag for deep tops

 ,flg_mf_deep        & ! stash flag for deep mass flux
 ,flg_mf_congest     & ! stash flag for congestus mass flux
 ,flg_mf_shall       & ! stash flag for shallow mass flux
 ,flg_mf_midlev      & ! stash flag for mid_level mass flux

 ,flg_dt_deep        & ! stash flag for deep dT
 ,flg_dt_congest     & ! stash flag for congestus dT
 ,flg_dt_shall       & ! stash flag for shallow dT
 ,flg_dt_midlev      & ! stash flag for mid_level dT

 ,flg_dq_deep        & ! stash flag for deep dq
 ,flg_dq_congest     & ! stash flag for congestus dq
 ,flg_dq_shall       & ! stash flag for shallow dq
 ,flg_dq_midlev      & ! stash flag for mid_level dq

 ,flg_du_deep        & ! stash flag for deep du
 ,flg_du_congest     & ! stash flag for congestus du
 ,flg_du_shall       & ! stash flag for shallow du
 ,flg_du_midlev      & ! stash flag for mid_level du

 ,flg_dv_deep        & ! stash flag for deep dv
 ,flg_dv_congest     & ! stash flag for congestus dv
 ,flg_dv_shall       & ! stash flag for shallow dv
 ,flg_dv_midlev        ! stash flag for mid_level dv

! Flag for Simpson & Wiggert w-equation diagnostic
logical ::                                                                     &
  flg_w_eqn            ! stash flag for w_eqn


! CoMorph - only diagnostics

logical ::                                                                     &
  flg_par_radius_up       & ! stash flag for updraught parcel radius
 ,flg_par_radius_dwn      & ! stash flag for downdraught parcel radius
 ,flg_freq_up             & ! stash flag for updraught on rho levels
 ,flg_freq_dwn            & ! stash flag for downdraught on rho levels
 ,flg_par_meanup_dtv      & ! stash flag for updraught parcel mean Tv excess
 ,flg_par_meandn_dtv      & ! stash flag for downdraught parcel mean Tv excess
 ,flg_par_meanup_rhl      & ! stash flag for updraught parcel mean RH
 ,flg_par_meandn_rhl      & ! stash flag for downdraught parcel mean RH
 ,flg_par_meanup_q        & ! stash flag for updraught parcel mean q
 ,flg_par_meandn_q        & ! stash flag for downdraught parcel mean q
 ,flg_par_meanup_qcl      & ! stash flag for updraught parcel mean qcl
 ,flg_par_meandn_qcl      & ! stash flag for downdraught parcel mean qcl
 ,flg_par_meanup_qcf      & ! stash flag for updraught parcel mean qcf
 ,flg_par_meandn_qcf      & ! stash flag for downdraught parcel mean qcf
 ,flg_par_meanup_qrain    & ! stash flag for updraught parcel mean qrain
 ,flg_par_meandn_qrain    & ! stash flag for downdraught parcel mean qrain
 ,flg_par_meanup_qgr      & ! stash flag for updraught parcel mean qgr
 ,flg_par_meandn_qgr      & ! stash flag for downdraught parcel mean qgr
 ,flg_par_meanup_qsnow    & ! stash flag for updraught parcel mean qsnow
 ,flg_par_meandn_qsnow    & ! stash flag for downdraught parcel mean qsnow
 ,flg_par_meanup_cfl      & ! stash flag for updraught parcel mean cfl
 ,flg_par_meandn_cfl      & ! stash flag for downdraught parcel mean cfl
 ,flg_par_meanup_cff      & ! stash flag for updraught parcel mean cff
 ,flg_par_meandn_cff      & ! stash flag for downdraught parcel mean cff
 ,flg_par_meanup_cfb      & ! stash flag for updraught parcel mean cfb
 ,flg_par_meandn_cfb      & ! stash flag for downdraught parcel mean cfb
 ,flg_par_meanup_u        & ! stash flag for updraught parcel mean u wind
 ,flg_par_meandn_u        & ! stash flag for downdraught parcel mean u wind
 ,flg_par_meanup_v        & ! stash flag for updraught parcel mean v wind
 ,flg_par_meandn_v        & ! stash flag for downdraught parcel mean v wind
 ,flg_par_meanup_w        & ! stash flag for updraught parcel mean w wind
 ,flg_par_meandn_w        & ! stash flag for downdraught parcel mean w wind
 ,flg_par_meanup_t        & ! stash flag for updraught parcel mean temperature
 ,flg_par_meandn_t        & ! stash flag for downdraught parcel mean temperature
 ,flg_par_coreup_dtv      & ! stash flag for updraught parcel core Tv excess
 ,flg_par_coredn_dtv      & ! stash flag for downdraught parcel core Tv excess
 ,flg_par_coreup_rhl      & ! stash flag for updraught parcel core RH
 ,flg_par_coredn_rhl        ! stash flag for downdraught parcel core RH

logical ::                & ! flags for turbulent perturbations
  flg_turb_pert_u         & ! stash flag for u wind
 ,flg_turb_pert_v         & ! stash flag for v wind
 ,flg_turb_pert_w         & ! stash flag for w wind
 ,flg_turb_pert_t         & ! stash flag for temperature
 ,flg_turb_pert_q         & ! stash flag for water vapour
 ,flg_turb_radius           ! stash flag for turbluent radius


logical ::                & ! flags for sub-region diags
  flg_dry_fraction        & ! stash flag for gridbox dry fraction
 ,flg_liq_fraction        & ! stash flag for gridbox liquid cloud fraction
 ,flg_mix_fraction        & ! stash flag for gridbox mixed phase fraction
 ,flg_icr_fraction        & ! stash flag for gridbox ice/rain fraction
 ,flg_dry_temp            & ! stash flag for gridbox dry temp
 ,flg_liq_temp            & ! stash flag for gridbox liquid cloud temp
 ,flg_mix_temp            & ! stash flag for gridbox mixed phase temp
 ,flg_icr_temp            & ! stash flag for gridbox ice/rain temp
 ,flg_dry_q_vap           & ! stash flag for gridbox dry vapour
 ,flg_liq_q_vap           & ! stash flag for gridbox liquid cloud vapour
 ,flg_mix_q_vap           & ! stash flag for gridbox mixed phase vapour
 ,flg_icr_q_vap           & ! stash flag for gridbox ice/rain vapour
 ,flg_dry_rhl             & ! stash flag for gridbox dry RH
 ,flg_liq_rhl             & ! stash flag for gridbox liquid cloud  RH
 ,flg_mix_rhl             & ! stash flag for gridbox mixed phase RH
 ,flg_icr_rhl               ! stash flag for gridbox ice/rain RH


character(len=*), parameter, private :: ModuleName='CV_STASH_FLG_MOD'

contains

!=========================================================================
! Full model version - sets flags according to user requirements
! i.e. stash settings held in array sf
!=========================================================================
subroutine set_convection_output_flags()

use model_domain_mod,   only: model_type, mt_single_column, mt_lfric
use stash_array_mod,    only: sf
use cosp_input_mod,     only: i_cosp_version
use cv_run_mod,         only: i_convection_vn, i_cv_comorph
use mphys_inputs_mod,   only: l_mcr_qcf2, graupel_option, no_graupel,          &
                              l_mcr_qrain
use comorph_constants_mod, only: l_cv_cloudfrac, l_par_core,                   &
                                 l_cv_rain, l_cv_snow, l_cv_graup

implicit none

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SET_CONVECTION_OUTPUT_FLAGS'

!---------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if ((model_type /= mt_single_column) .and.  (model_type /= mt_lfric)) then

  flg_conv_rain_3d = (sf(227,5) .and. l_apply_diag) .or. (i_cosp_version > 0)
  flg_conv_snow_3d = (sf(228,5) .and. l_apply_diag) .or. (i_cosp_version > 0)

  flg_up_flx = .true.  ! Always needed in atmos_physics2
  flg_up_flx_half = ( ( sf(249,5) .or. sf(246,5) ) .and. l_apply_diag )        &
               .or. ( flg_up_flx .and. i_convection_vn == i_cv_comorph )
  ! Note: flux on half-levels is needed in order to compute the flux
  ! on theta-levels if using the CoMorph convection scheme, since
  ! CoMorph only outputs the fluxes on rho-levels.

  flg_dwn_flx = .true.  ! Always needed in atmos_physics2
  flg_dwn_flx_half = ( flg_dwn_flx .and. i_convection_vn == i_cv_comorph ) .or.&
                       ( sf(616,5) .and. i_convection_vn == i_cv_comorph )

  ! Need to store downdraught mass flux rho on half levels if want on theta
  ! levels if using the CoMorph convection scheme.

  flg_entr_up  = sf(252,5) .and. l_apply_diag
  flg_detr_up  = sf(253,5) .and. l_apply_diag
  flg_entr_dwn = sf(254,5) .and. l_apply_diag
  flg_detr_dwn = sf(255,5) .and. l_apply_diag

  flg_uw_dp    = sf(258,5) .and. l_apply_diag
  flg_vw_dp    = sf(259,5) .and. l_apply_diag
  flg_uw_shall = sf(260,5) .and. l_apply_diag
  flg_vw_shall = sf(261,5) .and. l_apply_diag
  flg_uw_mid   = sf(263,5) .and. l_apply_diag
  flg_vw_mid   = sf(264,5) .and. l_apply_diag

  ! PC2 increments

  l_qcl_incr_cinh = ( sf(163,5) .or. sf(140,5) .or. sf(141,5) .or.             &
                      sf(150,5) .or. sf(151,5) ) .and. l_apply_diag
  l_qcf_incr_cinh = (sf(164,5) .or. (sf(142,5) .or. sf(143,5)))                &
                    .and. l_apply_diag
  l_bcf_incr_cinh = sf(172,5) .and. l_apply_diag
  l_cfl_incr_cinh = ( sf(173,5) .or. sf(146,5) .or. sf(147,5) .or.             &
                      sf(156,5) .or. sf(157,5) .or. sf(158,5) )                &
                    .and. l_apply_diag
  l_cff_incr_cinh = ( sf(174,5) .or. sf(148,5) .or. sf(149,5) )                &
                    .and. l_apply_diag
  l_bcf_incr_conv = sf(192,5) .and. l_apply_diag
  l_cfl_incr_conv = (sf(193,5) .or. (sf(156,5) .or. sf(157,5) .or. sf(158,5))) &
                    .and. l_apply_diag
  l_cff_incr_conv = sf(194,5) .and. l_apply_diag

  ! Convection increments

  l_theta_incr_conv = .false.  ! Theta increment only output by the SCM
  l_T_incr_conv   = ( sf(181,5) .or. sf(187,5) .or. sf(161,5) )                &
                    .and. l_apply_diag
  l_q_incr_conv   = ( sf(182,5) .or. sf(188,5) .or. sf(162,5) )                &
                    .and. l_apply_diag
  l_qcl_incr_conv = ( sf(183,5) .or. sf(150,5) .or. sf(151,5) .or. sf(152,5) ) &
                    .and. l_apply_diag
  l_qcf_incr_conv = sf(184,5) .and. l_apply_diag

  ! Increments to optional condensate species
  if ( i_convection_vn == i_cv_comorph .and. l_apply_diag ) then
    ! CoMorph convection scheme produces increments to any optional
    ! condensed water fields that are in use:
    l_qcf2_incr_conv   = l_mcr_qcf2 .and. sf(191,5)
    l_qgraup_incr_conv = ( graupel_option > no_graupel ) .and. sf(190,5)
    l_qrain_incr_conv  = l_mcr_qrain .and. sf(189,5)

  else
    ! Other convection schemes do not calculate increments for these fields:
    ! so not allowed to request them
    l_qcf2_incr_conv   = .false.
    l_qgraup_incr_conv = .false.
    l_qrain_incr_conv  = .false.
  end if

  l_T_conv_only   = ((l_fix_conv_diags_var .and. sf(187,5))                    &
                    .or. sf(161,5)) .and. l_apply_diag
  l_q_conv_only   = ((l_fix_conv_diags_var .and. sf(188,5))                    &
                    .or. sf(162,5)) .and. l_apply_diag


  ! 5A convection code

  flg_wqt_flux     = (sf(290,5) .or. sf(304,5) .or. sf(306,5) .or.             &
                      sf(429,5) .or. sf(430,5) .or. sf(431,5) .or. sf(432,5))  &
                     .and. l_apply_diag

  flg_wql_flux     = sf(291,5) .and. l_apply_diag

  flg_wthetal_flux = (sf(292,5) .or. sf(305,5) .or. sf(307,5) .or.             &
                      sf(425,5) .or. sf(426,5) .or. sf(427,5) .or. sf(428,5))  &
                     .and. l_apply_diag

  flg_wthetav_flux = sf(293,5) .and. l_apply_diag
  flg_deep_tops    = sf(319,5) .and. l_apply_diag
  flg_mf_deep      = sf(320,5) .and. l_apply_diag
  flg_mf_congest   = sf(321,5) .and. l_apply_diag

  flg_mf_shall     = (sf(322,5) .or. sf(417,5) .or. sf(418,5) .or.             &
                      sf(419,5) .or. sf(420,5)) .and. l_apply_diag
  flg_mf_midlev    = sf(323,5) .and. l_apply_diag
  flg_dt_deep      = sf(324,5) .and. l_apply_diag
  flg_dt_congest   = sf(325,5) .and. l_apply_diag
  flg_dt_shall     = (sf(326,5) .or. sf(409,5) .or. sf(410,5) .or.             &
                      sf(411,5) .or. sf(412,5)) .and. l_apply_diag
  flg_dt_midlev    = sf(327,5) .and. l_apply_diag
  flg_dq_deep      = sf(328,5) .and. l_apply_diag
  flg_dq_congest   = sf(329,5) .and. l_apply_diag
  flg_dq_shall     = (sf(330,5) .or. sf(413,5) .or. sf(414,5) .or.             &
                      sf(415,5) .or. sf(416,5)) .and. l_apply_diag
  flg_dq_midlev    = sf(331,5) .and. l_apply_diag
  flg_du_deep      = sf(332,5) .and. l_apply_diag
  flg_du_congest   = sf(333,5) .and. l_apply_diag
  flg_du_shall     = sf(334,5) .and. l_apply_diag
  flg_du_midlev    = sf(335,5) .and. l_apply_diag
  flg_dv_deep      = sf(336,5) .and. l_apply_diag
  flg_dv_congest   = sf(337,5) .and. l_apply_diag
  flg_dv_shall     = sf(338,5) .and. l_apply_diag
  flg_dv_midlev    = sf(339,5) .and. l_apply_diag

  ! This is also needed for the calculation of an updraught area (5,229)
  flg_w_eqn        = (sf(196,5) .or. sf(197,5) .or. sf(229,5))                 &
                     .and. l_apply_diag

  flg_dt_dd        = sf(198,5) .and. l_apply_diag
  flg_dq_dd        = sf(199,5) .and. l_apply_diag
  flg_du_dd        = sf(175,5) .and. l_apply_diag
  flg_dv_dd        = sf(176,5) .and. l_apply_diag

  flg_area_ud      = sf(229,5) .and. l_apply_diag
  flg_area_dd      = sf(230,5) .and. l_apply_diag

  if ( i_convection_vn == i_cv_comorph .and. l_apply_diag) then
    flg_par_radius_up  = sf(550,5)
    flg_par_radius_dwn = sf(551,5)
    flg_freq_up        = sf(552,5)
    ! Need mass flux
    if (flg_freq_up)  flg_up_flx_half = .true.
    flg_freq_dwn       = sf(553,5)
    if (flg_freq_dwn) flg_dwn_flx_half = .true.

    flg_par_meanup_dtv = sf(554,5)
    flg_par_meandn_dtv = sf(555,5)
    flg_par_meanup_rhl = sf(556,5)
    flg_par_meandn_rhl = sf(557,5)
    flg_par_meanup_q   = sf(558,5)
    flg_par_meandn_q   = sf(559,5)
    flg_par_meanup_qcl = sf(560,5)
    flg_par_meandn_qcl = sf(561,5)
    flg_par_meanup_qcf = sf(562,5)
    flg_par_meandn_qcf = sf(563,5)
    ! Conditional on CoMorph settings
    flg_par_meanup_qrain = sf(564,5) .and. l_cv_rain
    flg_par_meandn_qrain = sf(565,5) .and. l_cv_rain
    flg_par_meanup_qgr   = sf(566,5) .and. l_cv_graup
    flg_par_meandn_qgr   = sf(567,5) .and. l_cv_graup
    flg_par_meanup_qsnow = sf(568,5) .and. l_cv_snow
    flg_par_meandn_qsnow = sf(569,5) .and. l_cv_snow
    flg_par_meanup_cfl   = sf(570,5) .and. l_cv_cloudfrac
    flg_par_meandn_cfl   = sf(571,5) .and. l_cv_cloudfrac
    flg_par_meanup_cff   = sf(572,5) .and. l_cv_cloudfrac
    flg_par_meandn_cff   = sf(573,5) .and. l_cv_cloudfrac
    flg_par_meanup_cfb   = sf(574,5) .and. l_cv_cloudfrac
    flg_par_meandn_cfb   = sf(575,5) .and. l_cv_cloudfrac
    flg_par_meanup_u   = sf(576,5)
    flg_par_meandn_u   = sf(577,5)
    flg_par_meanup_v   = sf(578,5)
    flg_par_meandn_v   = sf(579,5)
    flg_par_meanup_w   = sf(580,5)
    flg_par_meandn_w   = sf(581,5)
    flg_par_meanup_t   = sf(582,5)
    flg_par_meandn_t   = sf(583,5)

    flg_par_coreup_dtv = sf(584,5) .and. l_par_core
    flg_par_coredn_dtv = sf(585,5) .and. l_par_core
    flg_par_coreup_rhl = sf(586,5) .and. l_par_core
    flg_par_coredn_rhl = sf(587,5) .and. l_par_core

    flg_turb_pert_u    = sf(594,5)
    flg_turb_pert_v    = sf(595,5)
    flg_turb_pert_w    = sf(596,5)
    flg_turb_pert_t    = sf(597,5)
    flg_turb_pert_q    = sf(598,5)
    flg_turb_radius    = sf(599,5)

    flg_dry_fraction   = sf(600,5)
    flg_liq_fraction   = sf(601,5)
    flg_mix_fraction   = sf(602,5)
    flg_icr_fraction   = sf(603,5)
    flg_dry_temp       = sf(604,5)
    flg_liq_temp       = sf(605,5)
    flg_mix_temp       = sf(606,5)
    flg_icr_temp       = sf(607,5)
    flg_dry_q_vap      = sf(608,5)
    flg_liq_q_vap      = sf(609,5)
    flg_mix_q_vap      = sf(610,5)
    flg_icr_q_vap      = sf(611,5)
    flg_dry_rhl        = sf(612,5)
    flg_liq_rhl        = sf(613,5)
    flg_mix_rhl        = sf(614,5)
    flg_icr_rhl        = sf(615,5)
  else
    flg_par_radius_up  = .false.
    flg_par_radius_dwn = .false.
    flg_freq_up        = .false.
    flg_freq_dwn       = .false.
    flg_par_meanup_dtv = .false.
    flg_par_meandn_dtv = .false.
    flg_par_meanup_rhl = .false.
    flg_par_meandn_rhl = .false.
    flg_par_meanup_q   = .false.
    flg_par_meandn_q   = .false.
    flg_par_meanup_qcl = .false.
    flg_par_meandn_qcl = .false.
    flg_par_meanup_qcf = .false.
    flg_par_meandn_qcf = .false.
    flg_par_meanup_qrain = .false.
    flg_par_meandn_qrain = .false.
    flg_par_meanup_qgr   = .false.
    flg_par_meandn_qgr   = .false.
    flg_par_meanup_qsnow = .false.
    flg_par_meandn_qsnow = .false.
    flg_par_meanup_cfl   = .false.
    flg_par_meandn_cfl   = .false.
    flg_par_meanup_cff   = .false.
    flg_par_meandn_cff   = .false.
    flg_par_meanup_cfb   = .false.
    flg_par_meandn_cfb   = .false.
    flg_par_meanup_u   = .false.
    flg_par_meandn_u   = .false.
    flg_par_meanup_v   = .false.
    flg_par_meandn_v   = .false.
    flg_par_meanup_w   = .false.
    flg_par_meandn_w   = .false.
    flg_par_meanup_t   = .false.
    flg_par_meandn_t   = .false.

    flg_par_coreup_dtv = .false.
    flg_par_coredn_dtv = .false.
    flg_par_coreup_rhl = .false.
    flg_par_coredn_rhl = .false.

    flg_turb_pert_u    = .false.
    flg_turb_pert_v    = .false.
    flg_turb_pert_w    = .false.
    flg_turb_pert_t    = .false.
    flg_turb_pert_q    = .false.
    flg_turb_radius    = .false.

    flg_dry_fraction   = .false.
    flg_liq_fraction   = .false.
    flg_mix_fraction   = .false.
    flg_icr_fraction   = .false.
    flg_dry_temp       = .false.
    flg_liq_temp       = .false.
    flg_mix_temp       = .false.
    flg_icr_temp       = .false.
    flg_dry_q_vap      = .false.
    flg_liq_q_vap      = .false.
    flg_mix_q_vap      = .false.
    flg_icr_q_vap      = .false.
    flg_dry_rhl        = .false.
    flg_liq_rhl        = .false.
    flg_mix_rhl        = .false.
    flg_icr_rhl        = .false.
  end if

else

  ! Setting of flags defined in module

  l_apply_diag     = .true.

  flg_conv_rain_3d = .true.
  flg_conv_snow_3d = .true.

  flg_up_flx       = .true.
  flg_up_flx_half  = .true.
  flg_dwn_flx      = .true.
  flg_dwn_flx_half = .true.
  flg_entr_up      = .true.
  flg_detr_up      = .true.
  flg_entr_dwn     = .true.
  flg_detr_dwn     = .true.

  flg_uw_dp        = .true.
  flg_vw_dp        = .true.
  flg_uw_shall     = .true.
  flg_vw_shall     = .true.
  flg_uw_mid       = .true.
  flg_vw_mid       = .true.

  l_qcl_incr_cinh  = .true.
  l_qcf_incr_cinh  = .true.
  l_bcf_incr_cinh  = .true.
  l_cfl_incr_cinh  = .true.
  l_cff_incr_cinh  = .true.
  l_bcf_incr_conv  = .true.
  l_cfl_incr_conv  = .true.
  l_cff_incr_conv  = .true.

  l_theta_incr_conv= .true.
  l_T_incr_conv    = .true.
  l_q_incr_conv    = .true.
  l_qcl_incr_conv  = .true.
  l_qcf_incr_conv  = .true.
  if ( i_convection_vn == i_cv_comorph ) then
    ! CoMorph convection scheme produces increments to any optional
    ! condensed water fields that are in use:
    l_qcf2_incr_conv   = l_mcr_qcf2
    l_qgraup_incr_conv = ( graupel_option > no_graupel )
    l_qrain_incr_conv  = l_mcr_qrain
  else
    ! Other convection schemes do not calculate increments for these fields:
    l_qcf2_incr_conv   = .false.
    l_qgraup_incr_conv = .false.
    l_qrain_incr_conv  = .false.
  end if

  l_t_conv_only = .false.         ! Not output by the SCM for some reason.
  l_q_conv_only = .false.

  flg_dt_dd = .true.
  flg_dq_dd = .true.

  ! 5A & 6A convection schemes
  flg_du_dd    = .true.
  flg_dv_dd    = .true.
  flg_area_ud  = .true.
  flg_area_dd  = .true.
  flg_wqt_flux = .true.
  flg_wql_flux = .true.

  flg_wthetal_flux = .true.
  flg_wthetav_flux = .true.
  flg_deep_tops    = .true.

  flg_mf_deep    = .true.
  flg_mf_congest = .true.
  flg_mf_shall   = .true.
  flg_mf_midlev  = .true.
  flg_dt_deep    = .true.
  flg_dt_congest = .true.
  flg_dt_shall   = .true.
  flg_dt_midlev  = .true.
  flg_dq_deep    = .true.
  flg_dq_congest = .true.
  flg_dq_shall   = .true.
  flg_dq_midlev  = .true.
  flg_du_deep    = .true.
  flg_du_congest = .true.
  flg_du_shall   = .true.
  flg_du_midlev  = .true.
  flg_dv_deep    = .true.
  flg_dv_congest = .true.
  flg_dv_shall   = .true.
  flg_dv_midlev  = .true.
  flg_w_eqn      = .true.

  ! Comorph convection scheme
  if ( i_convection_vn == i_cv_comorph ) then
    flg_par_radius_up  = .true.
    flg_par_radius_dwn = .true.
  else
    flg_par_radius_up  = .false.
    flg_par_radius_dwn = .false.
  end if
  flg_freq_up        = .false. ! not really needed for SCM as not potentially
  flg_freq_dwn       = .false. ! time averaging diagnostics
  ! Following come out in comorph_diags_scm_mod directly from CoMorph's
  ! diagnostic arrays so don't need space allocated on a flag
  flg_par_meanup_dtv = .false.
  flg_par_meandn_dtv = .false.
  flg_par_meanup_rhl = .false.
  flg_par_meandn_rhl = .false.
  flg_par_meanup_q   = .false.
  flg_par_meandn_q   = .false.
  flg_par_meanup_qcl = .false.
  flg_par_meandn_qcl = .false.
  flg_par_meanup_qcf = .false.
  flg_par_meandn_qcf = .false.
  flg_par_meanup_qrain = .false.
  flg_par_meandn_qrain = .false.
  flg_par_meanup_qgr   = .false.
  flg_par_meandn_qgr   = .false.
  flg_par_meanup_qsnow = .false.
  flg_par_meandn_qsnow = .false.
  flg_par_meanup_cfl   = .false.
  flg_par_meandn_cfl   = .false.
  flg_par_meanup_cff   = .false.
  flg_par_meandn_cff   = .false.
  flg_par_meanup_cfb   = .false.
  flg_par_meandn_cfb   = .false.
  flg_par_meanup_u   = .false.
  flg_par_meandn_u   = .false.
  flg_par_meanup_v   = .false.
  flg_par_meandn_v   = .false.
  flg_par_meanup_w   = .false.
  flg_par_meandn_w   = .false.
  flg_par_meanup_t   = .false.
  flg_par_meandn_t   = .false.

  flg_par_coreup_dtv = .false.
  flg_par_coredn_dtv = .false.
  flg_par_coreup_rhl = .false.
  flg_par_coredn_rhl = .false.

  flg_turb_pert_u    = .false.
  flg_turb_pert_v    = .false.
  flg_turb_pert_w    = .false.
  flg_turb_pert_t    = .false.
  flg_turb_pert_q    = .false.
  flg_turb_radius    = .false.

  flg_dry_fraction   = .false.
  flg_liq_fraction   = .false.
  flg_mix_fraction   = .false.
  flg_icr_fraction   = .false.
  flg_dry_temp       = .false.
  flg_liq_temp       = .false.
  flg_mix_temp       = .false.
  flg_icr_temp       = .false.
  flg_dry_q_vap      = .false.
  flg_liq_q_vap      = .false.
  flg_mix_q_vap      = .false.
  flg_icr_q_vap      = .false.
  flg_dry_rhl        = .false.
  flg_liq_rhl        = .false.
  flg_mix_rhl        = .false.
  flg_icr_rhl        = .false.

end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine set_convection_output_flags

end module cv_stash_flg_mod
