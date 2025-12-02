! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: science_fixes_mod
!
! This module declares 'short-term' temporary logicals used to protect
! science bug fixes that lead to significant alterations in science results.
! It is expected that these logicals will be short lived as the preference
! should be for all configurations to use the corrected code. But
! to maintain short term reproducibility of results across UM versions
! the fixes are protected by logicals until the fixes become the default
! in all model configurations and the logical is retired.
!
! All logicals below should have a review period attached to them for future
! retirement of both the logical and the broken code.
!
! ! ticket #xxxx
! logical :: fix_me = .false.   ! review again MM YYYY
!
! then add logical to namelist /temp_fixes/
! and add code to subroutine warn_temp_fixes to report when fix is not used
! ie when the logical = .false.
!
! -----------------------------------------------------------
! -----------------------------------------------------------

module science_fixes_mod

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
use errormessagelength_mod, only: errormessagelength
use umprintmgr, only: newline
use model_domain_mod, only: model_type, mt_lam
use missing_data_mod, only: imdi

implicit none

! ticket  5456
logical :: l_fix_conserv = .false.      ! Review Again April 2015
! ticket #5841
logical :: l_fix_arcl_eg_levs = .false. ! Review Again April 2015
! ticket #6514
logical :: l_pc2_homog_turb_q_neg = .false. ! Review again Sept 2015
! ticket #6357
logical :: l_methox_fix = .false.       ! Review Again January 2016

! Below this point tickets correspond to the Science Repository Service
!----------------------------------------------------------------------
! JULES science fixes are now contained in jules_science_fixes_mod within JULES
!******************************************************************************
! ticket #47
logical :: l_stph_rhcrit_unbias = .false.     ! Review Again January 2016
! ticket #911
logical :: l_eg_damp_height_lid = .false.     ! Review Again Oct 2016
! ticket #430
logical :: l_rm_hardwire_gas360 = .false.   ! Review again Sep 2016
! ticket #1167
logical :: l_fix_conv_precip_evap = .false.   ! Review Again April 2019
! ticket #1421
logical :: l_fix_ukca_impscav = .false.   ! Review again Jan 2018
! ticket #1638
logical :: l_fix_rp_shock_amp = .false.   ! Review again Jan 2018
! ticket #2077
logical :: l_fix_dyndiag = .false.        ! Review again Sep 2018
! ticket #2556
logical :: l_fix_riming = .false.         ! Review in Dec 2018
! ticket #3076
logical :: l_fix_ccb_cct = .false.        ! Review in July 2017
! ticket #3005
logical :: l_fix_zh = .false.             ! Review in Dec 2018
! ticket #2405
logical :: l_fix_nacl_density = .false.   ! Review in Dec 2018
! ticket #3011
logical :: l_fix_rcf_mlsnow_icefreemax = .false. ! Review in Dec 2020
! ticket #3681
logical :: l_fix_conv_diags_var = .false. ! Review in Jan 2020
! ticket #2545
logical :: l_fix_lsp_incs_to_spt = .false. ! Review in Feb 2019
! ticket #2070
logical :: l_fix_ec_gen_hgt = .false.     ! Review in Dec 2018
! ticket #1250 (review in ticket #3997)
logical :: l_fix_improve_drydep = .false. ! Review in Jan 2019
! ticket #5651 (needed for UKESM1.0)
logical :: l_fix_drydep_so2_water = .false. ! Review in Jan 2023
! ticket #4157 (review in ticket #4245)
logical :: l_fix_ukca_h2dd_x = .false.    ! Review in May 2019
! ticket #4383 (review in ticket #4416)
logical :: l_fix_neg_pvol_wat = .false.    ! Review in Jan 2020
! ticket #3245
logical :: l_fix_ukca_h2so4_ystore = .false. ! Review Oct 2020
! ticket #4501
logical :: l_fix_pc2_cnv_mix_phase = .false. ! Review Dec 2021
! ticket #5031
logical :: l_fix_tidy_rainfracs = .false. ! Review Jan 2021
! ticket #5080
logical :: l_fix_incloud_qcf = .false. ! Review Dec 2021
logical :: l_fix_mcr_frac_ice = .false. ! Review Dec 2021
! ticket #5118
logical :: l_fix_ukca_offox_h2o_fac = .false. ! Review Jan 2021
! ticket #5009
logical :: l_fix_ukca_cloud_frac = .false. ! Review Jan 2021
! ticket #5233
logical :: l_fix_ukca_activate_vert_rep  = .false.
! Ticket #5479
logical :: l_fix_gr_autoc = .false. ! Review Apr 2022
! ticket #5532
logical :: l_fix_ukca_activate_pdf = .false. ! Review Dec 2021
! ticket #5497
logical :: l_fix_true_latlon = .false. ! Review April 2021
! ticket #6210
logical :: l_use_q1p5m_in_cape_diag = .false. ! Review Mar 2023
! ticket #5472 and associated parameters
integer :: i_fix_mphys_drop_settle = imdi ! Review September 2021
integer, parameter :: no_fix     = 0
integer, parameter :: first_fix  = 1
integer, parameter :: second_fix = 2
! ticket #6091
logical :: l_pc2_checks_sdfix = .false. ! Review September 2021
! ticket #6252
logical :: l_improve_cv_cons = .false.      ! Review Mar 2023
! ticket #6088
logical :: l_improve_aero_drydep = .false.  ! Review Dec 2022
! ticket #6438
logical :: l_enforce_f03_compliance = .true. ! Review July 2023
! ticket #6174
logical :: l_fix_ukca_hygroscopicities = .false.
! ticket #7634
logical :: l_fix_ukca_water_content = .false.
! ticket #7232
logical :: l_fix_level_indexing_bimodal = .false. ! Review October 2024
!------------------------------------------------

namelist/temp_fixes/                                                           &
        l_fix_conserv, l_fix_arcl_eg_levs,                                     &
        l_pc2_homog_turb_q_neg,                                                &
        l_methox_fix, l_stph_rhcrit_unbias,                                    &
        l_eg_damp_height_lid,                                                  &
        l_rm_hardwire_gas360,                                                  &
        l_fix_conv_precip_evap, l_fix_ukca_impscav,                            &
        l_fix_rp_shock_amp, l_fix_dyndiag,                                     &
        l_fix_riming, l_fix_zh,                                                &
        l_fix_ccb_cct, l_fix_nacl_density,                                     &
        l_fix_rcf_mlsnow_icefreemax, l_fix_conv_diags_var,                     &
        l_fix_lsp_incs_to_spt, l_fix_ec_gen_hgt,                               &
        l_fix_improve_drydep, l_fix_drydep_so2_water,                          &
        l_fix_ukca_h2dd_x, l_fix_neg_pvol_wat,                                 &
        l_fix_pc2_cnv_mix_phase,                                               &
        l_fix_ukca_h2so4_ystore, l_fix_tidy_rainfracs,                         &
        l_fix_incloud_qcf, l_fix_mcr_frac_ice,                                 &
        l_fix_ukca_offox_h2o_fac, l_fix_ukca_cloud_frac, l_fix_gr_autoc,       &
        l_fix_ukca_activate_vert_rep, l_fix_ukca_activate_pdf,                 &
        l_fix_true_latlon, l_use_q1p5m_in_cape_diag, l_pc2_checks_sdfix,       &
        i_fix_mphys_drop_settle, l_improve_cv_cons, l_improve_aero_drydep,     &
        l_enforce_f03_compliance, l_fix_ukca_hygroscopicities,                 &
        l_fix_ukca_water_content, l_fix_level_indexing_bimodal

! -----------------------------------------------------------
! -----------------------------------------------------------
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='SCIENCE_FIXES_MOD'

contains

subroutine warn_temp_fixes()

use ereport_mod, only: ereport

implicit none

real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WARN_TEMP_FIXES'

integer :: ErrorStatus            ! Return code : 0 Normal Exit : >0 Error
character(len=errormessagelength) :: cmessage
                                  ! Error message if Errorstatus /=0

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0


! -----------------------------------------------------------
! -----------------------------------------------------------
! define whether the fix is appropriate to RCF, UM RUN or both
! and warn the user if the fix is not used.


#if defined(RECON)

if (.not. l_fix_ec_gen_hgt) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Recon run excludes ticket #2070 as l_fix_ec_gen_hgt=.false.'//     newline//&
  'Level height generation for ECMWF data will not be accurate.'

  call ereport(RoutineName, ErrorStatus, cmessage)
end if

#else

if ( .not. l_fix_conserv ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #5456 as l_fixed_conserv=.false.'//      newline//&
  'This affects the accuracy of tracers conservation'

  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if (.not. l_fix_arcl_eg_levs) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #5841 as l_fix_arcl_eg_levs=.false.'//   newline//&
  'This will affect any ENDGame model runs which use aerosol climatologies.'

  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if ( i_fix_mphys_drop_settle /= second_fix) then

  if ( i_fix_mphys_drop_settle == first_fix ) then
    ErrorStatus = -100
    cmessage    =                                                     newline//&
    'Model run is using the older fix of droplet settling. Please  '//newline//&
    'migrate to the newer fix by setting i_fix_mphys_drop_settle=2 '//newline//&
    '(second fix) in the namelist temp_fixes'

  else if ( i_fix_mphys_drop_settle == no_fix ) then
    ErrorStatus = -100
    cmessage    = 'Model run excludes ticket #5472 as '//             newline//&
    'i_fix_mphys_drop_settle == 0. The recommended fix is to set   '//newline//&
    'i_fix_mphys_drop_settle to 2 (second fix) in the namelist temp_fixes'

    call ereport(RoutineName, ErrorStatus, cmessage)

  else
    ! Default case for e.g. the user tweaking the namelist by hand.
    ErrorStatus = 100 ! cause run to abort
    cmessage = 'Integer variable i_fix_mphys_drop_settle has not   '//newline//&
    'been set in the namelist temp_fixes.'

    call ereport(RoutineName, ErrorStatus, cmessage)

  end if ! i_fix_mphys_drop_settle == first_fix
end if ! i_fix_mphys_drop_settle /= second_fix


if (.not. l_pc2_homog_turb_q_neg) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #6514 as l_pc2_homog_turb_q_neg=.false.'//        &
                                                                      newline//&
  'your results might contain negative humidities especially if'//    newline//&
  'forced_cu setting is >0.'

  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if (.not. l_methox_fix) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #6357 as l_methox_fix=.false.'//         newline//&
  'This will affect any model runs which use the basic'//             newline//&
  'methane oxidation scheme.'

  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if (.not. l_stph_rhcrit_unbias) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #47 as l_stph_rhcrit_unbias=.false.'//   newline//&
  'If using Random Parameters, the RHCRIT perturbations can cause'//  newline//&
  'cloud cover biases.'

  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if (.not. l_eg_damp_height_lid) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #911 as l_eg_damp_height_lid=.false.'//  newline//&
  'The user is required to manually input the reference height' //    newline//&
  '(via the variable damp_height) rather than use the model top.' //  newline//&
  'This will affect any ENDGame model run.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_rm_hardwire_gas360) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #430 as l_rm_hardwire_gas360=.false.'//  newline//&
  'This will hardwire trace gas interpolation to use 360day'//        newline//&
  'calendar even if gregorian calendar set'

  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_conv_precip_evap) then
  ErrorStatus = -100
  CMessage    = 'Model run excludes ticket #1167 as'               // newline//&
                ' l_fix_conv_precip_evap=.false.'                  // newline//&
                ' This will affect all runs with parameterised'    // newline//&
                ' precipitating convection; evaporation rates for' // newline//&
                ' convective precipitation will be underestimated.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_ukca_impscav) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #1421 as l_fix_ukca_impscav=.false.'//   newline//&
  'This will affect any model runs which include UKCA.'

  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_rp_shock_amp) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #1638 as l_fix_rp_shock_amp=.false.'//   newline//&
  'This will affect any model using the Random Parameters 2b Scheme'//newline//&
  'by doubling the size of the shock amplitude.'

  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_dyndiag) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #2077 as l_fix_dyndiag=.false.'//        newline//&
  'This will mess up shear-dominated PBLs slightly'

  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_riming ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #2556 as'//                newline//&
  ' l_fix_riming=.false.'//                                           newline//&
  'This could affect any model runs where l_shape_rime is .true.   '//newline//&
  'and may produce small numbers in the microphysics leading to    '//newline//&
  'the UM crashing in the solver.'

  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_zh) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #3005 as l_fix_zh=.false.'//             newline//&
  'This will particularly affect the boundary layer depth  '//        newline//&
  'diagnostic (25 or 3025) in model runs which have non-zero '//      newline//&
  'setting for forced_cu but also have a small impact on '//          newline//&
  'model evolution in all runs.'


  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_ccb_cct ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #3076 as'//                newline//&
  ' l_fix_ccb_cct=.false.  '//                                        newline//&
  'This affects any model runs using 5A or 6A convection schemes.  '//newline//&
  'Without this switch, the convection scheme may occasionally '//    newline//&
  'diagnose cloud-top and cloud-base levels inconsistently.  '//      newline//&
  'In the 5A scheme it may spuriously leave the cloud-top level '//   newline//&
  'unset even where it has set the cloud-base, leading to '//         newline//&
  'inconsistent ccrad fields.'

  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_nacl_density ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #2405 as'//                newline//&
  ' l_fix_nacl_density=.false.'//                                     newline//&
  'This affects the sea-salt emissions generated by the model for  '//newline//&
  'the UKCA + GLOMAP aerosol model'

  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_rcf_mlsnow_icefreemax ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #3011 as'//                newline//&
  ' l_fix_rcf_mlsnow_icefreemax=.false.'//                            newline//&
  'This affects the cap imposed on the mass of snow at ice points '// newline//&
  'reconfigured to ice-free points in the multilayer snow scheme. '

  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_conv_diags_var ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #3681 as'//                newline//&
  ' l_fix_conv_diags_var=.false.'//                                   newline//&
  'This affects the diagnostics that are passed from the'//           newline//&
  'convection scheme to the PF model. Without this switch'//          newline//&
  'the increments will incorrectly include contributions'//           newline//&
  'from cloud erosion, a temperature rather than potential'//         newline//&
  'temperature will be passed, and the mass flux will exclude'//      newline//&
  'points where it is reducing with height.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_lsp_incs_to_spt) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #2545 as l_fix_lsp_incs_to_spt=.false.'//newline//&
  'This will mean that the microphysics changes due to mixed phase'// newline//&
  'turbulence will not be seen in the SPT scheme.'
  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if (.not. l_fix_improve_drydep) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #1250 as'//                newline//&
  ' l_fix_improve_drydep=.false.'//                                   newline//&
  ' This will mean that dry deposition velocities are set to null'//  newline//&
  ' for HCl, HOCl, HBr, HOBr, H2SO4, MeOH and Sec_Org and that dry'// newline//&
  ' deposition velocities for 9 tiles are inconsistant with 13/17/27 tiles.'
  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if (.not. l_fix_drydep_so2_water) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #5651 as'//                newline//&
  ' l_fix_drydep_so2_water=.false.'//                                 newline//&
  ' This will mean that the surface resistance of water when '//      newline//&
  ' calculating the dry deposition of SO2 is too high for 13, 17 '//  newline//&
  ' and 27 tiles.'
  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if (.not.  l_fix_ukca_h2dd_x) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #4157 as'//                newline//&
  ' l_fix_ukca_h2dd_x=.false. .'//                                    newline//&
  ' This will mean that the wrong element is used for h2dd_c & h2dd_m.'
  call ereport(RoutineName, ErrorStatus, cmessage)
end if

if (.not. l_fix_neg_pvol_wat ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #4157 as'//                newline//&
  ' l_fix_neg_pvol_wat=.false. .'//                                   newline//&
  ' This will mean that field pvol_wat can contain a negative concentration.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not.  l_fix_ukca_h2so4_ystore) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #3245 as'//                newline//&
  ' l_fix_ukca_h2so4_ystore=.false. .'//                              newline//&
  ' This will mean that the UKCA H2SO4 tracer will not be updated correctly.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_pc2_cnv_mix_phase ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #4501 as'//                newline//&
  ' l_fix_pc2_cnv_mix_phase=.false.'//                                newline//&
  'Without this switch water will not be exactly conserved'//         newline//&
  'in PC2s liquid and ice water partition calculation'//              newline//&
  'within the 6a convection scheme.'//                                newline//&
  'Note that if l_cv_conserve_check=.true. a retrospective '//        newline//&
  'correction will enforce water conservation in convection.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_tidy_rainfracs ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #5031 as'//                newline//&
  ' l_fix_tidy_rainfracs=.false.'//                                   newline//&
  'This means that if l_mcr_qcf2 =.false. the rain fractions'//       newline//&
  'will not be updated from the melting of snow in lsp_tidy.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if ( .not. l_fix_incloud_qcf ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #5080 as'//                newline//&
  ' l_fix_incloud_qcf=.false. .'//                                    newline//&
  ' When true this fixes a bug in the PC2 code after inhomogeneous'// newline//&
  ' forcing of ice cloud by turbulence which checks to avoid'//       newline//&
  ' generating unreasonably high in-cloud ice water contents'//       newline//&
  ' qcf/cff.  Such occurences wrongly result in the ice-cloud'//      newline//&
  ' fraction being reset to near-zero instead of increased'//         newline//&
  ' as intended.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if ( .not. l_fix_mcr_frac_ice ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #5080 as'//                newline//&
  ' l_fix_mcr_frac_ice=.false. .'//                                   newline//&
  ' This switch fixes 2 minor issues in the treatment of ice-cloud'// newline//&
  ' fractions by the microphysics scheme, both of which can'//        newline//&
  ' result in the ice fraction cff going to near-zero while'//        newline//&
  ' retaining a significant grid-mean ice mass qcf.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_ukca_offox_h2o_fac ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #5118 as'//                newline//&
  ' l_fix_ukca_offox_h2o_fac=.false.'//                               newline//&
  'This means that water vapour will not be in vmr when'//            newline//&
  'using UKCA Offline-Oxidants.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_ukca_cloud_frac ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #5009 as'//                newline//&
  'l_fix_ukca_cloud_frac=.false.'//                                   newline//&
  'This means that values of the cloud fraction field used in'//      newline//&
  'aerosol chemistry will be those from the level below the'//        newline//&
  'correct level due to an array indexing bug.'//                     newline//&
  'Runs which do not use UKCA aerosol chemistry are unaffected.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_gr_autoc ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes graupel autoconversion fix from ticket'//       newline//&
  '#5479. This produces spurious graupel in sublimating ice'//        newline//&
  'clouds in any configurations which include prognostic graupel.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_ukca_activate_vert_rep ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #5233 as'//                newline//&
  'l_fix_ukca_activate_vert_rep=.false.'//                            newline//&
  'This means that droplets in contiguous vertical columns'//         newline//&
  'will be replicated in per cm^-3 instead of per kg with the'//      newline//&
  'latter being more physical. The cldflag'//                         newline//&
  'will also be incorrectly used a real'//                            newline//&
  'instead of an integer.'//                                          newline//&
  'Runs which do not use UKCA aerosol chemistry are unaffected.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_ukca_activate_pdf ) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #5532 as'//                newline//&
  'l_fix_ukca_activate_pdf=.false.'//                                 newline//&
  'This means that if the grid-cell mean updraft speed exceeds'//     newline//&
  'approx twice the width of the updraft PDF (defined by the  '//     newline//&
  'BL TKE) the full PDF will not get integrated over and the  '//     newline//&
  'cloud droplet number concentration will be wrong. Runs     '//     newline//&
  'which do not use UKCA or GLOMAP aerosol climatologies, and '//     newline//&
  'runs which do not use parameterized convection, are unaffected.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_true_latlon ) then
  ErrorStatus = -100
  cmessage    = 'Model run excludes a change from ticket #5497 as'//  newline//&
  'l_fix_true_latlon = .false.'//                                     newline//&
  'This means errors in the values of these arrays carried over'//    newline//&
  'from the New Dynamics will be present in the run.'//               newline//&
  'These errors effect the radiation scheme and B-grid diagnostics.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_use_q1p5m_in_cape_diag) then
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run includes temporary change from ticket #6120 as'//        newline//&
  'l_use_q1p5m_in_cape_diag = .false.'//                              newline//&
  'PWS convection diagnostics will use level 1 water vapour'//        newline//&
  'mixing ratio instead of screen-level humidity in surface-based'//  newline//&
  'parcel ascents.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_pc2_checks_sdfix ) then
  ErrorStatus = -100
  cmessage    = 'Model run excludes a change from ticket #6091 as'//  newline//&
  'l_pc2_checks_sdfix = .false.'//                                    newline//&
  'This means instances of total liquid-cloud cover in'//             newline//&
  'sub-saturated grid-boxes (which should be impossible) are left'//  newline//&
  'uncorrected by the pc2_checks routine if there is not enough'//    newline//&
  'cloud liquid water content to bring the grid-box to saturation'//  newline//&
  'by evaporating it.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_improve_cv_cons ) then
  ErrorStatus = -100
  cmessage    = 'Model run excludes a change from ticket #6252 as'//  newline//&
  'l_improve_cv_cons = .false.'//                                     newline//&
  'This means there will be errors in conservation in the'//          newline//&
  'convection scheme because of an incorrect conversion from'//       newline//&
  'specific humidities to mixing ratios'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_improve_aero_drydep ) then
  ErrorStatus = -100
  cmessage    = 'Model run excludes a change from ticket #6088 as'//  newline//&
  'l_improve_aero_drydep = .false.'//                                 newline//&
  'This means the incorrect roughness length is being used in the'//  newline//&
  'UKCA aerosol dry deposition, as well as a possible mismatch in'//  newline//&
  ' surface types seen by the scheme'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_ukca_hygroscopicities) then
  ErrorStatus = -100
  cmessage    = 'Model run excludes a change from ticket #6174 as'//  newline//&
  'l_fix_ukca_hygroscopicities = .false.'//                           newline//&
  'This means hygroscopicities used in aerosol activation will be'//  newline//&
  'outdated and, if UKCA-ACTIVATE is used, incorrectly averaged'//    newline//&
  'across the different composition types if more than one'//         newline//&
  'soluble species is used in a given mode.'//                        newline//&
  'These errors affect both UKCA and GLOMAP-CLIM runs.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_ukca_water_content) then
  ErrorStatus = -100
  cmessage    = 'Model run excludes a change from ticket #7634 as'//  newline//&
  'l_fix_ukca_water_content = .false.'//                              newline//&
  'This means that aerosol water content values will be '//           newline//&
  'incorrect.'//                                                      newline//&
  'These errors affect both UKCA and GLOMAP-CLIM runs.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

#endif

if (.not. l_enforce_f03_compliance ) then
  ErrorStatus = -100
  cmessage    = 'Run excludes a change from ticket #6438 as'//        newline//&
  'l_enforce_f03_compliance = .false.'//                              newline//&
  'This means that ancilliary file(s) identified as containing '//    newline//&
  'time mean fields which are not compliant with UMDP F3 will be '//  newline//&
  'permitted.  Interpolation to the correct time value cannot be '//  newline//&
  'guaranteed.'
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

if (.not. l_fix_level_indexing_bimodal) then
  ErrorStatus = -100
  cmessage    = 'Run excludes a change from ticket #7232 as'//        newline//&
  'l_fix_level_indexing_bimodal = .false.'//                          newline//&
  'This means that the correct indexing supplying pressure    '//     newline//&
  'levels to the bimodal cloud scheme will not be applied    '//      newline//&
  'meaning that the lowest level will be duplicated in error.    '//  newline//&
  '           '
  call ereport(RoutineName, ErrorStatus, CMessage)
end if

! -----------------------------------------------------------
! -----------------------------------------------------------


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine warn_temp_fixes


subroutine print_nlist_temp_fixes()

use umPrintMgr, only: umPrint, maxLineLen

implicit none
character(len=maxLineLen) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_TEMP_FIXES'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist temp_fixes', src=ModuleName)

write(lineBuffer,'(A,L1)') ' l_fix_conserv = ', l_fix_conserv
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_arcl_eg_levs = ', l_fix_arcl_eg_levs
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_pc2_homog_turb_q_neg = ',l_pc2_homog_turb_q_neg
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_methox_fix = ',l_methox_fix
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_stph_rhcrit_unbias = ', l_stph_rhcrit_unbias
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_eg_damp_height_lid = ',l_eg_damp_height_lid
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_rm_hardwire_gas360 = ',l_rm_hardwire_gas360
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_conv_precip_evap = ',l_fix_conv_precip_evap
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ukca_impscav = ',l_fix_ukca_impscav
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_rp_shock_amp = ',l_fix_rp_shock_amp
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_dyndiag = ',l_fix_dyndiag
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_riming = ',l_fix_riming
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_zh = ',l_fix_zh
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ccb_cct = ',l_fix_ccb_cct
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_nacl_density = ',l_fix_nacl_density
call umPrint(lineBuffer,src=ModuleName)
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_rcf_mlsnow_icefreemax = ',                  &
                             l_fix_rcf_mlsnow_icefreemax
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_conv_diags_var = ',l_fix_conv_diags_var
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_lsp_incs_to_spt = ',l_fix_lsp_incs_to_spt
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ec_gen_hgt = ',l_fix_ec_gen_hgt
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_improve_drydep', l_fix_improve_drydep
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_drydep_so2_water', l_fix_drydep_so2_water
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') 'l_fix_ukca_h2dd_x', l_fix_ukca_h2dd_x
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') 'l_fix_neg_pvol_wat', l_fix_neg_pvol_wat
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') 'l_fix_ukca_h2so4_ystore',                          &
                            l_fix_ukca_h2so4_ystore
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_pc2_cnv_mix_phase = ',                      &
                             l_fix_pc2_cnv_mix_phase
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_tidy_rainfracs = ', l_fix_tidy_rainfracs
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_incloud_qcf = ', l_fix_incloud_qcf
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_mcr_frac_ice = ', l_fix_mcr_frac_ice
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ukca_offox_h2o_fac = ',                     &
                             l_fix_ukca_offox_h2o_fac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ukca_cloud_frac = ', l_fix_ukca_cloud_frac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_gr_autoc = ', l_fix_gr_autoc
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ukca_activate_vert_rep = ',                 &
  l_fix_ukca_activate_vert_rep
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ukca_activate_pdf = ',                      &
                             l_fix_ukca_activate_pdf
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_true_latlon = ', l_fix_true_latlon
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_use_q1p5m_in_cape_diag = ',                     &
                             l_use_q1p5m_in_cape_diag
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,I0)') ' i_fix_mphys_drop_settle =',                       &
                             i_fix_mphys_drop_settle
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_pc2_checks_sdfix = ', l_pc2_checks_sdfix
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_improve_cv_cons = ', l_improve_cv_cons
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_improve_aero_drydep = ', l_improve_aero_drydep
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_enforce_f03_compliance = ',                     &
                             l_enforce_f03_compliance
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ukca_hygroscopicities = ',                  &
                             l_fix_ukca_hygroscopicities
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_ukca_water_content = ',                     &
                             l_fix_ukca_water_content
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,L1)') ' l_fix_level_indexing_bimodal =',                  &
                             l_fix_level_indexing_bimodal
call umPrint(lineBuffer,src=ModuleName)

call umPrint('- - - - - - end of namelist - - - - - -', src=ModuleName)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_temp_fixes


end module science_fixes_mod
