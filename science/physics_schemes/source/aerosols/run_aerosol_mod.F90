! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing tunable parameters and dependent logicals
!   used in the aerosol scheme
!
module run_aerosol_mod

!
! Description:
!   This module contains declarations for tunable parameters
!   used to set the emission heights/levels of manmade aerosols
!   and dependent Classic aerosol logicals
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: aerosols
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards UMDP 003
!
use missing_data_mod,      only: imdi
use ereport_mod,           only: ereport
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength

implicit none

!
! Control logicals for Classic aerosol options
! ============================================
!
! Sulphur Cycle scheme
logical :: l_sulpc_so2 = .false.  ! S Cycle: SO2 MMR included
! with emission options
logical :: l_so2_surfem = .false. ! SO2 Surface Emissions
logical :: l_so2_hilem = .false.  ! SO2 High Level Emissions
integer :: SO2_high_level = imdi   ! for chimney SO2 emission
logical :: l_so2_natem = .false.  ! SO2 Natural Emissions
! Dimethyl sulphide
logical :: l_sulpc_dms = .false.  ! S Cycle: DMS MMR included
logical :: l_dms_em = .false.     ! DMS Emissions
logical :: l_dms_em_inter = .false. ! Interactive DMS Emissions
! Switch to choose scheme to use for interactive sea-air exchange of DMS
integer :: i_dms_flux = imdi
! Ozone included for oxidation of DMS and SO2
logical :: l_sulpc_ozone = .false. ! S Cycle: Ozone oxidation
! Online oxidants from UKCA
logical :: l_sulpc_online_oxidants = .false. ! S Cycle : Online oxidants
! Depleted oxidants are passed back to UKCA (.false. is deemed unsafe)
logical :: l_sulpc_2_way_coupling = .true.  ! S Cycle : Depleted oxidants
! SO2+O3 reaction not buffered by NH3
logical :: l_sulpc_so2_o3_nonbuffered = .false. ! S Cycle: SO2+O3 reaction
! Ammonia
logical :: l_sulpc_nh3 = .false.  ! S Cycle: NH3 tracer included
logical :: l_nh3_em = .false.     ! S Cycle: NH3 emissions

! Soot scheme
logical :: l_soot = .false.  ! Soot scheme
! with emission options
logical :: l_soot_surem = .false. ! Soot Surface Emissions
logical :: l_soot_hilem = .false.  ! Soot High Level Emissions

! Biomass scheme
logical :: l_biomass = .false.  ! Biomass scheme
! with emission options
logical :: l_bmass_surem = .false. ! Biomass Surface Emissions
logical :: l_bmass_hilem = .false.  ! Biomass High Level Emissions

! OCFF (Organic Carbon from Fossil Fuels) scheme
logical :: l_ocff = .false.  ! OCFF scheme
! with emission options
logical :: l_ocff_surem = .false.  ! OCFF Surface Emissions
logical :: l_ocff_hilem = .false.  ! OCFF High Level Emissions

! Ammonium nitrate aerosol
logical :: l_nitrate = .false.     ! Ammonium nitrate aerosol

! Additional Constraints

! Use sulphate aerosol no. in S-cycle
logical :: l_use_sulphate_sulpc = .false.

! Use biomass aerosol no. in S-cycle
logical :: l_use_bmass_sulpc = .false.

! Use organic carbon fossil fuel aerosol no. in S-cycle
logical :: l_use_ocff_sulpc = .false.

! Use ammonium nitrate aerosol no. in S-cycle
logical :: l_use_nitrate_sulpc = .false.

! Use sea-salt aerosol no. in S-cycle
logical :: l_use_seasalt_sulpc = .false.

! Use sea-salt aerosol no. in PM diagnostics (unconditional choice)
logical :: l_use_seasalt_pm = .false.

! Use diurnal and weekly emission cycles for primary emissions in CLASSIC
logical :: l_temporal_emi = .false.

! Use variable injection heights for biomass burning emissions
logical :: l_bmass_hilem_variable = .false.

!
! Control logicals for Classic aerosol LBC options
! ================================================
!
! Sulphur Cycle scheme
logical :: l_so2_lbc = .false.  ! S Cycle LBCs included
! Dimethyl sulphide
logical :: l_dms_lbc = .false.  ! S Cycle: DMS LBC included
! Ammonia
logical :: l_nh3_lbc = .false.  ! S Cycle: NH3 LBC included
! Soot scheme
logical :: l_soot_lbc = .false.  ! Soot scheme LBCs included
! Biomass scheme
logical :: l_bmass_lbc = .false.  ! Biomass scheme LBCs included
! OCFF (Organic Carbon from Fossil Fuels) scheme
logical :: l_ocff_lbc = .false.  ! OCFF scheme LBCs included
! Ammonium nitrate aerosol
logical :: l_nitr_lbc = .false.     ! Nitrate scheme LBCs included

!
! Dependent logicals for Classic aerosol options
! ==============================================
!
logical :: l_so2 = .false.
logical :: l_so4_aitken = .false.
logical :: l_so4_accu = .false.
logical :: l_so4_diss = .false.
logical :: l_dms = .false.
logical :: l_nh3 = .false.

logical :: l_soot_new = .false.
logical :: l_soot_agd = .false.
logical :: l_soot_cld = .false.

logical :: l_bmass_new = .false.
logical :: l_bmass_agd = .false.
logical :: l_bmass_cld = .false.

logical :: l_ocff_new = .false.
logical :: l_ocff_agd = .false.
logical :: l_ocff_cld = .false.

logical :: l_nitr_acc  = .false.
logical :: l_nitr_diss = .false.

!
! Dependent logicals for Classic aerosol LBC options
! ==================================================
!
logical :: l_so4_aitken_lbc = .false.
logical :: l_so4_accu_lbc = .false.
logical :: l_so4_diss_lbc = .false.

logical :: l_soot_new_lbc = .false.
logical :: l_soot_agd_lbc = .false.
logical :: l_soot_cld_lbc = .false.

logical :: l_bmass_new_lbc = .false.
logical :: l_bmass_agd_lbc = .false.
logical :: l_bmass_cld_lbc = .false.

logical :: l_ocff_new_lbc = .false.
logical :: l_ocff_agd_lbc = .false.
logical :: l_ocff_cld_lbc = .false.

logical :: l_nitr_acc_lbc  = .false.
logical :: l_nitr_diss_lbc = .false.

!
! Heights, i.e. levels, used for manmade aerosol emissions
! ========================================================
!
      ! Aerosol Modelling - model level heights...

integer :: bmass_high_level_1 = imdi ! Lowest and highest
integer :: bmass_high_level_2 = imdi ! for biomass emissions.
integer :: soot_high_level = imdi  ! for chimney soot emission
integer :: ocff_high_level = imdi  ! for chimney OCFF emission

!
! Fixed values no longer in namelist
! ==================================
!
      ! improved nonhydrostatic weights in tracer1 calculation
logical, parameter :: L_tracer1_non_hydro = .true. ! not in namelist

! For Aero_Ctl (Sulph cycle or Soot) No.of times chem called per timestep
integer, parameter :: call_chem_freq = 1 ! Fixed; no longer in namelist

! =======================================================================
!
namelist/RUN_Aerosol/                                                          &
  l_sulpc_so2, l_so2_surfem, l_so2_hilem, l_so2_natem,                         &
  SO2_high_level,                                                              &
  l_sulpc_dms, l_dms_em, l_dms_em_inter, i_dms_flux,                           &
  l_sulpc_ozone, l_sulpc_online_oxidants, l_sulpc_2_way_coupling,              &
  l_sulpc_so2_o3_nonbuffered, l_sulpc_nh3, l_nh3_em,                           &
  l_soot, l_soot_surem, l_soot_hilem, soot_high_level,                         &
  l_biomass, l_bmass_surem, l_bmass_hilem, bmass_high_level_1,                 &
  bmass_high_level_2,                                                          &
  l_ocff, l_ocff_surem, l_ocff_hilem, ocff_high_level, l_nitrate,              &
  l_use_sulphate_sulpc, l_use_nitrate_sulpc, l_use_seasalt_sulpc,              &
  l_use_seasalt_pm, l_temporal_emi, l_use_bmass_sulpc, l_use_ocff_sulpc,       &
  l_so2_lbc, l_dms_lbc, l_nh3_lbc, l_soot_lbc, l_bmass_lbc,                    &
  l_ocff_lbc, l_nitr_lbc, l_bmass_hilem_variable

! DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='RUN_AEROSOL_MOD'

contains
!
! Internal subroutine to check that entries to RUN_aerosol are consistent
! =======================================================================
!
subroutine run_aerosol_check( )

use UM_ParCore, only: mype
use umPrintMgr, only: umPrint, umMessage
use ukca_option_mod, only: i_ukca_chem, l_ukca_classic_hetchem

use ukca_config_specification_mod, only:                                       &
  i_ukca_chem_trop,                                                            &
  i_ukca_chem_raq,                                                             &
  i_ukca_chem_tropisop,                                                        &
  i_ukca_chem_strattrop,                                                       &
  i_ukca_chem_strat

use rad_input_mod, only: l_use_sulpc_direct, l_use_sulpc_indirect_lw,          &
                         l_use_sulpc_indirect_sw, l_use_seasalt_indirect,      &
                         l_use_nitrate_direct, l_use_nitrate_indirect,         &
                         l_use_soot_direct, l_use_soot_indirect,               &
                         l_use_bmass_direct, l_use_bmass_indirect,             &
                         l_use_ocff_direct, l_use_ocff_indirect

use nlsizes_namelist_mod, only:                                                &
    model_levels

implicit none
!
!
character (len=errormessagelength)           :: cmessage
character (len=*), parameter  :: RoutineName = 'RUN_AEROSOL_CHECK'
integer :: icode
logical :: l_ukca_chem
real(kind=jprb) :: zhook_handle
!

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check there is a valid UKCA scheme active
!
l_ukca_chem = (i_ukca_chem==i_ukca_chem_trop)      .or.                        &
              (i_ukca_chem==i_ukca_chem_raq)       .or.                        &
              (i_ukca_chem==i_ukca_chem_tropisop)  .or.                        &
              (i_ukca_chem==i_ukca_chem_strattrop) .or.                        &
              (i_ukca_chem==i_ukca_chem_strat)

! Check that all user-defined schemes and emissions heights are sensible
!
! Sulphur scheme
!
icode = 0

if (l_ukca_classic_hetchem .and. .not. l_sulpc_so2) then
  if (mype == 0) then
    write(umMessage,'(A)')'Error -RUN_aerosol- Cannot run ' //                 &
                          'l_ukca_classic_hetchem without sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_sulpc_dms .and. .not. l_sulpc_so2) then
  if (mype == 0) then
    write(umMessage,'(A)')'Error -RUN_aerosol- sulpc_dms but not sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
else
  l_dms = l_sulpc_dms
end if

if (l_dms_lbc .and. .not. l_dms) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error -RUN_aerosol- dms_lbc but not dms'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_dms_em .and. .not. l_sulpc_dms) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error -RUN_aerosol- dms_em but not sulpc_dms'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_dms_em_inter .and. .not. l_dms_em) then
  if (mype == 0) then
    write(umMessage,'(A)')'Error -RUN_aerosol- dms_em_inter but not dms_em'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (i_dms_flux/=imdi .and. .not. l_dms_em_inter) then
  if (mype == 0) then
    write(umMessage,'(A)')'Error-RUN_aerosol-i_dms_flux but not dms_em_inter'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_sulpc_ozone .and. .not. l_sulpc_so2) then
  if (mype == 0) then
    write(umMessage,'(A)')'Error-RUN_aerosol-sulpc_ozone but not sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_sulpc_nh3 .and. .not. l_sulpc_so2) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error -RUN_aerosol- sulpc_nh3 but not sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
else
  l_nh3 = l_sulpc_nh3
end if

if (l_nh3_lbc .and. .not. l_nh3) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - nh3_lbc but not nh3'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_nh3_em .and. .not. l_sulpc_nh3) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - nh3_em but not sulpc_nh3'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (.not. l_sulpc_ozone .or. .not. l_sulpc_nh3) then
  if (l_sulpc_so2_o3_nonbuffered) then
    if (mype == 0) then
      write(umMessage,'(A)') 'Error -RUN_aerosol- sulpc_so2_o3_nonbuffered'
      call umPrint(umMessage,src='run_aerosol_check')
      write(umMessage,'(A)')'not allowed unless sulpc_ozone/sulpc_nh3 both true'
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  end if
end if

if (l_so2_surfem .and. .not. l_sulpc_so2) then
  if (mype == 0) then
    write(umMessage,'(A)')'Error -RUN_aerosol- so2_surfem but not sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_so2_hilem .and. .not. l_sulpc_so2) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error -RUN_aerosol- so2_hilem but not sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if
if (l_so2_hilem) then
  if (SO2_high_level == imdi) then
    if (mype == 0) then
      write(umMessage,'(A)')'Error-RUN_aerosol-SO2 emission but level not set'
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  else if (SO2_high_level < 1 .or. SO2_high_level > model_levels) then
    if (mype == 0) then
      write(umMessage,'(A,2I6)') 'Error -RUN_aerosol- Query SO2_high_level',   &
                          SO2_high_level,model_levels
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  end if
end if


if (l_so2_natem .and. .not. l_sulpc_so2) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error -RUN_aerosol- so2_natem but not sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_sulpc_so2) then
  ! Dependent logicals normally all needed if l_sulpc_so2 is true
  l_so2 = .true.
  l_so4_aitken = .true.
  l_so4_accu = .true.
  l_so4_diss = .true.
end if

if (l_sulpc_so2 .and. l_so2_lbc) then
  ! Dependent logicals normally all needed if l_so2_lbc also true
  l_so4_aitken_lbc = .true.
  l_so4_accu_lbc = .true.
  l_so4_diss_lbc = .true.
end if

if (l_sulpc_online_oxidants .and.                                              &
            (.not. l_ukca_chem .or. .not. l_sulpc_so2)) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error-RUN_aerosol-sulpc_online_oxidants'
    call umPrint(umMessage,src='run_aerosol_check')
    write(umMessage,'(A)') ' but need both ukca_chem and sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_use_sulphate_sulpc .and. .not. l_sulpc_so2) then
  if (mype == 0) then
    write(umMessage,'(A)')                                                     &
                  'Error-RUN_aerosol-use_sulphate_sulpc but not sulpc_so2'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

!
! SOOT scheme
!
if (l_soot) then
  ! Dependent logicals normally all needed if l_soot is true
  l_soot_new = .true.
  l_soot_agd = .true.
  l_soot_cld = .true.
end if
if (l_soot .and. l_soot_lbc) then
  ! Dependent logicals normally all needed if l_soot_lbc also true
  l_soot_new_lbc = .true.
  l_soot_agd_lbc = .true.
  l_soot_cld_lbc = .true.
end if
!
if (l_soot_surem .and. .not. l_soot) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - soot_surem but not soot'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_soot_hilem .and. .not. l_soot) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - soot_hilem but not soot'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (L_Soot_hilem) then
  if (Soot_high_level == imdi) then
    if (mype == 0) then
      write(umMessage,'(A)')'Error-RUN_aerosol-Soot emission but level not set'
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  else if (Soot_high_level < 1 .or. Soot_high_level > model_levels) then
    if (mype == 0) then
      write(umMessage,'(A,2I6)') 'Error -RUN_aerosol- Query Soot_high_level',  &
                          Soot_high_level,model_levels
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  end if
end if

!
! BIOMASS scheme
!
if (l_biomass) then
  ! Dependent logicals normally all needed if l_biomass is true
  l_bmass_new = .true.
  l_bmass_agd = .true.
  l_bmass_cld = .true.
end if
if (l_biomass .and. l_bmass_lbc) then
  ! Dependent logicals normally all needed if l_bmass_lbc also true
  l_bmass_new_lbc = .true.
  l_bmass_agd_lbc = .true.
  l_bmass_cld_lbc = .true.
end if
!
if (l_bmass_surem .and. .not. l_biomass) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error -RUN_aerosol- bmass_surem but not biomass'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_bmass_hilem .and. .not. l_biomass) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error -RUN_aerosol- bmass_hilem but not biomass'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

! If injection heights for biomass burning emissions are fixed then check that
! the first model layer to inject them has been set correctly by the user.
if (L_bmass_hilem .and. .not. (l_bmass_hilem_variable)) then
  if (bmass_high_level_1 == imdi) then
    if (mype == 0) then
      write(umMessage,'(A)')'Error-RUN_aerosol-bmass emission level 1 not set'
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  else if (bmass_high_level_1<1 .or. bmass_high_level_1>model_levels) then
    if (mype == 0) then
      write(umMessage,'(A,A,2I6)') 'Error - RUN_aerosol - ',                   &
           'Query bmass_high_level_1',bmass_high_level_1,model_levels
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  end if
end if
!
! If injection heights for biomass burning emissions are fixed then check that
! the last model layer to inject them has been set correctly by the user.
if (L_bmass_hilem .and. .not. (l_bmass_hilem_variable)) then
  if (bmass_high_level_2 == imdi) then
    if (mype == 0) then
      write(umMessage,'(A)')'Error-RUN_aerosol-bmass emission level 2 not set'
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  else if (bmass_high_level_2<1 .or. bmass_high_level_2>model_levels) then
    if (mype == 0) then
      write(umMessage,'(A,A,2I6)') 'Error - RUN_aerosol - ',                   &
           'Query bmass_high_level_2',bmass_high_level_2,model_levels
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  end if
end if

if (l_use_bmass_sulpc .and.                                                    &
              (.not. l_biomass .or. .not. l_use_sulphate_sulpc)) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - use_bmass_sulpc '
    call umPrint(umMessage,src='run_aerosol_check')
    write(umMessage,'(A)') 'but not biomass/use_sulphate_sulpc'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

!
! Organic Carbon Fossil Fuel scheme
!
if (l_ocff) then
  ! Dependent logicals normally all needed if l_ocff is true
  l_ocff_new = .true.
  l_ocff_agd = .true.
  l_ocff_cld = .true.
end if
if (l_ocff .and. l_ocff_lbc) then
  ! Dependent logicals normally all needed if l_ocff_lbc also true
  l_ocff_new_lbc = .true.
  l_ocff_agd_lbc = .true.
  l_ocff_cld_lbc = .true.
end if

if (l_ocff_surem .and. .not. l_ocff) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - ocff_surem but not ocff'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_ocff_hilem .and. .not. l_ocff) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - ocff_hilem but not ocff'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (L_ocff_hilem) then
  if (ocff_high_level == imdi) then
    if (mype == 0) then
      write(umMessage,'(A)')'Error-RUN_aerosol-ocff emission but level not set'
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  else if (ocff_high_level < 1 .or. ocff_high_level > model_levels) then
    if (mype == 0) then
      write(umMessage,'(A,2I6)')'Error-RUN_aerosol-Query ocff_high_level',     &
                          ocff_high_level,model_levels
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  end if
end if

if (l_use_ocff_sulpc .and. (.not. l_ocff .or. .not. l_use_sulphate_sulpc)) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - use_ocff_sulpc '
    call umPrint(umMessage,src='run_aerosol_check')
    write(umMessage,'(A)') 'but not ocff/use_sulphate_sulpc'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

!
! Ammonium nitrate aerosol scheme
!
if (l_nitrate .and. (.not. l_nh3 .or. .not. l_ukca_chem)) then
  if (mype == 0) then
    write(umMessage,'(A)')'Error -RUN_aerosol- nitrate but not nh3/ukca_chem'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if

if (l_nitrate) then
  ! Dependent logicals normally all needed if l_nitrate is true
  l_nitr_acc  = .true.
  l_nitr_diss = .true.
end if
if (l_nitrate .and. l_nitr_lbc) then
  ! Dependent logicals normally all needed if l_nitr_lbc also true
  l_nitr_acc_lbc  = .true.
  l_nitr_diss_lbc = .true.
end if

if (.not. l_nitrate .or. .not. l_use_sulphate_sulpc) then
  if (l_use_nitrate_sulpc) then
    if (mype == 0) then
      write(umMessage,'(A)')'Error-RUN_aerosol-use_nitrate_sulpc not allowed'
      call umPrint(umMessage,src='run_aerosol_check')
      write(umMessage,'(A)') 'unless nitrate/use_sulphate_sulpc both true'
      call umPrint(umMessage,src='run_aerosol_check')
    end if
    icode = icode + 1
  end if
end if

!
! Sea Salt scheme
!
if (l_use_seasalt_sulpc .and. .not. l_use_sulphate_sulpc) then
  if (mype == 0) then
    write(umMessage,'(A)') 'Error - RUN_aerosol - use_seasalt_sulpc '
    call umPrint(umMessage,src='run_aerosol_check')
    write(umMessage,'(A)') 'but not use_sulphate_sulpc'
    call umPrint(umMessage,src='run_aerosol_check')
  end if
  icode = icode + 1
end if


if (icode /= 0) then
  cmessage='*** Error(s) in levels or logicals set in RUN_Aerosol namelist'
  call ereport(RoutineName, icode, cmessage)
end if

! Also check validity of aerosol radiative effect choices here,
! as this needs to be done after run_aerosol has been read.
if (.not. l_sulpc_so2) then
  l_use_sulpc_direct = .false.
  l_use_sulpc_indirect_lw = .false.
  l_use_sulpc_indirect_sw = .false.
end if
if (.not. l_use_sulpc_indirect_lw .or. .not. l_use_sulpc_indirect_sw) then
  l_use_seasalt_indirect = .false.
end if
if (.not. l_nitrate) then
  l_use_nitrate_direct = .false.
  l_use_nitrate_indirect = .false.
else if (.not. l_use_sulpc_indirect_lw .or. .not. l_use_sulpc_indirect_sw) then
  l_use_nitrate_indirect = .false.
end if
if (.not. l_soot) then
  l_use_soot_direct = .false.
end if
l_use_soot_indirect = .false.  ! Always false
if (.not. l_biomass) then
  l_use_bmass_direct = .false.
  l_use_bmass_indirect = .false.
else if (.not. l_use_sulpc_indirect_lw .or. .not. l_use_sulpc_indirect_sw) then
  l_use_bmass_indirect = .false.
end if
if (.not. l_ocff) then
  l_use_ocff_direct = .false.
  l_use_ocff_indirect = .false.
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return

end subroutine run_aerosol_check

subroutine print_nlist_run_aerosol()
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_AEROSOL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_aerosol',                               &
    src='run_aerosol_mod')

write(lineBuffer,'(A,L1)')' l_sulpc_so2 = ',l_sulpc_so2
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_so2_surfem = ',l_so2_surfem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_so2_hilem = ',l_so2_hilem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_so2_natem = ',l_so2_natem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_so2_lbc = ',l_so2_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_sulpc_dms = ',l_sulpc_dms
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_dms_em = ',l_dms_em
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_dms_em_inter = ',l_dms_em_inter
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,I0)')' i_dms_flux = ',i_dms_flux
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_sulpc_ozone = ',l_sulpc_ozone
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_sulpc_online_oxidants = ',l_sulpc_online_oxidants
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_sulpc_nh3 = ',l_sulpc_nh3
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_nh3_em = ',l_nh3_em
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)') 'l_sulpc_so2_o3_nonbuffered=',                      &
                                                    l_sulpc_so2_o3_nonbuffered
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_soot = ',l_soot
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_soot_surem = ',l_soot_surem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_soot_hilem = ',l_soot_hilem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_biomass = ',l_biomass
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_bmass_surem = ',l_bmass_surem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_bmass_hilem = ',l_bmass_hilem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_ocff = ',l_ocff
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_ocff_surem = ',l_ocff_surem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_ocff_hilem = ',l_ocff_hilem
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_nitrate = ',l_nitrate
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_use_sulphate_sulpc = ',l_use_sulphate_sulpc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_use_nitrate_sulpc = ',l_use_nitrate_sulpc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_use_seasalt_sulpc = ',l_use_seasalt_sulpc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_use_seasalt_pm = ',l_use_seasalt_pm
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A18, L7)')' l_temporal_emi = ',l_temporal_emi
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_use_bmass_sulpc = ',l_use_bmass_sulpc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_use_ocff_sulpc = ',l_use_ocff_sulpc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,I0)')' so2_high_level = ',so2_high_level
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,I0)')' soot_high_level = ',soot_high_level
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,I0)')' bMass_high_level_1 = ',bMass_high_level_1
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,I0)')' bMass_high_level_2 = ',bMass_high_level_2
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,I0)')' ocff_high_level = ',ocff_high_level
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_so2_lbc = ',l_so2_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_dms_lbc = ',l_dms_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_so4_aitken_lbc = ',l_so4_aitken_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_so4_accu_lbc = ',l_so4_accu_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_so4_diss_lbc = ',l_so4_diss_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_nh3_lbc = ',l_nh3_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_soot_new_lbc = ',l_soot_new_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_soot_agd_lbc = ',l_soot_agd_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_soot_cld_lbc = ',l_soot_cld_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_bmass_new_lbc = ',l_bmass_new_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_bmass_agd_lbc = ',l_bmass_agd_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_bmass_cld_lbc = ',l_bmass_cld_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_ocff_new_lbc = ',l_ocff_new_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_ocff_agd_lbc = ',l_ocff_agd_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_ocff_cld_lbc = ',l_ocff_cld_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_nitr_acc_lbc = ',l_nitr_acc_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A,L1)')' l_nitr_diss_lbc = ',l_nitr_diss_lbc
call umPrint(lineBuffer,src='run_aerosol_mod')
write(lineBuffer,'(A26, L1)') ' l_bmass_hilem_variable = ',                    &
                                l_bmass_hilem_variable
call umPrint(lineBuffer,src='run_aerosol_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='run_aerosol_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_aerosol


end module run_aerosol_mod
