! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!-------------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!-------------------------------------------------------------------------------

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

module casim_activation_in_um_mod

use micro_main, only: aerophys, aeroact, aerochem, dustact, dustphys,          &
                      dustchem, aeroice, dustliq



implicit none

private

public :: examine_ukca_aerosol_column, activate_column_ukca

character(len=*), parameter, private ::                                        &
    ModuleName = 'CASIM_ACTIVATION_IN_UM_MOD'

contains

subroutine examine_ukca_aerosol_column(tracer_ukca, t_col, p_col, rho_col)

!==========================================================================
! Description:
! This subroutine takes the ukca tracers (code variable tracer_ukca) and
! Sets up the CASIM aerosol species ready for an activation call outside
! of the main CASIM call.
!==========================================================================

use um_types,              only: real_umphys
use atm_fields_bounds_mod, only: tdims
use nlsizes_namelist_mod,  only: tr_ukca
use casim_switches,        only: l_casim_warm_only

use mphys_switches,        only: aero_index
use lognormal_funcs,       only: MNtoRm
use mphys_parameters,      only: nz
use thresholds,            only: ccn_tidy, aeromass_small, aeronumber_small

use ukca_config_specification_mod,  only: glomap_variables

use ukca_mode_setup,       only: mode_ait_sol, mode_acc_sol,                   &
                                 mode_cor_sol, mode_acc_insol,                 &
                                 mode_cor_insol

use ukca_scavenging_mod,     only: nmr_index_um,mmr_index_um
use ukca_constants,          only: mmw
use water_constants_mod,     only: rho_water
use chemistry_constants_mod, only: boltzmann
use yomhook,                 only: lhook, dr_hook
use parkind1,                only: jprb, jpim

implicit none

! Subroutine arguments (all intent(in) for now).
real(kind=real_umphys), intent(in) :: tracer_ukca(0:tdims%k_end, tr_ukca)
real(kind=real_umphys), intent(in) :: t_col(1:tdims%k_end)
real(kind=real_umphys), intent(in) :: p_col(1:tdims%k_end)
real(kind=real_umphys), intent(in) :: rho_col(1:tdims%k_end)

! Local variables

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
logical, pointer :: component(:,:)
real,    pointer :: mm (:)
logical, pointer :: mode (:)
integer, pointer :: ncp
real,    pointer :: no_ions (:)
real,    pointer :: rhocomp (:)

! Aerosol density and air density
real(kind=real_umphys) :: density, aird(nz)

integer :: k ! loop counter

integer :: top_level ! top model level to work with

integer :: imode, jmode, i_cpt, icp ! Aerosol modes and counters

! Aerosol real variables (number, mass, numerator and volume).
real(kind=real_umphys) :: mode_N, mode_M, bk_numerator, volume, mode_Bk

character(len=*), parameter :: RoutineName='EXAMINE_UKCA_AEROSOL_COLUMN'

real(kind=real_umphys), parameter :: zero       = 0.0_real_umphys
real(kind=real_umphys), parameter :: rcrit_high = 999.0_real_umphys

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

real(kind=real_umphys) :: eps_1 ! epsilon of 1.0

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mm          => glomap_variables%mm
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
no_ions     => glomap_variables%no_ions
rhocomp     => glomap_variables%rhocomp

eps_1 = epsilon(1.0_real_umphys)
top_level = min(tdims%k_end, nz)

!--------------------------------------------------
! Set up Soluble Aerosol Species
!--------------------------------------------------

do k = 1, top_level
  density = aerochem(k) % density(1)
  aird(k) = p_col(k) / (t_col(k)*boltzmann)

  do imode = 1, aero_index%nccn

    mode_N       = zero
    mode_M       = zero
    bk_numerator = zero
    volume       = zero
    mode_Bk      = zero

    if (imode  ==  aero_index%i_aitken) then
      jmode = mode_ait_sol
    else if (imode  ==  aero_index%i_accum) then
      jmode = mode_acc_sol
    else if (imode  ==  aero_index%i_coarse) then
      jmode = mode_cor_sol
    end if

    ! Loop over aerosol modes to set up the aerosol physics
    if (mode(jmode)) then
      i_cpt  = nmr_index_um(jmode)
      mode_N = tracer_ukca(k, i_cpt) * aird(k) / rho_col(k)

      do icp = 1, ncp
        if ( component(jmode, icp) ) then

          i_cpt  = mmr_index_um(jmode,icp)
          mode_M = mode_M + tracer_ukca(k, i_cpt)

          bk_numerator = bk_numerator + no_ions(icp) * mmw *                   &
                         tracer_ukca(k, i_cpt) / (rho_water*mm(icp))

          volume = volume + tracer_ukca(k, i_cpt) / rhocomp(icp)
        end if ! component(jmode, icp)
      end do ! icp

      if (volume > zero) then
        mode_Bk = bk_numerator / volume
      else
        mode_Bk = zero
      end if

    end if ! mode(jmode)

    if (mode_N > ccn_tidy .and. mode_M > ccn_tidy * eps_1 .and. &
         mode_Bk > zero) then
      aerophys(k)%n(imode)  = mode_N
      aerophys(k)%m(imode)  = mode_M
      aerophys(k)%rd(imode) = MNtoRm( mode_M, mode_N, density,                 &
                                      aerophys(k)%sigma(imode) )
      aerochem(k)%bk(imode) = mode_Bk
    else
      ! Zero the processes.
      aerophys(k)%n(imode)  = zero
      aerophys(k)%m(imode)  = zero
      aerophys(k)%rd(imode) = zero
      aerochem(k)%bk(imode) = zero

    end if ! mode_N > ccn_tidy etc


  end do ! imode

  ! Set other aerosol properties for this level
  aeroact(k) % nact       = zero
  aeroact(k) % mact       = zero
  aeroact(k) % rcrit      = rcrit_high
  aeroact(k) % mact_mean  = zero
  aeroact(k) % nact2      = zero
  aeroact(k) % nact1      = zero
  aeroact(k) % mact1      = zero
  aeroact(k) % rcrit1     = rcrit_high
  aeroact(k) % nact2      = zero
  aeroact(k) % rcrit2     = rcrit_high
  aeroact(k) % mact2      = zero
  aeroact(k) % mact2_mean = zero
  aeroact(k) % mact1_mean = zero
  aeroact(k) % rcrit1     = rcrit_high
  aeroact(k) % nact2      = zero
  aeroact(k) % rcrit2     = rcrit_high
  aeroact(k) % mact2      = zero
  aeroact(k) % mact2_mean = zero
  aeroact(k) % mact1_mean = zero
  aeroact(k) % nratio1    = zero
  aeroact(k) % nratio2    = zero
end do ! loop over levels (k)

!--------------------------------------------------
! Set up Insoluble Aerosol Species
!--------------------------------------------------

if (.not. l_casim_warm_only ) then
  ! Activated dust
  do k = 1, nz
    density = dustchem(k) % density(1)
    ! Examine interstitial dust
    do imode = 1, aero_index % nin

      mode_N = zero
      mode_M = zero

      if (imode  ==  aero_index % i_coarse_dust) then
        jmode = mode_cor_insol
      end if
      if (imode  ==  aero_index % i_accum_dust) then
        jmode = mode_acc_insol
      end if
      if (mode(jmode)) then
        i_cpt  = nmr_index_um(jmode)
        mode_N = tracer_ukca(k, i_cpt) * aird(k) / rho_col(k)

        do icp = 1, ncp
          if ( component(jmode, icp) ) then
            i_cpt  = mmr_index_um(jmode, icp)
            mode_M = mode_M + tracer_ukca(k, i_cpt)
          end if
        end do

      end if
      if (mode_m < zero .or. mode_n < zero ) then
        if (mode_m < zero) mode_m = aeromass_small
        if (mode_n < zero) mode_n = aeronumber_small
      end if

      dustphys(k)%n(imode) = mode_N
      dustphys(k)%m(imode) = mode_M

      dustphys(k)%rd(imode) = MNtoRm(mode_M, mode_N, density,                  &
                                     dustphys(k)%sigma(imode) )
    end do

    ! Examine activated dust
    ! Initialize to zero/defaults
    dustact(k) % nact       = zero
    dustact(k) % mact       = zero
    dustact(k) % rcrit      = rcrit_high
    dustact(k) % mact_mean  = zero
    dustact(k) % nact1      = zero
    dustact(k) % nratio1    = zero
    dustact(k) % mact1      = zero
    dustact(k) % rcrit1     = rcrit_high
    dustact(k) % mact1_mean = zero
    dustact(k) % mact2      = zero
    dustact(k) % nact2      = zero
    dustact(k) % nratio2    = zero
    dustact(k) % rcrit2     = rcrit_high
    dustact(k) % mact2_mean = zero
    dustact(k) % mact3      = zero
    dustact(k) % nact3      = zero
    dustact(k) % nratio3    = zero
    dustact(k) % rcrit3     = rcrit_high
    dustact(k) % mact3_mean = zero

  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine examine_ukca_aerosol_column

!==============================================================================

subroutine activate_column_ukca( cloud_mass_post_qtbal,                        &
                                 cloud_mass_pre_ap2,                           &
                                 rho_col, t_col, p_col,                        &
                                 cf_liquid, cf_liquid_pre_ap2, w_tke,          &
                                 cloudnumber_col,                              &
                                 tracer_ukca_col )

!----------------------------------------------------------------------------
! Code Description:
! Changes to the liquid and ice cloud mass can occur due to processes outside
! of the CASIM microphysics. These are often in the slow physics terms
! such as convection and boundary layer. For double-moment CASIM simulations,
! the cloud number and ice number have to be set in order to ensure that
! CASIM does not reject this new liquid and ice cloud. For UKCA simulations,
! this also means checking and activating the aerosol within the model column
! before calling the 'activate' routine from the CASIM repository to
! calculate the cloud number in the same way it would normally be done in
! the main CASIM call.
!----------------------------------------------------------------------------

use atm_fields_bounds_mod, only: tdims
use timestep_mod,          only: timestep, recip_timestep
use um_types,              only: real_umphys
use nlsizes_namelist_mod,  only: tr_ukca
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use mphys_parameters,      only: nz
use mphys_switches,        only: aero_index
use thresholds,            only: ql_tidy, cfliq_small
use condensation,          only: dnccn_all, dmac_all, dnccnd_all, dmad_all
use activation,            only: activate

implicit none

! Subroutine arguments
real(kind=real_umphys), intent(in) :: cloud_mass_post_qtbal(1:tdims%k_end)

real(kind=real_umphys), intent(in) :: cloud_mass_pre_ap2(1:tdims%k_end)

real(kind=real_umphys), intent(in) :: rho_col(1:tdims%k_end)

real(kind=real_umphys), intent(in) :: t_col(1:tdims%k_end)

real(kind=real_umphys), intent(in) :: p_col(1:tdims%k_end)

real(kind=real_umphys), intent(in) :: cf_liquid(1:tdims%k_end)

real(kind=real_umphys), intent(in) :: cf_liquid_pre_ap2(1:tdims%k_end)

real(kind=real_umphys), intent(in) :: w_tke(1:tdims%k_end)

real(kind=real_umphys), intent(in) :: tracer_ukca_col(0:tdims%k_end, tr_ukca)

real(kind=real_umphys), intent(in out) :: cloudnumber_col(1:tdims%k_end)

! Local variables

real(kind=real_umphys) :: w_cloud_number

real(kind=real_umphys) :: delta_mass
real(kind=real_umphys) :: delta_cfliq
real(kind=real_umphys) :: dnumber, dnumber_d, dnumber_a
real(kind=real_umphys) :: dmac
real(kind=real_umphys) :: dmad

real(kind=real_umphys) :: smax,ait_cdnc,accum_cdnc,tot_cdnc,                   &
                          activated_cloud, activated_arg



real(kind=real_umphys), parameter :: zero = 0.0_real_umphys

integer :: k ! loop counter
integer :: top_level ! highest level to use.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle
character(len=*), parameter   :: RoutineName='ACTIVATE_COLUMN_UKCA'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------
! Initial set up of constants
!-----------------------------------
! Check for bounds - the top level cannot be greater than nz or
! tdims%k_end.
top_level = min(nz, tdims%k_end)

!------------------------------------------------------------------------------
! Call examine_ukca_aerosol_column. This uses the UKCA tracers to set up the
! aerosol amounts within CASIM for a single column. The internal aerosol
! amounts within CASIM are required for the activate routine to be able to
! set the cloud number based on the aerosol in each mode.
!------------------------------------------------------------------------------
call examine_ukca_aerosol_column( tracer_ukca_col, t_col, p_col, rho_col )

!------------------------------------------------------------------------
! Main loop over levels to call the activate routine from within CASIM.
!------------------------------------------------------------------------
do k = 1, top_level

   ! Work out if cloud mass or liquid cloud fraction has changed.
  delta_mass  = cloud_mass_post_qtbal(k) - cloud_mass_pre_ap2(k)
  delta_cfliq = cf_liquid(k) - cf_liquid_pre_ap2(k)

  ! Set local working cloud number for this level.
  w_cloud_number = cloudnumber_col(k)

  if ( delta_mass > zero .and. cloud_mass_post_qtbal(k) > ql_tidy .and.        &
       cf_liquid(k) > cfliq_small ) then

    !-----------------------------------------------------------------
    ! There has been an increase in liquid cloud due to non-CASIM
    ! physical processes in the model.
    !-----------------------------------------------------------------
    ! Call CASIM routine activate (operates over one grid box) in
    ! order to work out what the new cloud number is based on the
    ! aerosol and the cloud mass added by the other processes.
    !-----------------------------------------------------------------

    call activate(timestep, cloud_mass_post_qtbal(k), w_cloud_number,          &
                  w_tke(k), rho_col(k), dnumber, dmac, t_col(k), p_col(k),     &
                  cf_liquid(k), cf_liquid_pre_ap2(k),                          &
                  aerophys(k), aerochem(k), aeroact(k),                        &
                  dustphys(k), dustchem(k), dustliq(k),                        &
                  dnccn_all, dmac_all, dnumber_d, dmad,                        &
                  dnccnd_all, dmad_all,                                        &
                  smax,ait_cdnc,accum_cdnc, tot_cdnc,activated_arg,            &
                  activated_cloud)

    dnumber_a = dnumber

  else if ( cloud_mass_post_qtbal(k) < ql_tidy       .or.                      &
            cf_liquid(k)             < cfliq_small ) then

    !-----------------------------------------------------------------
    ! Non-CASIM processes have removed all of the liquid cloud or the
    ! liquid cloud fraction. Therefore, we need to remove the cloud
    ! number
    !-----------------------------------------------------------------

    dnumber = -w_cloud_number * recip_timestep

    w_cloud_number = zero

  else

    if (delta_cfliq < zero) then
      dnumber   = w_cloud_number * delta_cfliq / cf_liquid_pre_ap2(k)  *       &
                  recip_timestep
      dmac      = zero
      dnumber_a = dnumber
      dmac_all( aero_index % i_accum)  = dmac
      ! put back into accum
      dnccn_all( aero_index % i_accum) = dnumber_a
      dnumber_d  = zero
      dnccnd_all = zero
      dmad_all   = zero
    else
      dnumber    = zero
      dnumber_a  = zero ! No aerosol processing required

      dnccn_all  = zero ! we assume no change in number during evap
      dmac       = zero ! No aerosol processing required
      dmac_all   = zero ! No aerosol processing required
      dnumber_d  = zero ! No aerosol processing required
      dnccnd_all = zero
      dmad_all   = zero
    end if ! delta_cfliq < zero

  end if ! delta_mass > zero etc.

  ! Update cloud number in the column
  w_cloud_number     = max(w_cloud_number + (dnumber * timestep), zero)
  cloudnumber_col(k) = w_cloud_number

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine activate_column_ukca

end module casim_activation_in_um_mod
