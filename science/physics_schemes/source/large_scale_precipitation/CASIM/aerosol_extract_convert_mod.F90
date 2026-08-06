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

! This module contains the declarations and subroutines for the conversion
! and incorporation of aerosol from UKCA into the CASIM microphysics.
!
!
module aerosol_extract_convert_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='AEROSOL_EXTRACT_CONVERT_MOD'

contains

subroutine aerosol_extract_convert(  p_layer_centres, t_work, rhoCasim,        &
                                     tracer_ukca,                              &
                                     AitkenSolMass, AitkenSolNumber,           &
                                     AccumSolMass,  AccumSolNumber,            &
                                     CoarseSolMass, CoarseSolNumber,           &
                                     AccumDustMass, AccumDustNumber,           &
                                     CoarseDustMass, CoarseDustNumber,         &
                                     AitkenSolBk, AccumSolBk, CoarseSolBk )

! Subroutine to read aerosol fields from UKCA
! This routine returns the soluble aitken, accumulation, coarse,
! and insoluble accumulation, coarse mass and number.

use atm_fields_bounds_mod, only: tdims, tdims_s
use nlsizes_namelist_mod,  only: tr_ukca

use ukca_config_specification_mod, only: glomap_variables

use ukca_mode_setup,       only: mode_ait_sol, mode_acc_sol,                   &
                                 mode_cor_sol, mode_acc_insol,                 &
                                 mode_cor_insol

use ukca_scavenging_mod, only: nmr_index_um,mmr_index_um
use ukca_constants,      only: mmw
use water_constants_mod, only: rho_water
use chemistry_constants_mod, only: boltzmann
use missing_data_mod, only: rmdi

! Dr Hook Modules
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

implicit none

!-------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Pressure at layer centres [Pa]
real(kind=real_umphys), intent(in) ::                                          &
                     p_layer_centres( tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end,             &
                                                  0 : tdims%k_end )

! Local Working temperature [K]
real(kind=real_umphys), intent(in) ::                                          &
                    t_work(           tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end,             &
                                                  1 : tdims%k_end )

!Air density for CASIM [kg m-3]
real(kind=real_umphys), intent(in) ::                                          &
                    rhoCasim(                     1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! UKCA tracers []
real(kind=real_umphys), intent(in) ::                                          &
                    tracer_ukca(    tdims_s%i_start : tdims_s%i_end,           &
                                    tdims_s%j_start : tdims_s%j_end,           &
                                    tdims_s%k_start : tdims_s%k_end, tr_ukca )

! Accumulation Model Soluble Mass [kg kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     AccumSolMass(                1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Accumulation Model Soluble Number [kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     AccumSolNumber(              1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Aitken mode Soluble Mass [kg kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     AitkenSolMass(               1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Accumulation Model Soluble Number [kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     AitkenSolNumber(             1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Coarse Model Soluble Mass [kg kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     CoarseSolMass(               1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Coarse Model Soluble Number [kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     CoarseSolNumber(             1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Coarse Dust Mass [kg kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     CoarseDustMass(              1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Coarse Dust Number [kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     CoarseDustNumber(            1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Accumulation Mode Dust Mass [kg kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     AccumDustMass(               1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Accumulation Mode Dust Number [kg-1]
real(kind=real_umphys), intent(out) ::                                         &
                     AccumDustNumber(             1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end )

! Abdul-Razzak-Ghan parameters (volume weighted)
real(kind=real_umphys), intent(out) ::                                         &
                             AccumSolBk(          1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end)
real(kind=real_umphys), intent(out) ::                                         &
                             AitkenSolBk(         1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end)
real(kind=real_umphys), intent(out) ::                                         &
                             CoarseSolBk(         1 : tdims%k_end,             &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end)
!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
logical, pointer :: component(:,:)
real,    pointer :: mm (:)
logical, pointer :: mode (:)
integer, pointer :: ncp
real,    pointer :: no_ions (:)
real,    pointer :: rhocomp (:)

! Local working density [kg m-3]

real(kind=real_umphys) ::                                                      &
                AccumSolVolume(       1 : tdims%k_end,                         &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end)
real(kind=real_umphys) ::                                                      &
                AitkenSolVolume(      1 : tdims%k_end,                         &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end)
real(kind=real_umphys) ::                                                      &
                CoarseSolVolume(      1 : tdims%k_end,                         &
                                      tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end)
real(kind=real_umphys) :: aird
integer :: i, j, k ! Loop counters

integer :: imode,icp,i_cpt ! UKCA counters

real(kind=real_umphys), parameter :: real_eps = epsilon(1.0_real_umphys)

character(len=*), parameter :: RoutineName='AEROSOL_EXTRACT_CONVERT'

! Declarations for Dr Hook
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------
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

! sum up component masses in each mode for total mass
! and copy number to inputs for CASIM.

! initialise modal Mass arrays to zero (will contain sum of mass over all
!  components)

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED( tdims, AitkenSolMass, AitkenSolNumber, CoarseSolMass,            &
!$OMP         CoarseSolNumber, AccumDustMass, AccumDustNumber,                 &
!$OMP         CoarseDustMass,  CoarseDustNumber,  AccumSolMass,                &
!$OMP         AccumSolNumber, AitkenSolBk, AccumSolBk, CoarseSolBk,            &
!$OMP         AitkenSolVolume, AccumSolVolume, CoarseSolVolume,                &
!$OMP         tracer_ukca, p_layer_centres, no_ions, mm, rhocomp,              &
!$OMP         rhoCASIM, T_work, mode, component, ncp, mmr_index_um,            &
!$OMP         nmr_index_um)                                                    &
!$OMP private( i, j, k, aird, imode, icp, i_cpt )
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    do k = 1, tdims%k_end
      AitkenSolMass (k,i,j)=0.0
      AccumSolMass  (k,i,j)=0.0
      CoarseSolMass (k,i,j)=0.0
      AccumDustMass (k,i,j)=0.0
      CoarseDustMass(k,i,j)=0.0
      AitkenSolBk   (k,i,j)=0.0
      AccumSolBk    (k,i,j)=0.0
      CoarseSolBk   (k,i,j)=0.0
      AitkenSolVolume(k,i,j)=0.0
      AccumSolVolume(k,i,j)=0.0
      CoarseSolVolume(k,i,j)=0.0
    end do ! loop over k
  end do ! loop over i
end do  ! loop over j
!$OMP end do
!
! include modes 2 to 6 (Aitsol,accsol,corsol,Aitins,accins,corins)
do imode=mode_ait_sol,mode_cor_insol
  if (mode(imode)) then
    do icp=1,ncp
      i_cpt=0
      if (component(imode,icp)) then
        i_cpt=mmr_index_um(imode,icp)
        !The code in the CASIM activation scheme that is replaced by this interface:
        !start:
        ! Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rho_water)
        ! no_ions,rho_comp from ukca_mode_setup,massMole is 0.132
        !which is ammonium sulphate (18*2+96). Density is set to 1777 in
        ! mphys_constants.F90.
        !end
        ! The Bk value here does not assume ammonium sulphate but is a volume
        ! weighted average of the Bk values that would come from the UKCA components

        ! For each active mode calculate mass to transfer to CASIM, and Bk numerator
        ! and denominator
        select case (imode)

        case (mode_ait_sol)
!$OMP do SCHEDULE(STATIC)
          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              do k = 1, tdims%k_end
                AitkenSolMass(k,i,j)   = AitkenSolMass(k,i,j) +                &
                                         tracer_ukca(i,j,k,i_cpt)
                AitkenSolBk(k,i,j)     = AitkenSolBk(k,i,j) +no_ions(icp) *    &
                                         mmw * tracer_ukca(i,j,k,i_cpt) /      &
                                         (rho_water*mm(icp))
                AitkenSolVolume(k,i,j) = AitkenSolVolume(k,i,j) +              &
                                         tracer_ukca(i,j,k,i_cpt) / rhocomp(icp)
              end do ! loop over k
            end do ! loop over i
          end do  ! loop over j
!$OMP end do
        case (mode_acc_sol)
!$OMP do SCHEDULE(STATIC)
          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              do k = 1, tdims%k_end
                AccumSolMass(k,i,j)    = AccumSolMass(k,i,j) +                 &
                                         tracer_ukca(i,j,k,i_cpt)
                AccumSolBk(k,i,j)      = AccumSolBk(k,i,j) + no_ions(icp) *    &
                                         mmw * tracer_ukca(i,j,k,i_cpt) /      &
                                         (rho_water*mm(icp))
                AccumSolVolume(k,i,j)  = AccumSolVolume(k,i,j) +               &
                                         tracer_ukca(i,j,k,i_cpt)/rhocomp(icp)
              end do ! loop over k
            end do ! loop over i
          end do  ! loop over j
!$OMP end do
        case (mode_cor_sol)
!$OMP do SCHEDULE(STATIC)
          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              do k = 1, tdims%k_end
                CoarseSolMass(k,i,j)   = CoarseSolMass(k,i,j) +                &
                                         tracer_ukca(i,j,k,i_cpt)
                CoarseSolBk(k,i,j)     = CoarseSolBk(k,i,j) + no_ions(icp) *   &
                                         mmw * tracer_ukca(i,j,k,i_cpt) /      &
                                         (rho_water*mm(icp))
                CoarseSolVolume(k,i,j) = CoarseSolVolume(k,i,j) +              &
                                         tracer_ukca(i,j,k,i_cpt) / rhocomp(icp)
              end do ! loop over k
            end do ! loop over i
          end do  ! loop over j
!$OMP end do
        case (mode_acc_insol)
!$OMP do SCHEDULE(STATIC)
          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              do k = 1, tdims%k_end
                AccumDustMass (k,i,j) = AccumDustMass(k,i,j) +                 &
                                        tracer_ukca(i,j,k,i_cpt)
              end do ! loop over k
            end do ! loop over i
          end do  ! loop over j
!$OMP end do
        case (mode_cor_insol)
!$OMP do SCHEDULE(STATIC)
          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              do k = 1, tdims%k_end
                CoarseDustMass(k,i,j) = CoarseDustMass(k,i,j) +                &
                                        tracer_ukca(i,j,k,i_cpt)
              end do ! loop over k
            end do ! loop over i
          end do  ! loop over j
!$OMP end do
        end select
      end if ! if(component(imode,icp)
    end do ! loop over icp

    ! For each active mode calculate number to transfer to CASIM
    i_cpt = nmr_index_um(imode)
    select case (imode)

    case (mode_ait_sol)
!$OMP do SCHEDULE(STATIC)
      do k=1,tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            aird = p_layer_centres(i,j,k)/                                     &
                  ( T_work(i,j,k)*boltzmann )
            AitkenSolNumber(k,i,j)          =                                  &
                  tracer_ukca(i,j,k,i_cpt)*aird/rhoCASIM(k,i,j)
          end do ! loop over i
        end do ! loop over j
      end do  ! loop over k
!$OMP end do
    case (mode_acc_sol)
!$OMP do SCHEDULE(STATIC)
      do k=1,tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            aird = p_layer_centres(i,j,k)/                                     &
                  ( T_work(i,j,k)*boltzmann )
            AccumSolNumber(k,i,j)           =                                  &
                  tracer_ukca(i,j,k,i_cpt)*aird/rhoCASIM(k,i,j)
          end do ! loop over i
        end do ! loop over j
      end do  ! loop over k
!$OMP end do
    case (mode_cor_sol)
!$OMP do SCHEDULE(STATIC)
      do k=1,tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            aird = p_layer_centres(i,j,k)/                                     &
                  ( T_work(i,j,k)*boltzmann )
            CoarseSolNumber(k,i,j)          =                                  &
                  tracer_ukca(i,j,k,i_cpt)*aird/rhoCASIM(k,i,j)
          end do ! loop over i
        end do ! loop over j
      end do  ! loop over k
!$OMP end do
    case (mode_acc_insol)
!$OMP do SCHEDULE(STATIC)
      do k=1,tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            aird = p_layer_centres(i,j,k)/                                     &
                  ( T_work(i,j,k)*boltzmann )
            AccumDustNumber(k,i,j)          =                                  &
                  tracer_ukca(i,j,k,i_cpt)*aird/rhoCASIM(k,i,j)
          end do ! loop over i
        end do ! loop over j
      end do  ! loop over k
!$OMP end do
    case (mode_cor_insol)
!$OMP do SCHEDULE(STATIC)
      do k=1,tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            aird = p_layer_centres(i,j,k)/                                     &
                  ( T_work(i,j,k)*boltzmann )
            CoarseDustNumber(k,i,j)         =                                  &
                  tracer_ukca(i,j,k,i_cpt)*aird/rhoCASIM(k,i,j)
          end do ! loop over i
        end do  ! loop over j
      end do  ! loop over k
!$OMP end do
    end select
  end if ! if(mode(imode))
end do ! loop over imode

! calculate activation parameters Bk for soluble Aitken, accumulation,
! coarse modes. These modes are always used in UKCA, even if not always
! used in CASIM
! If there is no aerosol in a mode then CASIM will not carry out activation
! in that mode. In order to avoid carrying out spurious divides by zero here,
! we set the Bk to missing data if the volume is zero.
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    do k = 1, tdims%k_end
      if (abs(AitkenSolVolume(k,i,j)) < real_eps) then
        AitkenSolBk(k,i,j) = rmdi
      else
        AitkenSolBk(k,i,j) = AitkenSolBk(k,i,j)/AitkenSolVolume(k,i,j)
      end if
    end do
  end do
end do
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    do k = 1, tdims%k_end
      if (abs(AccumSolBk(k,i,j)) < real_eps) then
        AccumSolBk(k,i,j) = rmdi
      else
        AccumSolBk(k,i,j)  = AccumSolBk(k,i,j)/AccumSolVolume(k,i,j)
      end if
    end do
  end do
end do
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    do k = 1, tdims%k_end
      if (abs(CoarseSolBk(k,i,j)) < real_eps) then
        CoarseSolBk(k,i,j) = rmdi
      else
        CoarseSolBk(k,i,j) = CoarseSolBk(k,i,j)/CoarseSolVolume(k,i,j)
      end if
    end do
  end do
end do
!$OMP end do NOWAIT
!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine aerosol_extract_convert

end module aerosol_extract_convert_mod
