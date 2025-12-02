! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

!************************************
! Large-scale precipitation scheme.
! Heterogeneous nucleation temperature (tnuc)
! is defined as a function of dust number concentration.
! Developed at National Institute of Water and Atmospheric Research Ltd. (NIWA)
! Collaboration with UK MetOffice to improve the short-wave radiation
! biases in the model

!************************************
! subroutine lsp_prognostic_tnuc
!************************************

module lsp_prognostic_tnuc_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='LSP_PROGNOSTIC_TNUC_MOD'

contains

subroutine lsp_prognostic_tnuc(                                                &
                    dust_div1                                                  &
                   ,dust_div2                                                  &
                   ,dust_div3                                                  &
                   ,dust_div4                                                  &
                   ,dust_div5                                                  &
                   ,dust_div6                                                  &
                   ,tracer_ukca_iacc                                           &
                   ,tracer_ukca_icor                                           &
                   ,theta                                                      &
                   ,exner_theta_levels                                         &
                   ,p_theta_levels                                             &
                   ,dust_tot_nd, tnuc_new                                      &
                   ,arcldust_b1                                                &
                   ,arcldust_b2                                                &
                   ,arcldust_b3                                                &
                   ,arcldust_b4                                                &
                   ,arcldust_b5                                                &
                   ,arcldust_b6                                                &
                   )


use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
use conversions_mod, only: pi
use lsprec_mod, only: zero
use dust_parameters_mod, only: drep,rhop,l_twobin_dust,i_dust
use planet_constants_mod, only: r
use atm_fields_bounds_mod, only: tdims_s
use nlsizes_namelist_mod,  only: row_length, rows, model_levels
use mphys_ice_mod, only: thomo
use mphys_constants_mod, only:tnuc
use rad_input_mod,       only: l_use_arcldust
use errormessagelength_mod, only: errormessagelength
use ereport_mod, only: ereport
use umprintmgr,  only: newline
use chemistry_constants_mod, only: boltzmann
use ukca_option_mod, only: l_ukca_dust

implicit none

integer :: i,j,k,l
integer :: ErrorStatus
character(len=errormessagelength) :: comments

real(kind=real_umphys),intent(in) ::                                           &
    dust_div1(tdims_s%i_start:tdims_s%i_end,                                   &
              tdims_s%j_start:tdims_s%j_end,                                   &
              tdims_s%k_start:tdims_s%k_end)
               !dust mmr in div1
real(kind=real_umphys),intent(in) ::                                           &
    dust_div2(tdims_s%i_start:tdims_s%i_end,                                   &
              tdims_s%j_start:tdims_s%j_end,                                   &
              tdims_s%k_start:tdims_s%k_end)
               !dust mmr in div2
real(kind=real_umphys),intent(in) ::                                           &
    dust_div3(tdims_s%i_start:tdims_s%i_end,                                   &
              tdims_s%j_start:tdims_s%j_end,                                   &
              tdims_s%k_start:tdims_s%k_end)
               !dust mmr in div3
real(kind=real_umphys),intent(in) ::                                           &
    dust_div4(tdims_s%i_start:tdims_s%i_end,                                   &
              tdims_s%j_start:tdims_s%j_end,                                   &
              tdims_s%k_start:tdims_s%k_end)
               !dust mmr in div4
real(kind=real_umphys),intent(in) ::                                           &
    dust_div5(tdims_s%i_start:tdims_s%i_end,                                   &
              tdims_s%j_start:tdims_s%j_end,                                   &
              tdims_s%k_start:tdims_s%k_end)
               !dust mmr in div5
real(kind=real_umphys),intent(in) ::                                           &
    dust_div6(tdims_s%i_start:tdims_s%i_end,                                   &
              tdims_s%j_start:tdims_s%j_end,                                   &
              tdims_s%k_start:tdims_s%k_end)
               !dust mmr in div6

real(kind=real_umphys),intent(in) ::                                           &
    tracer_ukca_iacc(tdims_s%i_start:tdims_s%i_end,                            &
                     tdims_s%j_start:tdims_s%j_end,                            &
                     tdims_s%k_start:tdims_s%k_end)

real(kind=real_umphys),intent(in) ::                                           &
    tracer_ukca_icor(tdims_s%i_start:tdims_s%i_end,                            &
                     tdims_s%j_start:tdims_s%j_end,                            &
                     tdims_s%k_start:tdims_s%k_end)

real(kind=real_umphys),intent(in) ::                                           &
    theta    (tdims_s%i_start:tdims_s%i_end,                                   &
              tdims_s%j_start:tdims_s%j_end,                                   &
              tdims_s%k_start:tdims_s%k_end)

real(kind=real_umphys),intent(in) ::                                           &
    p_theta_levels (tdims_s%i_start:tdims_s%i_end,                             &
                    tdims_s%j_start:tdims_s%j_end,                             &
                    tdims_s%k_start:tdims_s%k_end)

real(kind=real_umphys), allocatable, intent(in out) :: dust_tot_nd(:,:,:)

real(kind=real_umphys), allocatable, intent(in out) :: tnuc_new(:,:,:)

real(kind=real_umphys),intent(in) ::                                           &
    exner_theta_levels (tdims_s%i_start:tdims_s%i_end,                         &
                        tdims_s%j_start:tdims_s%j_end,                         &
                        tdims_s%k_start:tdims_s%k_end)

real, pointer, intent(in) :: arcldust_b1(:,:,:)

real, pointer, intent(in) :: arcldust_b2(:,:,:)

real, pointer, intent(in) :: arcldust_b3(:,:,:)

real, pointer, intent(in) :: arcldust_b4(:,:,:)

real, pointer, intent(in) :: arcldust_b5(:,:,:)

real, pointer, intent(in) :: arcldust_b6(:,:,:)

!local variables
character(len=*), parameter :: RoutineName='LSP_PROGNOSTIC_TNUC'
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

real(kind=real_umphys) :: temp
!air temperature theta levels

real(kind=real_umphys) :: air_density_number_conc
real(kind=real_umphys) :: air_density
!density of air
real(kind=real_umphys) :: nd_factor2(2)
real(kind=real_umphys) :: nd_factor6(6)
! Factor for converting each dust division mass concentration
! to a number concentration

real(kind=real_umphys), parameter :: refdust = 758770.964
!This parameter can be tuned such that the new tnuc will act
!mostly on the Southern Ocean latitudes

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=================================================
!For calculating dust number density(particles/m3)
!from mmr (both 6-bin as well 2-bin dust schemes)
!=================================================
if (l_ukca_dust) then
  do k = 1, model_levels
    do j = 1, rows
      do i = 1, row_length
          ! Number mixing ratio in the um_tracers_ array expressed as
          ! the number of (aerosol) particles per molecule of air. To
          ! translate this to particles per m^3, we multiply by the
          ! number of air molecules per m^3.
        temp = theta(i,j,k)*exner_theta_levels(i,j,k)
        air_density_number_conc = p_theta_levels(i,j,k) / (boltzmann * temp)
        dust_tot_nd(i,j,k) =( tracer_ukca_iacc(i,j,k)                          &
                        + tracer_ukca_icor(i,j,k) )                            &
                        * air_density_number_conc
      end do
    end do
  end do
else if (l_use_arcldust) then !if using dust climatology
  do l = 1, 6
    nd_factor6(l) = 6.0/(pi*rhop*(drep(l)*drep(l)*drep(l)))
  end do
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP private(k, j, i, temp, air_density)                                      &
!$OMP SHARED(model_levels, rows, row_length, theta, exner_theta_levels, r,     &
!$OMP p_theta_levels, dust_tot_nd, nd_factor6, arcldust_b1, arcldust_b2,       &
!$OMP arcldust_b3, arcldust_b4, arcldust_b5, arcldust_b6)
  do k = 1, model_levels
    do j = 1, rows
      do i = 1, row_length
        temp = theta(i,j,k)*exner_theta_levels(i,j,k)
        air_density = p_theta_levels(i,j,k) / (r * temp)
        dust_tot_nd(i,j,k) = ( arcldust_b1(i,j,k) * nd_factor6(1)              &
                             + arcldust_b2(i,j,k) * nd_factor6(2)              &
                             + arcldust_b3(i,j,k) * nd_factor6(3)              &
                             + arcldust_b4(i,j,k) * nd_factor6(4)              &
                             + arcldust_b5(i,j,k) * nd_factor6(5)              &
                             + arcldust_b6(i,j,k) * nd_factor6(6))             &
                             * air_density
      end do
    end do
  end do
!$OMP end PARALLEL do
else if (i_dust==1) then  !if using prognostic dust
  if (l_twobin_dust) then  !i.e. if 2-bin dust
    do l = 1, 2
        ! For each dust bin, calculate the factor for converting
        ! its mass concentration to a number concentration
      nd_factor2(l) = 6.0/(pi*rhop*(drep(l)*drep(l)*drep(l)))
      ! mass of a single particle = 4/3 pi r^3 rho
      !                           = 4/3 pi (d/2)^3 rho
      !                           = 1/6 pi d^3 rho
    end do
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP private(k, j, i, temp, air_density)                                      &
!$OMP SHARED(model_levels, rows, row_length, theta, exner_theta_levels, r,     &
!$OMP p_theta_levels, dust_tot_nd, nd_factor2, dust_div1, dust_div2)
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
            ! air temperature on theta levels
          temp = theta(i,j,k)*exner_theta_levels(i,j,k)
          ! density of air
          air_density = p_theta_levels(i,j,k) / (r * temp)
          ! total dust number density ; sum of all 2 bins number
          ! concentrations  * air density.
          dust_tot_nd(i,j,k) = ( dust_div1(i,j,k) * nd_factor2(1)              &
                               + dust_div2(i,j,k) * nd_factor2(2))             &
                               * air_density
        end do
      end do
    end do
!$OMP end PARALLEL do
  else ! if 6-bin dust
    do l = 1, 6
      nd_factor6(l) = 6.0/(pi*rhop*(drep(l)*drep(l)*drep(l)))
    end do
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP private(k, j, i, temp, air_density)                                      &
!$OMP SHARED(model_levels, rows, row_length, theta, exner_theta_levels, r,     &
!$OMP p_theta_levels, dust_tot_nd, nd_factor6, dust_div1, dust_div2,           &
!$OMP dust_div3, dust_div4, dust_div5, dust_div6)
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          ! air temperature on theta levels
          temp = theta(i,j,k)*exner_theta_levels(i,j,k)
          ! density of air
          air_density = p_theta_levels(i,j,k) / (r * temp)
          ! total dust number density ; sum of all 6 bins number
          ! concentrations  * air density.
          dust_tot_nd(i,j,k) = ( dust_div1(i,j,k) * nd_factor6(1)              &
                               + dust_div2(i,j,k) * nd_factor6(2)              &
                               + dust_div3(i,j,k) * nd_factor6(3)              &
                               + dust_div4(i,j,k) * nd_factor6(4)              &
                               + dust_div5(i,j,k) * nd_factor6(5)              &
                               + dust_div6(i,j,k) * nd_factor6(6) )            &
                               * air_density
        end do
      end do
    end do
!$OMP end PARALLEL do
  end if ! for progn. dust
else
  ErrorStatus = 100
  comments = 'l_progn_tnuc has been turned on        '// newline//             &
       'with no dust selected.Please open the GUI and choose   '// newline//   &
       'i_dust = prognostic or l_use_arcldust = true'
  call ereport(RoutineName, ErrorStatus, comments)
end if ! climatology dust

#if defined(CRAYFTN_VERSION) && (CRAYFTN_VERSION < 15000000)
! The following calculation breaks bit-reproducibility between different
! processor decompositions if vectorised by Cray fortran for some reason.
! Disable vectorisation when compiling with older versions of Cray fortran
!DIR$ NOVECTOR
#endif
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP private(k, j, i)                                                         &
!$OMP SHARED(model_levels, rows, row_length, dust_tot_nd, tnuc_new)
do k = 1, model_levels
  do j = 1, rows
    do i = 1, row_length
      !===========================================================
      ! Heterogeneous nucleation temperature as a function of dust
      !===========================================================
      if (dust_tot_nd(i,j,k) > zero) then
        tnuc_new(i,j,k) = thomo + (tnuc - thomo) *                             &
                        (atan(5.0*log10(dust_tot_nd(i,j,k)/refdust))/pi + 0.5)
      else
        ! If there is no dust then we set tnuc_new to the homogenous
        ! freezing limit
        tnuc_new(i,j,k) = thomo
      end if
    end do
  end do
end do
!$OMP end PARALLEL do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_prognostic_tnuc
end module lsp_prognostic_tnuc_mod
