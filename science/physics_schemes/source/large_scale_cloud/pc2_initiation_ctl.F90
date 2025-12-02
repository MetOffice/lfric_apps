! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initiate cloud and liquid water within the PC2 Cloud Scheme.

module pc2_initiation_ctl_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_INITIATION_CTL_MOD'
contains

subroutine pc2_initiation_ctl (                                                &
! Dimensions of Rh crit array
  rhc_row_length, rhc_rows, zlcl_mixed, r_theta_levels,                        &

! Model switches
  l_mixing_ratio, l_wtrac,                                                     &

! SCM diagnostics switches
  nSCMDpkgs,L_SCMDiags,                                                        &

! Primary fields passed in/out
  t,q,qcl,qcf,qcf2,cf,cfl,cff,rhts,tlts,qtts,ptts,cf_area,                     &

! Primary fields passed in
  p,pstar,p_theta_levels,cumulus,rhcrit,                                       &

! Input Bimodal fields
  tgrad_bm,bl_w_var,tau_dec_bm,tau_hom_bm,tau_mph_bm,z_theta,                  &
  ri_bm, mix_len_bm, zh,zhsc,dzh,bl_type_7,                                    &

! Output increments
  calculate_increments,                                                        &
  t_work,q_work,qcl_work,qcf_work,qcf2_work,cf_work,cfl_work,cff_work,         &

! Output bimodal fields
  sskew, svar_turb, svar_bm,                                                   &

! Water tracer structure
  wtrac                                                                        &
 )

use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims, pdims_s, tdims_l


use cloud_inputs_mod,      only: i_cld_area, i_pc2_init_method, cloud_pc2_tol
use pc2_constants_mod,     only: acf_cusack, pc2init_bimodal,                  &
                                 cloud_rounding_tol

use s_scmop_mod,           only: default_streams,                              &
                                 t_inst, d_wet, d_all, scmdiag_pc2

use model_domain_mod, only: model_type, mt_single_column

use pc2_arcld_mod, only: pc2_arcld
use pc2_checks_mod, only: pc2_checks
use pc2_checks2_mod, only: pc2_checks2
use pc2_hom_arcld_mod, only: pc2_hom_arcld
use pc2_initiate_mod, only: pc2_initiate
use pc2_bm_initiate_mod, only: pc2_bm_initiate
use mphys_inputs_mod,   only: l_mcr_qcf2

use free_tracers_inputs_mod, only: n_wtrac
use water_tracers_mod,       only: wtrac_type
use wtrac_pc2_mod,           only: wtrac_pc2_store
use wtrac_pc2_phase_chg_mod, only: wtrac_pc2_phase_chg

implicit none

! Description:
!   Initiate a small amount of cloud fraction and liquid
!   water content for the PC2 Cloud Scheme. Check that moisture
!   variables are consistent with each other.
!
! Method:
!   See the PC2 documentation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Code Description:
!   Language: fortran 77 + common extensions
!   This code is written to UMDP3 v6 programming standards.

! Declarations:
! Arguments with intent in. ie: input variables.

! Model dimensions
integer ::                                                                     &
  rhc_row_length,  & ! Dimensions of RHcrit variable
  rhc_rows           ! Dimensions of RHcrit variable
!
logical ::                                                                     &
  cumulus(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end), & ! Convection is occurring.
  l_mixing_ratio                        ! Use mixing ratio formulation

logical :: l_wtrac   ! Controls water tracer calculations

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
real(kind=real_umphys) :: zlcl_mixed( pdims%i_start:pdims%i_end,               &
                                      pdims%j_start:pdims%j_end)

! Primary fields passed in/out
real(kind=real_umphys) ::                                                      &
  t(                 tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  q(                 tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  qcl(               tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  qcf(               tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  qcf2(              tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  cf(                tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  cfl(               tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  cff(               tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  cf_area(           tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end)

!cloud ice + snow
real(kind=real_umphys) ::  qcf_total(tdims%i_start:tdims%i_end,                &
                                tdims%j_start:tdims%j_end,                     &
                                            1:tdims%k_end)

! Primary fields passed in
real(kind=real_umphys) ::                                                      &
  p(                 pdims_s%i_start:pdims_s%i_end,                            &
                     pdims_s%j_start:pdims_s%j_end,                            &
                     pdims_s%k_start:pdims_s%k_end+1),                         &
  pstar(             pdims%i_start  :pdims%i_end,                              &
                     pdims%j_start  :pdims%j_end),                             &
  p_theta_levels(    pdims%i_start:pdims%i_end,                                &
                     pdims%j_start:pdims%j_end,                                &
                     pdims%k_start:pdims%k_end),                               &
  rhcrit(            rhc_row_length,                                           &
                     rhc_rows,                                                 &
                                 1:tdims%k_end),                               &
  rhts(              tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  tlts(              tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       TL at start of timestep
    qtts(              tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end),                             &
!       qT at start of timestep
    ptts(              pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end,                              &
                       pdims%k_start:pdims%k_end),                             &
!       Pressure at theta levels at start of timestep
  bl_w_var(          tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Vertical velocity variance
  z_theta(           tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Height of levels
  r_theta_levels(    tdims_l%i_start:tdims_l%i_end,                            &
                     tdims_l%j_start:tdims_l%j_end,                            &
                     tdims_l%k_start:tdims_l%k_end),                           &
!       Height of levels
  tgrad_bm(          tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Gradient of liquid potential temperature
  ri_bm(             tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Richardson Number for bimodal cloud scheme
  mix_len_bm(        tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Turbulent mixing length for bimodal cloud scheme
  zh(                tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end),                               &
!       Boundary-layer height for bimodal cloud scheme
  zhsc(              tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end),                               &
!       Decoupled layer height for bimodal cloud scheme
  dzh(               tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end),                               &
!       Inversion Thickness for bimodal cloud scheme
  bl_type_7(         tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end),                               &
!       Shear-driven boundary layer indicator
  tau_dec_bm(        tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Decorrelation time scale
  tau_mph_bm(        tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Phase-relaxation time scale
  tau_hom_bm(        tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end)
!       Turbulence homogenisation time scale

integer ::                                                                     &
  nSCMDpkgs             ! No of SCM diagnostics packages

logical ::                                                                     &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

logical, intent(in) :: calculate_increments
! are the increments required for e.g. diagnostics?

! Water tracer structure.
! (Note, this routine can be called from pc2_pressure_forcing which updates
! the q fields and from cloud_call_b4_conv which updates the q_star fields.
! For water tracers, this routine currently only works for the water tracer q
! fields. Hence, it is not possible to run with water tracers if
! l_cloud_call_b4_conv = T.)
type(wtrac_type), intent(in out) :: wtrac(n_wtrac)

! Local variables

! Output increment diagnostics
real(kind=real_umphys) ::                                                      &
  t_work(            tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  q_work(            tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  qcl_work(          tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  qcf_work(          tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  qcf2_work(         tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  cf_work(           tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  cfl_work(          tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  cff_work(          tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
  p_layer_boundaries(pdims%i_start:pdims%i_end,                                &
                     pdims%j_start:pdims%j_end,                                &
                                 0:pdims%k_end),                               &
!       Pressure at layer boundaries. Same as p except at
!       bottom level = pstar, and at top = 0.
  p_layer_centres(   pdims%i_start:pdims%i_end,                                &
                     pdims%j_start:pdims%j_end,                                &
                                 0:pdims%k_end),                               &
!       Pressure at layer centres. Same as p_theta_levels
!       except bottom level = pstar, and at top = 0.
  sskew(             tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Skewness of mixture of two s-distributions on levels
  svar_turb(         tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end),                               &
!       Variance of turbulence-based uni-modal s-distribution on levels
  svar_bm(           tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end,                                &
                                 1:tdims%k_end)
!       Variance of mixture of two s-distributions in the bi-modal scheme

integer ::                                                                     &
  i,j,k,iScm,jScm,                                                             &
!       Loop counters
    large_levels,                                                              &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((model_levels - 2)*levels_per_level) + 2
    levels_per_level
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops


character(len=*), parameter ::  RoutineName = 'PC2_INITIATION_CTL'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

! External Functions:

!- End of header

! Extra diagnostics from bimodal initiation (only output for SCM).
real(kind=real_umphys), allocatable :: entzone(:,:,:)
real(kind=real_umphys), allocatable :: sl_modes(:,:,:,:)
real(kind=real_umphys), allocatable :: qw_modes(:,:,:,:)
real(kind=real_umphys), allocatable :: rh_modes(:,:,:,:)
real(kind=real_umphys), allocatable :: sd_modes(:,:,:,:)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (calculate_increments) then

!$OMP PARALLEL DEFAULT(SHARED) private(i,j,k)

!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

        ! Work fields are set to starting fields so diagnostics
        ! can be calculated

        q_work(i,j,k)   = q(i,j,k)
        qcl_work(i,j,k) = qcl(i,j,k)
        qcf_work(i,j,k) = qcf(i,j,k)
        cf_work(i,j,k)  = cf(i,j,k)
        cfl_work(i,j,k) = cfl(i,j,k)
        cff_work(i,j,k) = cff(i,j,k)
        if ( l_mcr_qcf2 ) then
          qcf2_work(i,j,k) = qcf2(i,j,k)
        end if

      end do
    end do
  end do
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        t_work(i,j,k)   = t(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
!$OMP end PARALLEL
end if

! Call checking routine
! Pass field arrays without halo cells.
call pc2_checks(p_theta_levels,                                                &
    t, cf, cfl, cff, q, qcl, qcf,                                              &
    l_mixing_ratio,                                                            &
    tdims%i_len, tdims%j_len, tdims%k_end,                                     &
    tdims%halo_i, tdims%halo_j, tdims%halo_i, tdims%halo_j, qcf2, wtrac)

!-----------------------------------------------------------------------
!       SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------

if ( i_cld_area == acf_cusack ) then
  ! Setup variables used by the area cloud fraction routines

  ! set p at layer boundaries.
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      p_layer_boundaries(i,j,0) = pstar(i,j)
      p_layer_centres(i,j,0) = pstar(i,j)
    end do
  end do
  do k = pdims%k_start, pdims%k_end - 1
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        p_layer_boundaries(i,j,k) = p(i,j,k+1)
        p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
      end do
    end do
  end do
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      p_layer_boundaries(i,j,pdims%k_end) = 0.0
      p_layer_centres(i,j,pdims%k_end) =                                       &
                       p_theta_levels(i,j,pdims%k_end)
    end do
  end do

  ! Determine number of sublevels for vertical gradient area cloud
  ! Want an odd number of sublevels per level: 3 is hardwired in do loops
  levels_per_level = 3
  large_levels = ((pdims%k_end - 2)*levels_per_level) + 2

end if  ! ( i_cld_area == acf_cusack )

if ( calculate_increments .and. i_pc2_init_method == pc2init_bimodal ) then
  ! Allocate SCM diagnostics for bimodal cloud-scheme
  allocate( entzone  ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                       1:tdims%k_end ) )
  allocate( sl_modes ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                       1:tdims%k_end, 3 ) )
  allocate( qw_modes ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                       1:tdims%k_end, 3 ) )
  allocate( rh_modes ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                       1:tdims%k_end, 3 ) )
  allocate( sd_modes ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                       1:tdims%k_end, 3 ) )
else
  allocate( entzone(1,1,1) )
  allocate( sl_modes(1,1,1,1) )
  allocate( qw_modes(1,1,1,1) )
  allocate( rh_modes(1,1,1,1) )
  allocate( sd_modes(1,1,1,1) )
end if

! Store water values prior to initiation routine call
if (l_wtrac) call wtrac_pc2_store(q, qcl)

! Call initiation routine

if (i_pc2_init_method == pc2init_bimodal) then

  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qcf_total(i,j,k) = qcf(i,j,k)

        if (l_mcr_qcf2) then
          !make a total qcf to pass to bimodal cloud scheme.
          qcf_total(i,j,k) = qcf_total(i,j,k) + qcf2(i,j,k)
        end if

      end do
    end do
  end do

  call pc2_bm_initiate(p_theta_levels,cumulus,tgrad_bm,bl_w_var,               &
      tau_dec_bm,tau_hom_bm,tau_mph_bm,ri_bm, mix_len_bm,                      &
      zh,zhsc,dzh,bl_type_7,                                                   &
      tdims%k_end,zlcl_mixed,r_theta_levels,z_theta,t,cf,cfl,cff,              &
      q,qcl,qcf_total,sskew,svar_turb,svar_bm,entzone,                         &
      sl_modes, qw_modes, rh_modes, sd_modes,                                  &
      calculate_increments, l_mixing_ratio)


else

  ! using area cloud parameterisation
  if (i_cld_area == acf_cusack) then

    call pc2_arcld(p_layer_centres,p_layer_boundaries,                         &
     cumulus,rhcrit,                                                           &
     rhc_row_length,rhc_rows,zlcl_mixed,                                       &
     large_levels,levels_per_level,cf_area,                                    &
     t,cf,cfl,cff,q,qcl,qcf,rhts,tlts,qtts,ptts,l_mixing_ratio)

  else !i_cld_area

    call pc2_initiate(p_theta_levels,cumulus,rhcrit,                           &
      tdims%k_end, rhc_row_length,rhc_rows,zlcl_mixed,r_theta_levels,          &
      t,cf,cfl,cff,q,qcl,rhts,l_mixing_ratio)

  end if !i_cld_area

end if

! Update water tracers for initiation
if (l_wtrac) then
  call wtrac_pc2_phase_chg(tdims, q, qcl, 'pc2_initiate', wtrac = wtrac)
end if

!CASIM not doing scm at the moment
!-----------------------------------------------------------------------
!       SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
! Update initinc arrays to hold net increment from initiation


! Call second checking routine
! Optionally tidy-up small liquid cloud fractions.  Note that the subsequent
! call to pc2_checks similarly tidies up with a different threshold
! cloud_rounding_tol.  So no need to call pc2_checks2 unless the threshold
! cloud-fraction used is higher than cloud_rounding_tol.
if ( cloud_pc2_tol > cloud_rounding_tol ) then

  if (l_wtrac) call wtrac_pc2_store(q, qcl)

  call pc2_checks2(p_theta_levels,rhcrit,                                      &
      rhc_row_length,rhc_rows,                                                 &
      t, cf, cfl, cff, q, qcl, l_mixing_ratio)

  if (l_wtrac) then
    call wtrac_pc2_phase_chg(tdims, q, qcl, 'pc2_checks2', wtrac = wtrac)
  end if

end if

! Call first checking routine again
! Pass field arrays without halo cells.
call pc2_checks(p_theta_levels,                                                &
    t, cf, cfl, cff, q, qcl, qcf,                                              &
    l_mixing_ratio,                                                            &
    tdims%i_len, tdims%j_len, tdims%k_end,                                     &
    tdims%halo_i, tdims%halo_j, tdims%halo_i, tdims%halo_j, qcf2, wtrac)

! using area cloud parameterisation
if (i_cld_area == acf_cusack) then
  call pc2_hom_arcld(p_layer_centres,p_layer_boundaries,                       &
     large_levels,levels_per_level,                                            &
     cf_area,t,cf,cfl,cff,q,qcl,qcf,                                           &
     l_mixing_ratio)
end if    ! i_cld_area

if (calculate_increments) then
  ! Update work array to hold net increment from the above routines

!$OMP PARALLEL DEFAULT(SHARED) private(i,j,k)

!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        q_work(i,j,k)   = q(i,j,k)   - q_work(i,j,k)
        qcl_work(i,j,k) = qcl(i,j,k) - qcl_work(i,j,k)
        qcf_work(i,j,k) = qcf(i,j,k) - qcf_work(i,j,k)
        cf_work(i,j,k)  = cf(i,j,k)  - cf_work(i,j,k)
        cfl_work(i,j,k) = cfl(i,j,k) - cfl_work(i,j,k)
        cff_work(i,j,k) = cff(i,j,k) - cff_work(i,j,k)
        if (l_mcr_qcf2) then
          qcf2_work(i,j,k) = qcf2(i,j,k) - qcf2_work(i,j,k)
        end if
      end do
    end do
  end do
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        t_work(i,j,k)   = t(i,j,k)   - t_work(i,j,k)
      end do
    end do
  end do
!$OMP end do
!$OMP end PARALLEL
end if

!-----------------------------------------------------------------------
!       SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------

! End of routine initiation_ctl

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_initiation_ctl
end module pc2_initiation_ctl_mod
