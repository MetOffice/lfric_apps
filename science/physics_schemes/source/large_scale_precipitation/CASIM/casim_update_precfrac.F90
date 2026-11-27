! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine casim_update_precfrac
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

module casim_update_precfrac_mod

implicit none

character(len=*), parameter, private :: ModuleName='CASIM_UPDATE_PRECFRAC_MOD'

contains

! Subroutine to update the prognostic precip fraction consistent with the
! rain and graupel mass increments after calling CASIM
subroutine casim_update_precfrac( nlayers, dz_casim, rho_casim,                &
                                  rainfall_3d, graupfall_3d,                   &
                                  qc_casim, qi_casim, qs_casim,                &
                                  qr_casim, qg_casim,                          &
                                  cfliq_casim, cfice_casim, precfrac_casim,    &
                                  dqc_casim, dqi_casim, dqs_casim,             &
                                  dqr_casim, dqg_casim,                        &
                                  dcfliq_casim, dcfice_casim )

use variable_precision, only: wp
use timestep_mod, only: timestep

implicit none

! Number of model-levels
integer, intent(in) :: nlayers

! Model-level thicknesses and dry-mass density
real(kind=wp), dimension(nlayers), intent(in) :: dz_casim, rho_casim

! Fall-fluxes of rain and graupel / kg m-2 s-1
real(kind=wp), dimension(nlayers), intent(in) :: rainfall_3d, graupfall_3d

! Condensate species mixing-ratios on input to CASIM
real(kind=wp), dimension(nlayers), intent(in) :: qc_casim, qi_casim, qs_casim, &
                                                 qr_casim, qg_casim

! Ice and liquid cloud-fractions
real(kind=wp), dimension(nlayers), intent(in) :: cfliq_casim, cfice_casim

! Prognostic precip fraction to be updated
real(kind=wp), dimension(nlayers), intent(inout) :: precfrac_casim

! CASIM increments to condensate species mixing-ratios
real(kind=wp), dimension(nlayers), intent(in) :: dqc_casim,dqi_casim,dqs_casim,&
                                                 dqr_casim, dqg_casim

! Increments to liquid and ice cloud-fractions
real(kind=wp), dimension(nlayers), intent(in) :: dcfliq_casim, dcfice_casim

! Layer-mass = rho * dz
real(kind=wp), dimension(nlayers) :: rhodz

! Precip masses:
real(kind=wp) :: prec_k     ! Pre-existing at level k
real(kind=wp) :: prec_k_f   ! After fall-in mass added
real(kind=wp) :: prec_k_f_c ! After fall-in and cloud sources added
real(kind=wp) :: prec_cl    ! Source from liquid-cloud (autoconversn)
real(kind=wp) :: prec_cf    ! Source from ice-cloud (melting)
real(kind=wp) :: prec_accl  ! Source from liquid-cloud (accretion)
real(kind=wp) :: prec_accf  ! Source from ice-cloud (riming)
real(kind=wp) :: prec_fall  ! Mass falling in from above

! Precip fraction at k
real(kind=wp) :: precfrac_k_f    ! After fall-in mass added
real(kind=wp) :: precfrac_k_f_c  ! After fall-in mass and cloud added

! Min limit on fractions, for safety
real(kind=wp) :: min_frac

! Fall-speed * dt / dz
real(kind=wp) :: v_dt_rdz
! Estimated end-of-timestep precip mixing-ratio for comparison with qcl, qcf
real(kind=wp) :: q_precip
! Fraction of cloud-sources from accretion vs autoconversion
real(kind=wp) :: accfac
! Fraction of cloud-sources from liquid vs ice cloud
real(kind=wp) :: liqfac
! Mid-point cloud condensate mixing-ratios used to compute liqfac
real(kind=wp) :: qc_tmp
real(kind=wp) :: qi_tmp
real(kind=wp) :: qs_tmp
! Liquid and ice cloud fractions used in precip fraction update
real(kind=wp) :: cfl_src
real(kind=wp) :: cff_src

! Loop counter
integer :: k

! 0.0, 0.5, 1.0 in native precision
real(kind=wp), parameter :: zero = 0.0_wp
real(kind=wp), parameter :: half = 0.5_wp
real(kind=wp), parameter :: one  = 1.0_wp

! Miniscule number for check to avoid div-by-zero
real(kind=wp), parameter :: min_float = tiny(prec_k)
real(kind=wp), parameter :: sqrt_min_float = sqrt(tiny(prec_k))


! Precompute layer-mass = rho * dz
do k = 1, nlayers
  rhodz(k) = rho_casim(k) * dz_casim(k)
end do

! Perform a downwards sweep, moving the precip fraction downwards with
! the precip mass...
do k = nlayers-1, 1, -1

  ! 1) Update level k precfrac due to flux falling-in from above

  ! Set min limit on the fractions to be combined
  !  = 0.01 times the largest source fraction.
  ! This is to avoid the fraction going stupidly-small due to instances
  ! of tiny amounts of precip or cloud mass with no fraction, e.g.
  ! due to numerical noise in the transport scheme.
  min_frac = max( 0.01 * max( precfrac_casim(k), precfrac_casim(k+1) ),        &
                  min_float )

  ! Pre-existing mass of "precip" (ignore negative values)
  prec_k = rhodz(k) * ( max(qr_casim(k),zero)                                  &
                      + max(qg_casim(k),zero) )
  ! Mass of precip falling-in from above (ignore negative values)
  prec_fall = timestep * ( max(rainfall_3d(k+1),zero)                          &
                         + max(graupfall_3d(k+1),zero) )
  ! Add-on mass of precip falling from above
  prec_k_f = prec_k + prec_fall
  ! Compute updated precfrac after combining with fall from above:
  ! 1/sqrt(frac) = sum( m/sqrt(frac) ) / sum( m )
  ! => frac = ( sum( m ) / sum( m/sqrt(frac) ) )**2
  precfrac_k_f = ( prec_k_f / max(                                             &
       prec_k    / sqrt( max(precfrac_casim(k),   min_frac) )                  &
     + prec_fall / sqrt( max(precfrac_casim(k+1), min_frac) ),                 &
                                   min_float ) )**2

  ! 2) Update level k precfrac due to precip created by cloud processes

  ! Compute net source of rain + graupel from CASIM on level k,
  ! accounting for the divergence of the fall-flux
  prec_cl = max( rhodz(k)*( dqr_casim(k) + dqg_casim(k) )                      &
               + timestep*( rainfall_3d(k)  - rainfall_3d(k+1)                 &
                          + graupfall_3d(k) - graupfall_3d(k+1) ),             &
                 zero )   ! Note: ignoring sinks of precip for now
  ! Add cloud sources onto precip mass at k
  prec_k_f_c = prec_k_f + prec_cl

  ! Partition between liquid and ice cloud sources as a function of
  ! the cloud condensate masses.
  ! Taking centred discretisation for the cloud masses:
  qc_tmp = max( qc_casim(k) + half*dqc_casim(k), zero )
  qi_tmp = max( qi_casim(k) + half*dqi_casim(k), zero )
  qs_tmp = max( qs_casim(k) + half*dqs_casim(k), zero )
  liqfac = qc_tmp / max( qc_tmp + qi_tmp + qs_tmp, min_float )
  prec_cf = (one-liqfac) * prec_cl
  prec_cl = liqfac       * prec_cl

  ! Partition precip produced by cloud into contributions from
  ! autoconversion / melting (whose area is the whole cloud-fraction)
  ! vs contributions from accretion / riming (whose area is just the
  ! overlap between the cloud and the existing precip fraction).
  ! Partition as a function of the ratio of cloud-mass to precip-mass,
  ! so that:
  ! if q_cloud >> q_precip : assume autoconversion / melting dominate
  ! if q_cloud << q_precip : assume accretion / riming dominate
  !
  ! First estimate q_precip after fall-in and fall-out only:
  ! (assume fall-out gives decay by fraction 1/(1 + v dt/dz))
  ! fall_flux = v rho q_eqm
  ! => v dt/dz = fall_flux dt / ( q_eqm rho dz)
  ! Take centred discretization for q_eqm:
  v_dt_rdz = timestep * ( max(rainfall_3d(k),zero)                             &
                        + max(graupfall_3d(k),zero) )                          &
           / max( prec_k + rhodz(k) * half                                     &
                           * ( dqr_casim(k) + dqg_casim(k) ),                  &
                  sqrt_min_float )
  q_precip = prec_k_f / ( rhodz(k) * ( one + v_dt_rdz ) )
  ! Fraction of mixing-ratio that is "precip" vs "cloud"
  accfac = q_precip / max( max(qc_casim(k),zero) + max(qi_casim(k),zero)       &
                         + max(qs_casim(k),zero) + q_precip, min_float )
  ! Partition sources of precip from cloud in proportion
  prec_accl = accfac * prec_cl
  prec_accf = accfac * prec_cf
  prec_cl = (one-accfac) * prec_cl
  prec_cf = (one-accfac) * prec_cf

  ! Compute updated precfrac after combining with cloud-sources:
  ! 1/sqrt(frac) = sum( m/sqrt(frac) ) / sum( m )
  ! => frac = ( sum( m ) / sum( m/sqrt(frac) ) )**2
  ! Take centred discretisation for the cloud-fractions; especially need
  ! to account for increase of cff by ice-fall.
  cfl_src = cfliq_casim(k) + half*dcfliq_casim(k)
  cff_src = cfice_casim(k) + half*dcfice_casim(k)
  min_frac = max( min_frac, 0.01*max( cfl_src, cff_src ) )
  precfrac_casim(k) = ( prec_k_f_c / max(                                      &
       prec_k_f  / sqrt( max(precfrac_k_f,               min_frac) )           &
     + prec_cl   / sqrt( max(cfl_src,                    min_frac) )           &
     + prec_cf   / sqrt( max(cff_src,                    min_frac) )           &
     + prec_accl / sqrt( max(min(cfl_src, precfrac_k_f), min_frac) )           &
     + prec_accf / sqrt( max(min(cff_src, precfrac_k_f), min_frac) ),          &
                                       min_float ) )**2
  ! Assuming accretion / riming occur in area with max overlap between
  ! cloud-fraction and precip-fraction = min(cf, precfrac)

end do  ! k = nlayers-1, 1, -1

! Final check; reset precfrac to zero where all precip has evaporated
do k = 1, nlayers
  if ( .not. ( qr_casim(k) + dqr_casim(k) > zero .or.                          &
               qg_casim(k) + dqg_casim(k) > zero ) ) then
    precfrac_casim(k) = zero
  end if
end do


return
end subroutine casim_update_precfrac

end module casim_update_precfrac_mod
