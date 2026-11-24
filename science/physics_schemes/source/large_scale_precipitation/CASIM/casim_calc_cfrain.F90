! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine casim_calc_cfrain
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

module casim_calc_cfrain_mod

implicit none

contains

! Subroutine to make a pre-estimate of the updated rain / graupel fraction to
! pass into CASIM, when using the prognostic precip fraction.
! We use the same method as is used to update the precip fraction after
! CASIM (subroutine casim_update_precfrac), but have to make additional
! approximations in place of the terms which haven't been calculated yet
! (the precip fall-speed / flux, and the source-terms from cloud).
subroutine casim_calc_cfrain( nlayers, dz_casim, rho_casim,                    &
                              qc_casim, qi_casim, qs_casim, qr_casim, qg_casim,&
                              cfliq_casim, cfice_casim, precfrac_casim,        &
                              cfrain_casim, cfgr_casim )

use variable_precision, only: wp
use timestep_mod, only: timestep

implicit none

! Number of model-levels
integer, intent(in) :: nlayers

! Model-level thicknesses and dry-mass density
real(kind=wp), dimension(nlayers), intent(in) :: dz_casim, rho_casim

! Condensate species mixing-ratios on input to CASIM
real(kind=wp), dimension(nlayers), intent(in) :: qc_casim, qi_casim, qs_casim, &
                                                 qr_casim, qg_casim

! Ice and liquid cloud-fractions, and prognostic precip fraction
real(kind=wp), dimension(nlayers), intent(in) :: cfliq_casim, cfice_casim,     &
                                                 precfrac_casim

! Rain and graupel fractions to be calculated
real(kind=wp), dimension(nlayers), intent(out) :: cfrain_casim, cfgr_casim

! Layer-mass = rho * dz
real(kind=wp), dimension(nlayers) :: rhodz

! Precip mass falling from above
real(kind=wp) :: prec_fall

! Precip masses:
real(kind=wp) :: prec_k     ! Pre-existing at level k
real(kind=wp) :: prec_k_f   ! After fall-in mass added
real(kind=wp) :: prec_k_f_c ! After fall-in and cloud sources added
real(kind=wp) :: prec_cl    ! Source from liquid-cloud (autoconversn)
real(kind=wp) :: prec_cf    ! Source from ice-cloud (melting)
real(kind=wp) :: prec_accl  ! Source from liquid-cloud (accretion)
real(kind=wp) :: prec_accf  ! Source from ice-cloud (riming)

! Precip fraction at k after fall-in mass added
real(kind=wp) :: precfrac_k_f

! Min limit on fractions, for safety
real(kind=wp) :: min_frac

! Timestep / time-scales for decay of liquid and ice into "precip"
real(kind=wp) :: dt_rtau_cl
real(kind=wp) :: dt_rtau_cf

! Fall-distance of precip over a timestep = fall-speed * dt
real(kind=wp) :: v_dt
! Fall-speed * dt / dz
real(kind=wp) :: v_dt_rdz
! Estimated end-of-timestep precip mixing-ratio for comparison with qcl, qcf
real(kind=wp) :: q_precip
! Fraction of cloud-sources from accretion vs autoconversion
real(kind=wp) :: accfac

! Loop counter
integer :: k

! 0.0, 1.0 in native precision
real(kind=wp), parameter :: zero = 0.0_wp
real(kind=wp), parameter :: one  = 1.0_wp

! Miniscule number for check to avoid div-by-zero
real(kind=wp), parameter :: min_float = tiny(prec_k)


! Set decay factors dt/tau for liquid and ice cloud into precip
dt_rtau_cl = timestep / 7200.0  ! 2 hours
dt_rtau_cf = timestep / 7200.0  ! 2 hours

! Set fall-distance of precip over current timestep
v_dt = timestep * 5.0  ! 5 ms-1

! Precompute layer-mass = rho * dz
do k = 1, nlayers
  rhodz(k) = rho_casim(k) * dz_casim(k)
end do

! Initialise fall-flux to zero
prec_fall = zero
! Initialise updated precip fraction to zero at the model-lid
cfrain_casim(nlayers) = zero

! Vertically-integrate downwards...
do k = nlayers-1, 1, -1

  ! 1) Update level k precfrac due to flux falling-in from above

  ! Set min limit on the fractions to be combined
  !  = 0.01 times the largest source fraction.
  ! This is to avoid the fraction going stupidly-small due to instances
  ! of tiny amounts of precip or cloud mass with no fraction, e.g.
  ! due to numerical noise in the transport scheme.
  min_frac = max( 0.01 * max( precfrac_casim(k), cfrain_casim(k+1) ),          &
                  min_float )

  ! Pre-existing mass of "precip" (ignore negative values)
  prec_k = rhodz(k) * ( max(qr_casim(k),zero)                                  &
                      + max(qg_casim(k),zero) )
  ! Add-on mass of precip falling from above
  prec_k_f = prec_k + prec_fall
  ! Compute updated precfrac after combining with fall from above:
  ! 1/sqrt(frac) = sum( m/sqrt(frac) ) / sum( m )
  ! => frac = ( sum( m ) / sum( m/sqrt(frac) ) )**2
  precfrac_k_f = ( prec_k_f / max(                                             &
                prec_k    / sqrt( max(precfrac_casim(k), min_frac) )           &
              + prec_fall / sqrt( max(cfrain_casim(k+1), min_frac) ),          &
                                   min_float ) )**2

  ! 2) Update level k precfrac due to precip created by cloud processes

  ! Total guess for the mass of precip that might be produced by the cloud;
  ! Model as decay of cloud into precip over a fixed time-scale
  prec_cl = rhodz(k) * dt_rtau_cl * max(qc_casim(k),zero)
  prec_cf = rhodz(k) * dt_rtau_cf * ( max(qi_casim(k),zero)                    &
                                    + max(qs_casim(k),zero) )
  ! Add cloud sources onto precip mass at k
  prec_k_f_c = prec_k_f + prec_cl + prec_cf

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
  v_dt_rdz = v_dt/dz_casim(k)
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
  min_frac = max( min_frac, 0.01*max( cfliq_casim(k), cfice_casim(k) ) )
  cfrain_casim(k) = ( prec_k_f_c / max(                                        &
       prec_k_f  / sqrt( max(precfrac_k_f,                      min_frac) )    &
     + prec_cl   / sqrt( max(cfliq_casim(k),                    min_frac) )    &
     + prec_cf   / sqrt( max(cfice_casim(k),                    min_frac) )    &
     + prec_accl / sqrt( max(min(cfliq_casim(k), precfrac_k_f), min_frac) )    &
     + prec_accf / sqrt( max(min(cfice_casim(k), precfrac_k_f), min_frac) ),   &
                                            min_float ) )**2
  ! Assuming accretion / riming occur in area with max overlap between
  ! cloud-fraction and precip-fraction = min(cf, precfrac)

  ! 3) Set fall-flux passed down to next level

  ! The full mass of precip passing into the current level is
  ! partitioned into the fraction remaining on this level:
  !   1 / ( 1 + v dt/dz )
  ! and the fraction falling down to the next level:
  !   1 - 1 / ( 1 + v dt/dz )
  ! = v dt/dz / ( 1 + v dt/dz )
  prec_fall = prec_k_f_c * v_dt_rdz / ( one + v_dt_rdz )

end do  ! k = nlayers-1, 1, -1

! Set graupel fraction = rain fraction
do k = 1, nlayers
  cfgr_casim(k) = cfrain_casim(k)
end do


return
end subroutine casim_calc_cfrain

end module casim_calc_cfrain_mod
