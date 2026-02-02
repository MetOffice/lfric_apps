! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_mean_q_cl_with_core_mod

implicit none

contains

! Subroutine to impose a minimum allowed liquid cloud mixing-ratio
! in the parcel mean properties, due to the mean encompassing the
! parcel core (in the case where the core contains nonzero liquid-cloud,
! the mean must contain liquid cloud as well).
! When the parcel mean q_cl is inflated in this way, the parcel mean
! properties will become subsaturated.  This is due to in-plume
! moisture variability making part of the parcel saturated even
! though its mean properties are subsaturated.
subroutine calc_mean_q_cl_with_core( n_points, n_points_super, n_fields_tot,   &
                                     core_q_cl,                                &
                                     core_mean_ratio, par_mean_fields )

use comorph_constants_mod, only: real_cvprec, zero, one, sqrt_min_delta,       &
                                 l_cv_cloudfrac, i_mean_q_cl,                  &
                                 i_mean_q_cl_full
use fields_type_mod, only: i_temperature, i_q_vap, i_q_cl, i_q_cf,             &
                           i_qc_first, i_qc_last, i_cf_liq, i_cf_bulk

use set_cp_tot_mod, only: set_cp_tot
use lat_heat_mod, only: lat_heat_incr, i_phase_change_con
use set_par_cloudfrac_mod, only: set_par_cloudfrac

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of the primary fields super-arrays
integer, intent(in) :: n_points_super

! Total number of primary fields
integer, intent(in) :: n_fields_tot

! Liquid cloud mixing ratio of the parcel core
real(kind=real_cvprec), intent(in) :: core_q_cl(n_points)

! Ratio of core buoyancy / mean buoyancy
! (sets the power in the assumed in-parcel power-law pdf)
real(kind=real_cvprec), intent(in) :: core_mean_ratio(n_points)

! Parcel mean fields
real(kind=real_cvprec), intent(in out) :: par_mean_fields                      &
                                          ( n_points_super, n_fields_tot )

! Condensation increment due to sub-plume moisture variability
real(kind=real_cvprec) :: dq_cond(n_points)
! Total heat capacity for phase-change
real(kind=real_cvprec) :: cp_tot(n_points)

! List of indices of points where calculations needed
integer :: nc
integer :: index_ic(n_points)

! Loop counter
integer :: ic, ic2

! Value very slightly below 1, for numerical safety-check
real(kind=real_cvprec), parameter :: almost_one = one - sqrt_min_delta


! Find points where, as it stands, the parcel core is saturated,
! but the parcel edge is subsaturated
nc = 0
do ic = 1, n_points
  ! The implied edge qcl will be negative if core_qcl > cmr * mean_qcl
  !IF ( core_q_cl(ic) > par_mean_fields(ic,i_q_cl) * core_mean_ratio(ic) ) THEN
  if ( core_q_cl(ic) / core_mean_ratio(ic) > par_mean_fields(ic,i_q_cl) ) then
    nc = nc + 1
    index_ic(nc) = ic
  end if
end do

if ( nc > 0 ) then
  ! Only anything else to do if any points found...

  if ( i_mean_q_cl == i_mean_q_cl_full ) then
    ! Parcel mean is assumed to be entirely cloudy if the parcel core
    ! contains any cloud.  This requires that the extrapolated edge
    ! value of q_cl is at least zero.

    ! Calculate amount of liquid cloud that needs to be added.
    ! This formula is the minimum q_cl the mean needs to have to
    ! avoid negative q_cl at the parcel edge:
    do ic2 = 1, nc
      ic = index_ic(ic2)
      dq_cond(ic) = core_q_cl(ic) / core_mean_ratio(ic)                        &
                  - par_mean_fields(ic,i_q_cl)
      ! Impose limit on dq_cond; must not exceed the available vapour
      dq_cond(ic) = min( dq_cond(ic), almost_one*par_mean_fields(ic,i_q_vap) )
    end do

    !ELSE ...
    ! Other options not implemented yet,
  end if  ! ( i_mean_q_cl )

  ! Add on condensation and latent heat release at these points
  call set_cp_tot( n_points, n_points_super,                                   &
                   par_mean_fields(:,i_q_vap),                                 &
                   par_mean_fields(:,i_qc_first:i_qc_last),                    &
                   cp_tot )
  call lat_heat_incr( n_points, nc, i_phase_change_con,                        &
                      cp_tot, par_mean_fields(:,i_temperature),                &
                      index_ic=index_ic, dq=dq_cond )
  do ic2 = 1, nc
    ic = index_ic(ic2)
    par_mean_fields(ic,i_q_cl)  = par_mean_fields(ic,i_q_cl)  + dq_cond(ic)
    par_mean_fields(ic,i_q_vap) = par_mean_fields(ic,i_q_vap) - dq_cond(ic)
  end do

  ! If using in-parcel cloud-fractions
  if ( l_cv_cloudfrac ) then

    ! Reset the in-parcel cloud-fractions using the modified q_cl
    call set_par_cloudfrac( n_points, n_points_super,                          &
                            par_mean_fields(:,i_q_cl),                         &
                            par_mean_fields(:,i_q_cf),                         &
                            par_mean_fields(:,i_cf_liq:i_cf_bulk) )

  end if  ! ( l_cv_cloudfrac )

end if  ! ( nc > 0 )


return
end subroutine calc_mean_q_cl_with_core

end module calc_mean_q_cl_with_core_mod
