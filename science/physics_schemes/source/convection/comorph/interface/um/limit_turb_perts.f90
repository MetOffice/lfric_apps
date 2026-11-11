! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module limit_turb_perts_mod

use um_types, only: real_umphys

implicit none

contains

! Routine applies a minimum limit to the turbulent w variance to ensure
! comorph's turbulence-based perturbations don't become stupidly large.
! This is only needed as a temporary fix to get round problems with the
! boundary-layer scheme predicting near-zero w-variance at points with
! large turbulent fluxes; since the perturbation scales with flux/sqrt(w_var),
! this leads to unrealistically large perturbations.
subroutine limit_turb_perts( z_theta, z_rho, p_layer_centres, temperature,     &
                             ftl, fqw, fu_rh, fv_rh, w_var_rh )

use atm_fields_bounds_mod, only: tdims, pdims
use nlsizes_namelist_mod, only: bl_levels
use qsat_mod, only: qsat_wat_mix

implicit none

! Model-level heights above surface
real(kind=real_umphys), intent(in) ::                                          &
                    z_theta        ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    z_rho          ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     pdims%k_start:pdims%k_end )

! Pressure at theta-levels
real(kind=real_umphys), intent(in) ::                                          &
                    p_layer_centres( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     tdims%k_start:tdims%k_end )

! Temperature
real(kind=real_umphys), intent(in) ::                                          &
                    temperature    ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )

! Turbulent fluxes on rho-levels
real(kind=real_umphys), intent(in) ::                                          &
                    ftl            ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     bl_levels )
real(kind=real_umphys), intent(in) ::                                          &
                    fqw            ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     bl_levels )
real(kind=real_umphys), intent(in) ::                                          &
                    fu_rh          ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     bl_levels )
real(kind=real_umphys), intent(in) ::                                          &
                    fv_rh          ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     bl_levels )

! Vertical velocity variance on rho-levels, to be increased where needed
real(kind=real_umphys), intent(in out) ::                                      &
                        w_var_rh   ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     bl_levels )

! Saturation vapour mixing-ratio
real(kind=real_umphys) ::                                                      &
        qsat                       ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:bl_levels )
! qsat interpolated onto current rho-level
real(kind=real_umphys) ::                                                      &
        qsat_rh                    ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )
! Required fractional increase in w_var
real :: w_var_ratio                ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     1:bl_levels )

! Interpolation weight
real(kind=real_umphys) :: interp

! Minimum allowed w standard deviation and variance
real(kind=real_umphys) :: min_w
real(kind=real_umphys) :: min_w_var

! Variables for printing info about excessive parcel perturbation check
integer :: n_inc_w_var
real(kind=real_umphys) :: mean_inc_w_var
real(kind=real_umphys) :: max_inc_w_var
integer :: istat

! Max allowed perturbation scale for water mixing-ratio
real(kind=real_umphys) :: max_qw

! Max allowed perturbation scales for temperature and winds
real(kind=real_umphys), parameter :: max_tl = 2.0
real(kind=real_umphys), parameter :: max_uv = 5.0

! Switch for whether to print info about points where limit on w_var imposed
logical, parameter :: l_print_turb_perts_exceeded = .false.

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE( i, j, k, interp, qsat_rh, max_qw, min_w, min_w_var )            &
!$OMP SHARED( bl_levels, tdims, pdims, qsat, temperature, p_layer_centres,     &
!$OMP         z_rho, z_theta, ftl, fqw, fu_rh, fv_rh, w_var_ratio, w_var_rh )

! Calculate qsat (max allowed perturbation to q scales with qsat).
! Note that making max allowed q pert scale with q runs into problems
! when q goes negative in some runs!  Using qsat to avoid this...
!$OMP DO SCHEDULE(STATIC)
do k = 1, bl_levels
  call qsat_wat_mix( qsat(:,:,k), temperature(:,:,k), p_layer_centres(:,:,k),  &
                     tdims%i_len, tdims%j_len )
end do
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
do k = 1, bl_levels
  if ( k==1 ) then
    ! Bottom rho-level; just use qsat from 1st theta-level
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        qsat_rh(i,j) = qsat(i,j,k)
      end do
    end do
  else
    ! Other rho-levels; interpolate between neighbouring theta-levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        interp = ( z_rho(i,j,k) - z_theta(i,j,k-1) )                           &
               / ( z_theta(i,j,k) - z_theta(i,j,k-1) )
        qsat_rh(i,j) = (1.0-interp) * qsat(i,j,k-1) + interp * qsat(i,j,k)
      end do
    end do
  end if
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ! Max allowed q perturbation scales with qsat
      max_qw = 0.5 * qsat_rh(i,j)
      ! Find min turbulent w std that ensures all perturbations
      ! will be below their specified upper bounds:
      min_w = maxval([ abs(ftl(i,j,k)) / max_tl,                               &
                       abs(fqw(i,j,k)) / max_qw,                               &
                       abs(fu_rh(i,j,k)) / max_uv,                             &
                       abs(fv_rh(i,j,k)) / max_uv ])
      ! Convert w std to variance
      min_w_var = min_w * min_w
      ! Calculate ratio of min w-variance over existing w-variance
      w_var_ratio(i,j,k) = min_w_var / w_var_rh(i,j,k)
      ! Impose limit on w_var_rh
      w_var_rh(i,j,k) = max( w_var_rh(i,j,k), min_w_var )
    end do
  end do
end do  ! k = 1, bl_levels
!$OMP END DO NOWAIT

!$OMP END PARALLEL

return
end subroutine limit_turb_perts

end module limit_turb_perts_mod
