! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Area cloud parameterisation for use with PC2 Cloud Scheme.

module pc2_arcld_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_ARCLD_MOD'
contains

subroutine pc2_arcld(                                                          &
!      Pressure related fields
 p_layer_centres, p_layer_boundaries, cumulus, rhcrit,                         &
!      Array dimensions
 rhc_row_length, rhc_rows,zlcl_mixed,                                          &
 large_levels, levels_per_level,                                               &
!      Prognostic Fields
 cf_area, t, cf, cfl, cff, q, qcl, qcf, rhts,                                  &
 tlts, qtts, ptts,                                                             &
!      Logical control
 l_mixing_ratio)

use planet_constants_mod,  only: lcrcp
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims, pdims_l
use level_heights_mod,     only: r_theta_levels

use qsat_mod, only: qsat_wat, qsat_wat_mix

use pc2_initiate_mod, only: pc2_initiate
implicit none

! Description:
!   Cusack-like vertical interpolation onto 3 sub-levels to calculate
!   cloud fraction and condensate in PC2 initiation. Also area cloud
!   parameterisation for use with radiation scheme.
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

!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
integer ::                                                                     &
                      !, intent(in)
 rhc_row_length, rhc_rows,                                                     &
!       Dimensions of the rhcrit variable.
   large_levels,                                                               &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((levels - 2)*levels_per_level) + 2
   levels_per_level
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops

logical ::                                                                     &
                      !, intent(in)
 cumulus(           tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end),                                &
!       Is this a boundary layer cumulus point
   l_mixing_ratio
!       Use mixing ratio formulation

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
real(kind=real_umphys) :: zlcl_mixed(                                          &
                   pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 p_layer_centres(   pdims%i_start:pdims%i_end,                                 &
                    pdims%j_start:pdims%j_end,                                 &
                    0:pdims%k_end),                                            &
!       Pressure at all points, on theta levels (Pa).
   p_layer_boundaries(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end,                               &
                      0:pdims%k_end),                                          &
!       Pressure at all points, on u,v levels (Pa).
   rhcrit(            rhc_row_length,                                          &
                      rhc_rows,                                                &
                                  1:tdims%k_end),                              &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
   cff(               tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Ice cloud fraction (no units)
   qcf(               tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end)
!       Cloud ice content at processed levels (kg water per kg air).

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 t(                 tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end,                                 &
                                1:tdims%k_end),                                &
!       Temperature (K)
   cf(                tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Total cloud fraction (no units)
   cfl(               tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Liquid cloud fraction (no units)
   q(                 tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Vapour content (kg water per kg air)
   qcl(               tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Liquid content (kg water per kg air)
   rhts(              tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Variable carrying initial RHT wrt TL from start of timestep
   tlts(              tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       TL at start of timestep
   qtts(              tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       qT at start of timestep
   ptts(              pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end,                               &
                      pdims%k_start:pdims%k_end),                              &
!       Pressure at theta levels at start of timestep
   cf_area(           tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end)
!       Area cloud fraction

! --------------------------------------------------------------------
! Local variables
! ---------------------------------------------------------------------
integer :: i,j,k   ! Loop counters: k   - vertical level index.
!                                     i,j - horizontal field index.
integer :: k_index ! Extra loop counter for large arrays.

real(kind=real_umphys) ::                                                      &
  inverse_level,                                                               &
                   ! Set to (1. / levels_per_level)
  qt_norm_next,                                                                &
                   ! Temporary space for qT_norm
  stretcher,                                                                   &
  delta_p          ! Layer pressure thickness * inverse_level

real(kind=real_umphys) ::                                                      &
  qsl(              tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end),                                &
!       Saturated specific humidity for temp TL or T.
    qcl_latest(       tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Cloud liquid content at processed levels (kg water per kg air).
    tl(               tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Liquid temperature (TL) (K).
    qt_norm(          tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end),                              &
!       Total water content normalized to qSAT_WAT.
    rhcrit_large(rhc_row_length,rhc_rows,large_levels),                        &
!
!       Values of quantities on large_levels
    p_large(          pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end,                               &
                      large_levels),                                           &
    r_large(          pdims_l%i_start:pdims_l%i_end,                           &
                      pdims_l%j_start:pdims_l%j_end,                           &
                      0:large_levels),                                         &
!
    t_large(          tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    q_large(          tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    qcl_large(        tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    cf_large(         tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    cfl_large(        tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    cff_large(        tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    rhts_large(       tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    tlts_large(       tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    qtts_large(       tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      large_levels),                                           &
!
    ptts_large(       pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end,                               &
                      large_levels)

logical ::                                                                     &
 linked(            tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end,                                 &
                                1:tdims%k_end)
!       True for sub-layers that have similar supersaturation properties

!  Local parameters and other physical constants------------------------
real(kind=real_umphys), parameter :: drat_thresh =3.0e-1
!       Test for continuity of sub-levels
real(kind=real_umphys), parameter :: tol_test  =1.0e-11
!       Tolerance for non-zero humidities

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_ARCLD'

! ---------------------------------------------------------------------
! Code starts
! ---------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
inverse_level = 1.0 / levels_per_level

! Create new arrays for TL and current qcl

do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      tl(i,j,k)         = t(i,j,k) - lcrcp*qcl(i,j,k)
      qcl_latest(i,j,k) = qcl(i,j,k)
    end do !i
  end do !j
end do !k

! Test for continuity between adjacent layers based on supersaturation
! (qt - qsl) / qsl : as we take differences the - qsl term drops out.

if ( l_mixing_ratio ) then
  call qsat_wat_mix( qsl, tl(:,:,1), p_layer_centres(:,:,1),                   &
          tdims%i_len,tdims%j_len)
else
  call qsat_wat( qsl, tl(:,:,1), p_layer_centres(:,:,1),                       &
          tdims%i_len,tdims%j_len)
end if

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    qt_norm(i,j) =(q(i,j,1)+qcl_latest(i,j,1)+qcf(i,j,1))/qsl(i,j)
  end do
end do

if ( l_mixing_ratio ) then
  call qsat_wat_mix( qsl, tl(:,:,2), p_layer_centres(:,:,2),                   &
          tdims%i_len,tdims%j_len)
else
  call qsat_wat( qsl, tl(:,:,2), p_layer_centres(:,:,2),                       &
          tdims%i_len,tdims%j_len)
end if

! Do nothing to top and bottom layers
do j = 1, rhc_rows
  do i = 1, rhc_row_length
    rhcrit_large(i,j,1)            = rhcrit(i,j,1)
    rhcrit_large(i,j,large_levels) = rhcrit(i,j,tdims%k_end)
  end do
end do

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end

    r_large   (i,j,0) = r_theta_levels(i,j,0)

    p_large   (i,j,1) = p_layer_centres(i,j,1)
    r_large   (i,j,1) = r_theta_levels(i,j,1)
    t_large   (i,j,1) = t(i,j,1)
    q_large   (i,j,1) = q(i,j,1)
    qcl_large (i,j,1) = qcl_latest(i,j,1)
    cf_large  (i,j,1) = cf(i,j,1)
    cfl_large (i,j,1) = cfl(i,j,1)
    cff_large (i,j,1) = cff(i,j,1)
    tlts_large(i,j,1) = tlts(i,j,1)
    qtts_large(i,j,1) = qtts(i,j,1)
    ptts_large(i,j,1) = ptts(i,j,1)
    rhts_large(i,j,1) = rhts(i,j,1)

    p_large   (i,j,large_levels) = p_layer_centres(i,j,pdims%k_end)
    r_large   (i,j,large_levels) = r_theta_levels(i,j,pdims%k_end)
    t_large   (i,j,large_levels) = t(i,j,tdims%k_end)
    q_large   (i,j,large_levels) = q(i,j,tdims%k_end)
    qcl_large (i,j,large_levels) = qcl_latest(i,j,tdims%k_end)
    cf_large  (i,j,large_levels) = cf(i,j,tdims%k_end)
    cfl_large (i,j,large_levels) = cfl(i,j,tdims%k_end)
    cff_large (i,j,large_levels) = cff(i,j,tdims%k_end)
    tlts_large(i,j,large_levels) = tlts(i,j,tdims%k_end)
    qtts_large(i,j,large_levels) = qtts(i,j,tdims%k_end)
    ptts_large(i,j,large_levels) = ptts(i,j,tdims%k_end)
    rhts_large(i,j,large_levels) = rhts(i,j,tdims%k_end)

    ! Test for continuity (assumed if linked is .true.)
    qt_norm_next=(q(i,j,2)+qcl_latest(i,j,2)+qcf(i,j,2))/qsl(i,j)
    linked(i,j,1) =                                                            &
          (drat_thresh >= abs(qt_norm(i,j) - qt_norm_next))
    qt_norm(i,j) = qt_norm_next

  end do !i
end do ! j

do k = 2, (tdims%k_end - 1)
  k_index = 3 + (levels_per_level * (k-2))

  if ( l_mixing_ratio ) then
    call qsat_wat_mix( qsl, tl(:,:,k+1), p_layer_centres(:,:,k+1),             &
            tdims%i_len,tdims%j_len)
  else
    call qsat_wat( qsl, tl(:,:,k+1), p_layer_centres(:,:,k+1),                 &
            tdims%i_len,tdims%j_len)
  end if

  ! Select associated rhcrit values
  do j = 1, rhc_rows
    do i = 1, rhc_row_length
      rhcrit_large(i,j,k_index-1) = rhcrit(i,j,k)
      rhcrit_large(i,j,k_index)   = rhcrit(i,j,k)
      rhcrit_large(i,j,k_index+1) = rhcrit(i,j,k)
      ! don't interpolate r_theta_levels
    end do
  end do
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      r_large(i,j,k_index-1) = r_theta_levels(i,j,k)
      r_large(i,j,k_index)   = r_theta_levels(i,j,k)
      r_large(i,j,k_index+1) = r_theta_levels(i,j,k)
    end do
  end do

  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Test for continuity (assumed if linked = .true.)
      qt_norm_next  = (q(i,j,(k+1)) + qcl_latest(i,j,(k+1))                    &
                    + qcf(i,j,(k+1))) / qsl(i,j)
      linked(i,j,k) =                                                          &
          (drat_thresh >= abs(qt_norm(i,j) - qt_norm_next))
      qt_norm(i,j)  = qt_norm_next
      !
      ! Select interpolated pressure levels
      delta_p = (p_layer_boundaries(i,j,(k-1)) -                               &
                 p_layer_boundaries(i,j,k))    * inverse_level
      if (p_layer_centres(i,j,k) >=                                            &
         (p_layer_boundaries(i,j,k) + delta_p)) then
        p_large(i,j,k_index) = p_layer_centres(i,j,k)
      else
        p_large(i,j,k_index) = 0.5*(p_layer_boundaries(i,j,k) +                &
                               p_layer_boundaries(i,j,(k-1)))
      end if
      p_large(i,j,(k_index-1)) = p_large(i,j,k_index)+delta_p
      p_large(i,j,(k_index+1)) = p_large(i,j,k_index)-delta_p
      !
      ! Select variable values at layer centres
      t_large   (i,j,k_index) = t(i,j,k)
      q_large   (i,j,k_index) = q(i,j,k)
      qcl_large (i,j,k_index) = qcl_latest(i,j,k)
      cf_large  (i,j,k_index) = cf(i,j,k)
      cfl_large (i,j,k_index) = cfl(i,j,k)
      cff_large (i,j,k_index) = cff(i,j,k)
      tlts_large(i,j,k_index) = tlts(i,j,k)
      qtts_large(i,j,k_index) = qtts(i,j,k)
      ptts_large(i,j,k_index) = ptts(i,j,k)
      rhts_large(i,j,k_index) = rhts(i,j,k)

      ! Calculate increment in variable values, pressure interpolation
      ! NB: Using X_large(i,j,(k_index+1)) as store for X increments
      if ( linked(i,j,(k-1)) ) then
        if ( linked(i,j,k) ) then
          !               Interpolate from level k-1 to k+1
          stretcher = delta_p /                                                &
         (p_layer_centres(i,j,k-1)-p_layer_centres(i,j,k+1))
          t_large(i,j,(k_index+1)) = stretcher *                               &
         (t(i,j,(k+1)) - t(i,j,(k-1)))
          q_large(i,j,(k_index+1)) = stretcher *                               &
         (q(i,j,(k+1)) - q(i,j,(k-1)))
          qcl_large(i,j,(k_index+1)) = stretcher *                             &
         (qcl_latest(i,j,(k+1)) - qcl_latest(i,j,(k-1)))
          tlts_large(i,j,(k_index+1)) = stretcher *                            &
         (tlts(i,j,(k+1)) - tlts(i,j,(k-1)))
          qtts_large(i,j,(k_index+1)) = stretcher *                            &
         (qtts(i,j,(k+1)) - qtts(i,j,(k-1)))
          ptts_large(i,j,(k_index+1)) = stretcher *                            &
         (ptts(i,j,(k+1)) - ptts(i,j,(k-1)))
        else
          !               Interpolate from level k-1 to k
          stretcher = delta_p /                                                &
         (p_layer_centres(i,j,k-1) - p_large(i,j,k_index))
          t_large(i,j,(k_index+1)) = stretcher *                               &
         (t_large(i,j,k_index) - t(i,j,(k-1)))
          q_large(i,j,(k_index+1)) = stretcher *                               &
         (q_large(i,j,k_index) - q(i,j,(k-1)))
          qcl_large(i,j,(k_index+1)) = stretcher *                             &
         (qcl_large(i,j,k_index) - qcl_latest(i,j,(k-1)))
          tlts_large(i,j,(k_index+1)) = stretcher *                            &
         (tlts_large(i,j,k_index) - tlts(i,j,(k-1)))
          qtts_large(i,j,(k_index+1)) = stretcher *                            &
         (qtts_large(i,j,k_index) - qtts(i,j,(k-1)))
          ptts_large(i,j,(k_index+1)) = stretcher *                            &
         (ptts_large(i,j,k_index) - ptts(i,j,(k-1)))
        end if

      else
        if ( linked(i,j,k) ) then
          !               Interpolate from level k to k+1
          stretcher = delta_p /                                                &
          (p_large(i,j,k_index) - p_layer_centres(i,j,k+1))
          t_large(i,j,(k_index+1)) = stretcher *                               &
         (t(i,j,(k+1)) - t_large(i,j,k_index))
          q_large(i,j,(k_index+1)) = stretcher *                               &
         (q(i,j,(k+1)) - q_large(i,j,k_index))
          qcl_large(i,j,(k_index+1)) = stretcher *                             &
         (qcl_latest(i,j,(k+1)) - qcl_large(i,j,k_index))
          tlts_large(i,j,(k_index+1)) = stretcher *                            &
         (tlts(i,j,(k+1)) - tlts_large(i,j,k_index))
          qtts_large(i,j,(k_index+1)) = stretcher *                            &
         (qtts(i,j,(k+1)) - qtts_large(i,j,k_index))
          ptts_large(i,j,(k_index+1)) = stretcher *                            &
         (ptts(i,j,(k+1)) - ptts_large(i,j,k_index))
        else
          !               No interpolation, freeze at level k
          t_large   (i,j,(k_index+1)) = 0.0
          q_large   (i,j,(k_index+1)) = 0.0
          qcl_large (i,j,(k_index+1)) = 0.0
          tlts_large(i,j,(k_index+1)) = 0.0
          qtts_large(i,j,(k_index+1)) = 0.0
          ptts_large(i,j,(k_index+1)) = 0.0
        end if

      end if

      ! Protect against q or qcl going negative (T would imply blow-up anyway)
      if (q_large(i,j,k_index)  <                                              &
                 (abs(q_large(i,j,(k_index+1)))+tol_test))                     &
                      q_large(i,j,(k_index+1)) = 0.0
      if (qcl_large(i,j,k_index)  <                                            &
                   (abs(qcl_large(i,j,(k_index+1)))+tol_test))                 &
                        qcl_large(i,j,(k_index+1)) = 0.0
      if (qtts_large(i,j,k_index) <                                            &
                  (abs(qtts_large(i,j,(k_index+1)))+tol_test))                 &
                       qtts_large(i,j,(k_index+1)) = 0.0

      ! Select variable values at level below layer centre
      t_large   (i,j,(k_index-1)) = t_large(i,j,k_index) -                     &
                                      t_large(i,j,(k_index+1))
      q_large   (i,j,(k_index-1)) = q_large(i,j,k_index) -                     &
                                      q_large(i,j,(k_index+1))
      qcl_large (i,j,(k_index-1)) = qcl_large(i,j,k_index)-                    &
                                      qcl_large(i,j,(k_index+1))
      cf_large  (i,j,(k_index-1)) = cf(i,j,k)
      cfl_large (i,j,(k_index-1)) = cfl(i,j,k)
      cff_large (i,j,(k_index-1)) = cff(i,j,k)
      tlts_large(i,j,(k_index-1)) = tlts_large(i,j,k_index) -                  &
                                      tlts_large(i,j,(k_index+1))
      qtts_large(i,j,(k_index-1)) = qtts_large(i,j,k_index) -                  &
                                      qtts_large(i,j,(k_index+1))
      ptts_large(i,j,(k_index-1)) = ptts_large(i,j,k_index) -                  &
                                      ptts_large(i,j,(k_index+1))

      ! Select variable values at level above layer centre
      ! NB: CEASE using X_large(i,j,(k_index+1)) as store for X increments.
      t_large   (i,j,(k_index+1)) = t_large(i,j,(k_index+1)) +                 &
                                      t_large(i,j,k_index)
      q_large   (i,j,(k_index+1)) = q_large(i,j,(k_index+1)) +                 &
                                      q_large(i,j,k_index)
      qcl_large (i,j,(k_index+1)) = qcl_large(i,j,(k_index+1)) +               &
                                      qcl_large(i,j,k_index)
      cf_large  (i,j,(k_index+1)) = cf(i,j,k)
      cfl_large (i,j,(k_index+1)) = cfl(i,j,k)
      cff_large (i,j,(k_index+1)) = cff(i,j,k)
      tlts_large(i,j,(k_index+1)) = tlts_large(i,j,(k_index+1)) +              &
                                      tlts_large(i,j,k_index)
      qtts_large(i,j,(k_index+1)) = qtts_large(i,j,(k_index+1)) +              &
                                      qtts_large(i,j,k_index)
      ptts_large(i,j,(k_index+1)) = ptts_large(i,j,(k_index+1)) +              &
                                      ptts_large(i,j,k_index)

    end do
  end do

  ! Calculate RH above and below layer centres
  if ( l_mixing_ratio ) then
    call qsat_wat_mix(qsl, tlts_large(:,:,k_index-1),                          &
                          ptts_large(:,:,k_index-1),tdims%i_len,tdims%j_len)
  else
    call qsat_wat(qsl, tlts_large(:,:,k_index-1),                              &
                          ptts_large(:,:,k_index-1),tdims%i_len,tdims%j_len)
  end if

  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      rhts_large(i,j,k_index-1)=qtts_large(i,j,k_index-1)                      &
                                  /qsl(i,j)
    end do
  end do

  if ( l_mixing_ratio ) then
    call qsat_wat_mix(qsl, tlts_large(:,:,k_index+1),                          &
                          ptts_large(:,:,k_index+1),tdims%i_len,tdims%j_len)
  else
    call qsat_wat(qsl, tlts_large(:,:,k_index+1),                              &
                          ptts_large(:,:,k_index+1),tdims%i_len,tdims%j_len)
  end if

  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      rhts_large(i,j,k_index+1)=qtts_large(i,j,k_index+1)                      &
                                  /qsl(i,j)
    end do
  end do

end do

call pc2_initiate(p_large,cumulus,rhcrit_large,                                &
  large_levels, rhc_row_length,rhc_rows,zlcl_mixed,r_large,                    &
  t_large,cf_large,cfl_large,cff_large,q_large,qcl_large,                      &
  rhts_large,l_mixing_ratio)

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    t         (i,j,1)           = t_large(i,j,1)
    q         (i,j,1)           = q_large(i,j,1)

    cf_area   (i,j,1)           = cf_large(i,j,1)
    cf        (i,j,1)           = cf_large(i,j,1)

    qcl_latest(i,j,1)           = qcl_large(i,j,1)
    cfl       (i,j,1)           = cfl_large(i,j,1)

    t         (i,j,tdims%k_end) = t_large(i,j,large_levels)
    q         (i,j,tdims%k_end) = q_large(i,j,large_levels)

    cf_area   (i,j,tdims%k_end) = cf_large(i,j,large_levels)
    cf        (i,j,tdims%k_end) = cf_large(i,j,large_levels)

    qcl_latest(i,j,tdims%k_end) = qcl_large(i,j,large_levels)
    cfl       (i,j,tdims%k_end) = cfl_large(i,j,large_levels)
  end do
end do

! Output variables for remaining layers
do k = 2, (tdims%k_end - 1)

  k_index = 3 + (levels_per_level * (k-2))

  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Area cloud fraction is maximum of sub-layer cloud fractions
      cf_area(i,j,k) =                                                         &
            max(cf_large(i,j,k_index),                                         &
           (max(cf_large(i,j,(k_index+1)),                                     &
                cf_large(i,j,(k_index-1)))) )

      ! Bulk cloud fraction is mean of sub-layer cloud fractions : strictly
      ! this is a pressure weighted mean being used to approximate a volume
      ! mean. Over the depth of a layer the difference should not be large.
      cf(i,j,k) = inverse_level *                                              &
         ( cf_large(i,j,(k_index-1)) +                                         &
           cf_large(i,j, k_index)    +                                         &
           cf_large(i,j,(k_index+1)) )

      ! The pressure weighted mean of qcf is the input qcf: do not update.

      ! Qcl is the pressure weighted mean of qcl from each sub-layer.
      qcl_latest(i,j,k) = inverse_level *                                      &
    ( qcl_large(i,j,(k_index-1)) +                                             &
      qcl_large(i,j,k_index) + qcl_large(i,j,(k_index+1)) )

      ! Liq. cloud fraction is mean of sub-layer cloud fractions : strictly
      ! this is a pressure weighted mean being used to approximate a volume
      ! mean. Over the depth of a layer the difference should not be large.
      cfl(i,j,k) = inverse_level *                                             &
    ( cfl_large(i,j,(k_index-1)) +                                             &
      cfl_large(i,j,k_index)     +                                             &
      cfl_large(i,j,(k_index+1)) )

      ! Update Q
      ! Update T
      ! Move qcl_latest into qcl.
      q(i,j,k)   = q(i,j,k) + qcl(i,j,k) - qcl_latest(i,j,k)
      t(i,j,k)   = t(i,j,k) - (qcl(i,j,k)*lcrcp) +                             &
                              (qcl_latest(i,j,k) * lcrcp)
      qcl(i,j,k) = qcl_latest(i,j,k)

    end do
  end do

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_arcld
! ======================================================================
end module pc2_arcld_mod
