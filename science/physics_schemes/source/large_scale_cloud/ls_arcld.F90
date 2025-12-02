! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Large-scale Area (Vertical Gradient) Cloud Scheme.
module ls_arcld_mod

use ls_cld_mod, only: ls_cld

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'LS_ARCLD_MOD'
contains

subroutine ls_arcld(                                                           &
!      Pressure related fields
  p_layer_centres, rhcrit, p_layer_boundaries,                                 &
!      Array dimensions
  rhc_row_length, rhc_rows, bl_levels,                                         &
  levels_per_level, large_levels,                                              &
!      Needed for LS_ACF_Brooks
  fv_cos_theta_latitude,                                                       &
!      Convection diagnosis information (only used for A05_4A)
  ntml, cumulus, l_mixing_ratio,                                               &
!      Prognostic Fields
  qcf_latest,                                                                  &
  t_latest, q_latest, qcl_latest,                                              &
!      Various cloud fractions
  area_cloud_fraction, bulk_cloud_fraction,                                    &
  cloud_fraction_liquid, cloud_fraction_frozen,                                &
  error_code)

use planet_constants_mod, only: lcrcp
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim
use atm_fields_bounds_mod,only: tdims, pdims, tdims_s
use cloud_inputs_mod,     only: i_cld_area
use pc2_constants_mod,    only: acf_off, acf_cusack, acf_brooks

use qsat_mod, only: qsat_wat, qsat_wat_mix

use ls_acf_brooks_mod, only: ls_acf_brooks

implicit none

! Purpose:
!   This subroutine calculates liquid and ice cloud fractional cover
!   for use with the enhanced precipitation microphysics scheme. It
!   also returns area and bulk cloud fractions for use in radiation.

! Method:
!   Statistical cloud scheme separates input moisture into specific
!   humidity and cloud liquid water. Temperature calculated from liquid
!   water temperature. Cloud fractions calculated from statistical
!   relation between cloud fraction and cloud liquid/ice water content.
!   Area cloud fraction calculated by subdividing layers, calculating
!   cloud on each sublayer and taking maximum (use mean for bulk cloud).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.29

!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------

! arguments with intent in. ie: input variables.

integer, intent(in) ::                                                         &
  bl_levels,                                                                   &
!       No. of boundary layer levels
    rhc_row_length, rhc_rows,                                                  &
!       Horizontal dimensions of rh_crit diagnostic.
    levels_per_level,                                                          &
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops
    large_levels,                                                              &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((tdims%k_end - 2)*levels_per_level) + 2
    ntml(               pdims%i_start:pdims%i_end,                             &
                        pdims%j_start:pdims%j_end)
!       Height of diagnosed BL top

logical, intent(in) ::                                                         &
 l_mixing_ratio
!       True if using mixing rations

logical, intent(in) ::                                                         &
  cumulus(            tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end)
!       Logical: indicator of convection

real(kind=real_umphys), intent(in) ::                                          &
 qcf_latest(          tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                                  1:tdims%k_end),                              &
!       Cloud ice content at processed levels (kg water per kg air).
   p_layer_centres(      pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,                            &
                                     0:pdims%k_end),                           &
!       Pressure at all points, on theta levels (Pa).
   rhcrit(               rhc_row_length,                                       &
                         rhc_rows,                                             &
                         1:tdims%k_end),                                       &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for a given
!       set of levels.
   p_layer_boundaries(   pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,                            &
                                     0:pdims%k_end),                           &
!       Pressure at all points, on u,v levels (Pa).
   fv_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                        &
                         tdims_s%j_start:tdims_s%j_end)
!       Finite volume cos(lat)

! arguments with intent in/out. ie: input variables changed on output.

real(kind=real_umphys), intent(in out) ::                                      &
 q_latest(             tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end),                             &
!       On input : Vapour + liquid water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
   t_latest(             tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                                     1:tdims%k_end)
!       On input : Liquid water temperature (TL) (K).
!       On output: Temperature at processed levels (K).

! arguments with intent out. ie: output variables.

!     Error Status:
integer, intent(out) :: error_code  ! 0 if OK; 1 if bad arguments.

real(kind=real_umphys), intent(out) ::                                         &
 qcl_latest(           tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end),                             &
!       Cloud liquid content at processed levels (kg water per kg air).
   area_cloud_fraction(  tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                                     1:tdims%k_end),                           &
!       Area cloud fraction at processed levels (decimal fraction).
   bulk_cloud_fraction(  tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                                     1:tdims%k_end),                           &
!       Cloud fraction at processed levels (decimal fraction).
   cloud_fraction_liquid(tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                                     1:tdims%k_end),                           &
!       Liquid cloud fraction at processed levels (decimal fraction).
   cloud_fraction_frozen(tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                                     1:tdims%k_end)
!       Frozen cloud fraction at processed levels (decimal fraction).

!  Local parameters and other physical constants------------------------
real(kind=real_umphys) :: drat_thresh ! Test for continuity of sub-levels
real(kind=real_umphys) :: tol_test ! Tolerance for non-zero humidities
parameter (drat_thresh=3.0e-1, tol_test=1.0e-11)

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
real(kind=real_umphys) ::                                                      &
  qt_norm_next,                                                                &
                     ! Temporary space for qT_norm
  stretcher,                                                                   &

  inverse_level,                                                               &
                     ! Set to (1. / levels_per_level)
  delta_p        ! Layer pressure thickness * inverse_level

!  (b) Others.
integer :: i,j,k      ! Loop counters: k - vertical level index.
!                                       i,j - horizontal field index.
integer :: k_index    ! Extra loop counter for large arrays.

!  Local dynamic arrays-------------------------------------------------
!    11 blocks of real workspace are required.
real(kind=real_umphys) ::                                                      &
  qsl(                 tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end),                             &
!        Saturated specific humidity for temp TL or T.
    qt_norm(             tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end),                           &
!        Total water content normalized to qSAT_WAT.
    rhcrit_large(rhc_row_length,rhc_rows,large_levels),                        &

    p_large(             tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end, large_levels),             &
    t_large(             tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end, large_levels),             &
    q_large(             tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end, large_levels),             &
    qcf_large(           tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end, large_levels),             &
    cloud_fraction_large(tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end, large_levels),             &
    qcl_large(           tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end, large_levels),             &
    cloud_fraction_liquid_large(                                               &
                         tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end, large_levels),             &
    cloud_fraction_frozen_large(                                               &
                         tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end, large_levels)

integer ::                                                                     &
  ntml_large(          tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end)

logical ::                                                                     &
 linked(               tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end)
!       True for sub-layers that have similar supersaturation properties

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LS_ARCLD'


!- End of Header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
! ----------------------------------------------------------------------
error_code=0

! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------

if (i_cld_area == acf_off) then
  !     As before
  call ls_cld( p_layer_centres(1,1,1), rhcrit,                                 &
               tdims%k_end, bl_levels,                                         &
               rhc_row_length, rhc_rows,                                       &
               ntml,cumulus,                                                   &
               l_mixing_ratio,t_latest, bulk_cloud_fraction,                   &
               q_latest, qcf_latest, qcl_latest,                               &
               cloud_fraction_liquid,                                          &
               cloud_fraction_frozen, error_code )

  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        area_cloud_fraction(i,j,k)=bulk_cloud_fraction(i,j,k)
      end do
    end do
  end do

else if (i_cld_area == acf_cusack) then
  !       Vertical gradient area cloud option

  inverse_level = 1.0 / levels_per_level

  ! Test for continuity between adjacent layers based on supersaturation
  ! (qt - qsl) / qsl : as we take differences the - qsl term drops out.
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED( rhc_rows, rhc_row_length, rhcrit_large, rhcrit, tdims,           &
!$OMP         large_levels, p_large, t_large, q_large, qcf_large,              &
!$OMP         linked, qt_norm, p_layer_centres, t_latest, q_latest,            &
!$OMP         qcf_latest, qsl, l_mixing_ratio)                                 &
!$OMP private(i, j, qt_norm_next )
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end

    if ( l_mixing_ratio ) then
      call qsat_wat_mix(qsl(:,j),t_latest(:,j,1),                              &
                            p_layer_centres(:,j,1), tdims%i_len)
    else
      call qsat_wat(qsl(:,j),t_latest(:,j,1),                                  &
                            p_layer_centres(:,j,1), tdims%i_len)
    end if

    do i = tdims%i_start, tdims%i_end
      qt_norm(i,j) =(q_latest(i,j,1)+qcf_latest(i,j,1))/qsl(i,j)
    end do

    if ( l_mixing_ratio ) then
      call qsat_wat_mix(qsl(:,j),t_latest(:,j,2),                              &
                            p_layer_centres(:,j,2), tdims%i_len)
    else
      call qsat_wat(qsl(:,j),t_latest(:,j,1),                                  &
                            p_layer_centres(:,j,1), tdims%i_len)
    end if

  end do
!$OMP end do NOWAIT

  ! Do nothing to top and bottom layers
!$OMP do SCHEDULE(STATIC)
  do j = 1, rhc_rows
    do i = 1, rhc_row_length
      rhcrit_large(i,j,1) = rhcrit(i,j,1)
      rhcrit_large(i,j,large_levels)= rhcrit(i,j,tdims%k_end)
    end do
  end do
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      p_large(i,j,1) = p_layer_centres(i,j,1)
      t_large(i,j,1) = t_latest(i,j,1)
      q_large(i,j,1) = q_latest(i,j,1)
      qcf_large(i,j,1) = qcf_latest(i,j,1)

      p_large(i,j,large_levels) =                                              &
        p_layer_centres(i,j,tdims%k_end)

      t_large(i,j,large_levels) = t_latest(i,j,tdims%k_end)
      q_large(i,j,large_levels) = q_latest(i,j,tdims%k_end)
      qcf_large(i,j,large_levels) = qcf_latest(i,j,tdims%k_end)
      ! Test for continuity (assumed if linked is .true.)
      qt_norm_next=(q_latest(i,j,2)+qcf_latest(i,j,2))/qsl(i,j)
      linked(i,j,1) =                                                          &
            (drat_thresh  >=  abs(qt_norm(i,j) - qt_norm_next))
      qt_norm(i,j) = qt_norm_next
    end do
  end do
!$OMP end do
!$OMP end PARALLEL

  !Parameters used: drat_thresh
!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED(levels_per_level, t_latest,                                      &
!$OMP  p_layer_centres, l_mixing_ratio, rhc_rows,                              &
!$OMP  rhc_row_length, rhcrit_large, rhcrit, inverse_level,                    &
!$OMP  p_layer_boundaries, p_large, t_large, q_large, q_latest,                &
!$OMP  qcf_large, qcf_latest, linked, qt_norm_next, qt_norm, tdims,            &
!$OMP  ntml_large, ntml )                                                      &
!$OMP  private(i, j, k, k_index, stretcher, delta_p, qsl)

!$OMP do SCHEDULE(STATIC,1) ORDERED
  do k = 2, (tdims%k_end - 1)

    if ( l_mixing_ratio ) then
      call qsat_wat_mix(qsl, t_latest(:,:,(k+1)),                              &
                            p_layer_centres(:,:,(k+1)),                        &
                            tdims%i_len,tdims%j_len)
    else
      call qsat_wat(qsl, t_latest(:,:,(k+1)),                                  &
                            p_layer_centres(:,:,(k+1)),                        &
                            tdims%i_len,tdims%j_len)
    end if

!$OMP ORDERED
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        ! Test for continuity (assumed if linked = .true.)
        qt_norm_next=(q_latest(i,j,(k+1))+qcf_latest(i,j,(k+1)))               &
                     / qsl(i,j)
        linked(i,j,k) =                                                        &
            (drat_thresh  >=  abs(qt_norm(i,j) - qt_norm_next))
        qt_norm(i,j) = qt_norm_next
      end do
    end do
!$OMP end ORDERED

  end do
!$OMP end do
  !Implicit barrier

!$OMP do SCHEDULE(STATIC,1)
  do k = 2, (tdims%k_end - 1)
    k_index = 3 + (levels_per_level * (k-2))

    ! Select associated rhcrit values
    do j = 1, rhc_rows
      do i = 1, rhc_row_length
        rhcrit_large(i,j,k_index-1) = rhcrit(i,j,k)
        rhcrit_large(i,j,k_index)   = rhcrit(i,j,k)
        rhcrit_large(i,j,k_index+1) = rhcrit(i,j,k)
      end do
    end do

    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

        ! Select interpolated pressure levels
        delta_p = (p_layer_boundaries(i,j,(k-1)) -                             &
                   p_layer_boundaries(i,j,k))    * inverse_level
        if (p_layer_centres(i,j,k)  >=                                         &
           (p_layer_boundaries(i,j,k) + delta_p)) then
          p_large(i,j,k_index) = p_layer_centres(i,j,k)
        else
          p_large(i,j,k_index) = 0.5*(p_layer_boundaries(i,j,k) +              &
                                 p_layer_boundaries(i,j,(k-1)))
        end if
        p_large(i,j,(k_index-1)) = p_large(i,j,k_index)+delta_p
        p_large(i,j,(k_index+1)) = p_large(i,j,k_index)-delta_p

        ! Select variable values at layer centres
        t_large(i,j,k_index) = t_latest(i,j,k)
        q_large(i,j,k_index) = q_latest(i,j,k)
        qcf_large(i,j,k_index) = qcf_latest(i,j,k)

        ! Calculate increment in variable values, pressure interpolation
        ! NB: Using X_large(i,j,(k_index+1)) as store for X increments

        if ( linked(i,j,(k-1)) ) then
          if ( linked(i,j,k) ) then
            !               Interpolate from level k-1 to k+1
            stretcher = delta_p /                                              &
           (p_layer_centres(i,j,k-1)-p_layer_centres(i,j,k+1))
            t_large(i,j,(k_index+1)) = stretcher *                             &
           (t_latest(i,j,(k+1)) - t_latest(i,j,(k-1)))
            q_large(i,j,(k_index+1)) = stretcher *                             &
           (q_latest(i,j,(k+1)) - q_latest(i,j,(k-1)))
            qcf_large(i,j,(k_index+1)) = stretcher *                           &
           (qcf_latest(i,j,(k+1)) - qcf_latest(i,j,(k-1)))
          else
            !               Interpolate from level k-1 to k
            stretcher = delta_p /                                              &
           (p_layer_centres(i,j,k-1) - p_large(i,j,k_index))
            t_large(i,j,(k_index+1)) = stretcher *                             &
           (t_large(i,j,k_index) - t_latest(i,j,(k-1)))
            q_large(i,j,(k_index+1)) = stretcher *                             &
           (q_large(i,j,k_index) - q_latest(i,j,(k-1)))
            qcf_large(i,j,(k_index+1)) = stretcher *                           &
           (qcf_large(i,j,k_index) - qcf_latest(i,j,(k-1)))
          end if

        else
          if ( linked(i,j,k) ) then
            !               Interpolate from level k to k+1
            stretcher = delta_p /                                              &
            (p_large(i,j,k_index) - p_layer_centres(i,j,k+1))
            t_large(i,j,(k_index+1)) = stretcher *                             &
           (t_latest(i,j,(k+1)) - t_large(i,j,k_index))
            q_large(i,j,(k_index+1)) = stretcher *                             &
           (q_latest(i,j,(k+1)) - q_large(i,j,k_index))
            qcf_large(i,j,(k_index+1)) = stretcher *                           &
           (qcf_latest(i,j,(k+1)) - qcf_large(i,j,k_index))
          else
            !               No interpolation, freeze at level k
            t_large(i,j,(k_index+1)) = 0.0
            q_large(i,j,(k_index+1)) = 0.0
            qcf_large(i,j,(k_index+1)) = 0.0
          end if

        end if

        ! Protect against q or qcf going negative (T would imply blow-up anyway)
        if (q_large(i,j,k_index)  <                                            &
                   (abs(q_large(i,j,(k_index+1)))+tol_test))                   &
                        q_large(i,j,(k_index+1)) = 0.0
        if (qcf_large(i,j,k_index)  <                                          &
                     (abs(qcf_large(i,j,(k_index+1)))+tol_test))               &
                          qcf_large(i,j,(k_index+1)) = 0.0

        ! Select variable values at level below layer centre
        t_large(i,j,(k_index-1)) = t_large(i,j,k_index) -                      &
                                   t_large(i,j,(k_index+1))
        q_large(i,j,(k_index-1)) = q_large(i,j,k_index) -                      &
                                   q_large(i,j,(k_index+1))
        qcf_large(i,j,(k_index-1))=qcf_large(i,j,k_index)-                     &
                                   qcf_large(i,j,(k_index+1))

        ! Select variable values at level above layer centre
        ! NB: CEASE using X_large(i,j,(k_index+1)) as store for X increments.
        t_large(i,j,(k_index+1)) = t_large(i,j,(k_index+1)) +                  &
                                   t_large(i,j,k_index)
        q_large(i,j,(k_index+1)) = q_large(i,j,(k_index+1)) +                  &
                                   q_large(i,j,k_index)
        qcf_large(i,j,(k_index+1))=qcf_large(i,j,(k_index+1)) +                &
                                   qcf_large(i,j,k_index)
      end do
    end do

  end do
!$OMP end do

  ! Create an array of NTML values adjusted to the large levels
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ntml_large(i,j) = 3+(levels_per_level*(ntml(i,j)-1))
    end do
  end do
!$OMP end PARALLEL

  call ls_cld( p_large, rhcrit_large,                                          &
               large_levels, bl_levels,                                        &
               rhc_row_length, rhc_rows,                                       &
               ntml_large, cumulus, l_mixing_ratio,                            &
               t_large, cloud_fraction_large,                                  &
               q_large, qcf_large, qcl_large,                                  &
               cloud_fraction_liquid_large,                                    &
               cloud_fraction_frozen_large, error_code )

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(levels_per_level, area_cloud_fraction,                            &
!$OMP        cloud_fraction_large, bulk_cloud_fraction,                        &
!$OMP        inverse_level, qcl_latest, qcl_large,                             &
!$OMP        cloud_fraction_liquid, cloud_fraction_liquid_large,               &
!$OMP        cloud_fraction_frozen, cloud_fraction_frozen_large,               &
!$OMP        q_latest, t_latest, tdims, lcrcp, t_large, q_large,               &
!$OMP        large_levels )                                                    &
!$OMP private(i, j, k, k_index)
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      t_latest(i,j,1) = t_large(i,j,1)
      q_latest(i,j,1) = q_large(i,j,1)

      area_cloud_fraction(i,j,1) = cloud_fraction_large(i,j,1)
      bulk_cloud_fraction(i,j,1) = cloud_fraction_large(i,j,1)

      qcl_latest(i,j,1) = qcl_large(i,j,1)
      cloud_fraction_liquid(i,j,1) =                                           &
        cloud_fraction_liquid_large(i,j,1)
      cloud_fraction_frozen(i,j,1) =                                           &
          cloud_fraction_frozen_large(i,j,1)

      t_latest(i,j,tdims%k_end) = t_large(i,j,large_levels)
      q_latest(i,j,tdims%k_end) = q_large(i,j,large_levels)

      area_cloud_fraction(i,j,tdims%k_end) =                                   &
          cloud_fraction_large(i,j,large_levels)
      bulk_cloud_fraction(i,j,tdims%k_end) =                                   &
          cloud_fraction_large(i,j,large_levels)

      qcl_latest(i,j,tdims%k_end)=qcl_large(i,j,large_levels)
      cloud_fraction_liquid(i,j,tdims%k_end) =                                 &
          cloud_fraction_liquid_large(i,j,large_levels)
      cloud_fraction_frozen(i,j,tdims%k_end) =                                 &
          cloud_fraction_frozen_large(i,j,large_levels)
    end do
  end do
!$OMP end do

  ! Output variables for remaining layers

!$OMP do SCHEDULE(STATIC, 1)
  do k = 2, (tdims%k_end - 1)
    k_index = 3 + (levels_per_level * (k-2))

    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        ! Area cloud fraction is maximum of sub-layer cloud fractions
        area_cloud_fraction(i,j,k) =                                           &
        max( cloud_fraction_large(i,j,k_index),                                &
             (max(cloud_fraction_large(i,j,(k_index+1)),                       &
                  cloud_fraction_large(i,j,(k_index-1)))) )

        ! Bulk cloud fraction is mean of sub-layer cloud fractions : strictly
        ! this is a pressure weighted mean being used to approximate a volume
        ! mean. Over the depth of a layer the difference should not be large.
        bulk_cloud_fraction(i,j,k) = inverse_level *                           &
           ( cloud_fraction_large(i,j,(k_index-1)) +                           &
             cloud_fraction_large(i,j, k_index)    +                           &
             cloud_fraction_large(i,j,(k_index+1)) )

        ! The pressure weighted mean of qcf is the input qcf: do not update.

        ! Qcl is the pressure weighted mean of qcl from each sub-layer.
        qcl_latest(i,j,k) = inverse_level *                                    &
      ( qcl_large(i,j,(k_index-1)) +                                           &
        qcl_large(i,j,k_index) + qcl_large(i,j,(k_index+1)) )

        ! Liq. cloud fraction is mean of sub-layer cloud fractions : strictly
        ! this is a pressure weighted mean being used to approximate a volume
        ! mean. Over the depth of a layer the difference should not be large.
        cloud_fraction_liquid(i,j,k) = inverse_level *                         &
      ( cloud_fraction_liquid_large(i,j,(k_index-1)) +                         &
        cloud_fraction_liquid_large(i,j,k_index)     +                         &
        cloud_fraction_liquid_large(i,j,(k_index+1)) )

        ! Froz cloud fraction is mean of sub-layer cloud fractions : strictly
        ! this is a pressure weighted mean being used to approximate a volume
        ! mean. Over the depth of a layer the difference should not be large.
        cloud_fraction_frozen(i,j,k) = inverse_level *                         &
      ( cloud_fraction_frozen_large(i,j,(k_index-1)) +                         &
        cloud_fraction_frozen_large(i,j,k_index)     +                         &
        cloud_fraction_frozen_large(i,j,(k_index+1)) )

        ! Transform q_latest from qT(vapour + liquid) to specific humidity.
        ! Transform T_latest from TL(vapour + liquid) to temperature.
        q_latest(i,j,k) = q_latest(i,j,k) - qcl_latest(i,j,k)
        t_latest(i,j,k) = t_latest(i,j,k) +                                    &
                            (qcl_latest(i,j,k) * lcrcp)
      end do
    end do

  end do
!$OMP end do
!$OMP end PARALLEL

else if (i_cld_area == acf_brooks) then

  !       As before, update variables that would have been updated
  !       without area cloud fraction on
  call ls_cld( p_layer_centres(1,1,1), rhcrit,                                 &
             tdims%k_end, bl_levels,                                           &
             rhc_row_length, rhc_rows,                                         &
             ntml, cumulus, l_mixing_ratio,                                    &
             t_latest, bulk_cloud_fraction,                                    &
             q_latest, qcf_latest, qcl_latest,                                 &
             cloud_fraction_liquid,                                            &
             cloud_fraction_frozen, error_code )

  !      Calculate the area cloud fraction

  call ls_acf_brooks (                                                         &
           fv_cos_theta_latitude,                                              &
           bulk_cloud_fraction, cloud_fraction_liquid,                         &
           cloud_fraction_frozen, cumulus,                                     &
           area_cloud_fraction )

end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_arcld
! ======================================================================
end module ls_arcld_mod
