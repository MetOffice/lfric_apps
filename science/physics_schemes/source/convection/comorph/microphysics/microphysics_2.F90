! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module microphysics_2_mod

implicit none

contains

! Subroutine to do final calculations for BOTEMS
! (Back-Of-The-Envelope Microphysical Scheme)
!
! This includes any microphysical processes that need to be done
! after the implicit solution of phase-changes.
! Current only does autoconversion of liquid cloud to rain
subroutine microphysics_2( n_points, n_points_super, nc, index_ic,             &
                           delta_t, vert_len, wf_cond, q_cond,                 &
                           l_diags, moist_proc_diags,                          &
                           n_points_diag, n_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, zero, one, two, four,            &
                                 n_cond_species, i_cond_cl, i_cond_rain,       &
                                 indi_thresh, autoc_opt, autoc_linear,         &
                                 autoc_quadratic, q_cl_auto, coef_auto
use moist_proc_diags_type_mod, only: moist_proc_diags_type

implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points in the condensed water species super-array
! (this might be reused on subsequent calls with bigger size
!  than needed, to save having to re-allocate it)
integer, intent(in) :: n_points_super

! Number of points where each hydrometeor species is non-zero
integer, intent(in) :: nc(n_cond_species)
! Indices of those points
integer, intent(in) :: index_ic(n_points,n_cond_species)

! Time interval for converting process rates to increments.
real(kind=real_cvprec), intent(in) :: delta_t(n_points)

! Vertical length-scale of the parcel.
real(kind=real_cvprec), intent(in) :: vert_len(n_points)

! Fall-speed of each hydrometeor species
real(kind=real_cvprec), intent(in) :: wf_cond                                  &
                                   ( n_points, n_cond_species )

! Mixing ratios of condensed water species
real(kind=real_cvprec), intent(in out) :: q_cond                               &
                             ( n_points_super, n_cond_species )

! Master switch for whether or not to calculate any diagnostics
logical, intent(in) :: l_diags
! Structure containing diagnostic flags etc.
type(moist_proc_diags_type), intent(in) :: moist_proc_diags
! Number of points in the diagnostics super-array
integer, intent(in) :: n_points_diag
! Total number of diagnostics in the super-array
integer, intent(in) :: n_diags
! Super-array to store all the diagnostics
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                         ( n_points_diag, n_diags )


! Increment to q_rain due to autoconversion
real(kind=real_cvprec) :: dq_auto(n_points)

! Temporary store for ratio fall-speed / height-interval
real(kind=real_cvprec) :: wf_over_lz

! Coefficients of quadratic equation
real(kind=real_cvprec) :: a_quad, b_quad, c_quad

! Solution to quadratic equation for final q_cl
real(kind=real_cvprec) :: q_cl_out

! Loop counters
integer :: ic, ic2, i_super


! Autoconversion of liquid cloud to rain

! If any liquid cloud present
if ( nc(i_cond_cl) > 0 ) then

  ! If calculation required at majority of points
  if ( real(nc(i_cond_cl),real_cvprec) > indi_thresh                           &
                                         * real(n_points,real_cvprec) ) then
    ! Full-field calculation

    ! Implicit discretisation of autoconversion and fall-out:
    ! Two choices, linear or quadratic
    if ( autoc_opt==autoc_linear ) then
      ! ----------------------------------------------------------------
      ! Linear, with threshold
      !
      ! d_q_cl_auto = delta_t coef_auto ( q_cl_out - q_cl_auto )
      ! d_q_cl_fall = delta_t q_cl_out wf/lz
      !
      ! Therefore
      !
      ! q_cl_out = q_cl_in - delta_t (
      !                coef_auto ( q_cl_out - q_cl_auto )
      !              + q_cl_out wf/lz   )
      !
      ! Rearranging:
      !
      ! q_cl_out - q_cl_auto = q_cl_in - q_cl_auto - delta_t (
      !                coef_auto ( q_cl_out - q_cl_auto )
      !              + ( q_cl_out - q_cl_auto + q_cl_auto ) wf/lz   )
      !
      ! => ( q_cl_out - q_cl_auto ) ( 1 + delta_t ( coef_auto + wf/lz ) )
      !     = ( q_cl_in - q_cl_auto ) - delta_t wf/lz q_cl_auto
      !
      ! => ( q_cl_out - q_cl_auto )
      !    = ( ( q_cl_in - q_cl_auto ) - delta_t wf/lz q_cl_auto )
      !    / ( 1 + delta_t ( coef_auto + wf/lz ) )
      !
      ! Then multiply this by delta_t * coef_auto to get
      ! the autoconversion increment
      !
      ! Note: the solution only makes sense where the numerator
      ! ( q_cl_in - q_cl_auto ) - delta_t wf/lz q_cl_auto
      ! is positive (otherwise, the solution with no autoconversion
      ! yields q_cl < q_cl_auto, and we should have no autoconversion
      ! as we're below the threshold).

      do ic = 1, n_points
        ! Store ratio fall-speed / height-interval
        wf_over_lz = wf_cond(ic,i_cond_cl) / vert_len(ic)
        ! Compute autoconversion increment, using the above formula
        dq_auto(ic) = delta_t(ic) * coef_auto                                  &
          * max( ( q_cond(ic,i_cond_cl) - q_cl_auto )                          &
               - delta_t(ic) * wf_over_lz * q_cl_auto, zero )                  &
          / ( one + delta_t(ic) * ( coef_auto + wf_over_lz ) )
      end do

    else if ( autoc_opt==autoc_quadratic ) then
      ! -----------------------------------------------------------
      ! Quadratic, no threshold
      !
      ! d_q_cl_auto = delta_t * coef_auto * q_cl_out^2
      ! d_q_cl_fall = delta_t * q_cl_out * wf/lz
      !
      ! Therefore
      !
      ! q_cl_out = q_cl_in - delta_t ( coef_auto q_cl_out^2
      !                                  + q_cl_out wf/lz   )
      !
      ! Rearranging:
      !
      ! delta_t*coef_auto * q_cl_out^2 + (1+delta_t*wf/lz) * q_cl_out
      !                                - q_cl_in = 0
      !
      ! => quadratic in q_cl_out = (-b +/- sqrt(b^2-4ac) )/2a
      ! Only +ve root will be +ve
      ! => q_cl_out = ( -(1+delta_t*wf/lz) +
      !          sqrt( (1+delta_t*wf/lz)^2 + 4*delta_t*coef_auto*q_cl_in ) ) /
      !          (2*delta_t*coef_auto)
      !
      ! Then autoconversion increment is delta_t coef_auto * q_cl_out^2
      !
      ! Note: the solution only makes sense where the numerator
      ! is positive.

      do ic = 1, n_points
        ! Store coefficients in the quadratic formula
        a_quad = delta_t(ic) * coef_auto
        b_quad = one + delta_t(ic) * wf_cond(ic,i_cond_cl) / vert_len(ic)
        c_quad = - q_cond(ic,i_cond_cl)
        ! Compute the solution to the quadratic
        q_cl_out = max( zero,                                                  &
                 ( -b_quad + sqrt( b_quad*b_quad - four*a_quad*c_quad ) )      &
                  /( two*a_quad ) )
        ! Find increment consistent with this
        dq_auto(ic) = delta_t(ic) * coef_auto * q_cl_out * q_cl_out
      end do

    end if  ! ( autoc_opt )

    ! Add increment to q_cl and q_rain
    do ic = 1, n_points
      q_cond(ic,i_cond_cl) = q_cond(ic,i_cond_cl)                              &
                           - dq_auto(ic)
      q_cond(ic,i_cond_rain) = q_cond(ic,i_cond_rain)                          &
                             + dq_auto(ic)
    end do

    if ( l_diags ) then
      ! Save autoconversion increment diagnostics, if requested
      if ( moist_proc_diags % diags_cond(i_cond_cl)%pt                         &
           % dq_aut % flag ) then
        ! Extract super-array address
        i_super = moist_proc_diags % diags_cond(i_cond_cl)%pt                  &
                  % dq_aut % i_super
        ! Copy increment
        do ic = 1, n_points
          diags_super(ic,i_super) = diags_super(ic,i_super)                    &
                                  - dq_auto(ic)
        end do
      end if
      if ( moist_proc_diags % diags_cond(i_cond_rain)%pt                       &
           % dq_aut % flag ) then
        ! Extract super-array address
        i_super = moist_proc_diags % diags_cond(i_cond_rain)%pt                &
                  % dq_aut % i_super
        ! Copy increment
        do ic = 1, n_points
          diags_super(ic,i_super) = diags_super(ic,i_super)                    &
                                  + dq_auto(ic)
        end do
      end if
    end if  ! ( l_diags )

    ! Liquid-cloud present at only a small fraction of points
  else
    ! Compressed version of exactly the same calculation

    if ( autoc_opt==autoc_linear ) then
      ! linear, with threshold q_cl_auto

      do ic2 = 1, nc(i_cond_cl)
        ic = index_ic(ic2,i_cond_cl)
        wf_over_lz = wf_cond(ic,i_cond_cl) / vert_len(ic)
        dq_auto(ic) = delta_t(ic) * coef_auto                                  &
          * max( ( q_cond(ic,i_cond_cl) - q_cl_auto )                          &
               - delta_t(ic) * wf_over_lz * q_cl_auto, zero )                  &
          / ( one + delta_t(ic) * ( coef_auto + wf_over_lz ) )
      end do

    else if ( autoc_opt==autoc_quadratic ) then
      ! quadratic, with no threshold q_cl

      do ic2 = 1, nc(i_cond_cl)
        ic = index_ic(ic2,i_cond_cl)
        a_quad = delta_t(ic) * coef_auto
        b_quad = one + delta_t(ic) * wf_cond(ic,i_cond_cl) / vert_len(ic)
        c_quad = - q_cond(ic,i_cond_cl)
        q_cl_out = max( zero,                                                  &
                 ( -b_quad + sqrt( b_quad*b_quad - four*a_quad*c_quad ) )      &
                  /( two*a_quad ) )
        dq_auto(ic) = delta_t(ic) * coef_auto * q_cl_out * q_cl_out
      end do

    end if  ! ( autoc_opt )

    do ic2 = 1, nc(i_cond_cl)
      ic = index_ic(ic2,i_cond_cl)
      q_cond(ic,i_cond_cl) = q_cond(ic,i_cond_cl)                              &
                           - dq_auto(ic)
      q_cond(ic,i_cond_rain) = q_cond(ic,i_cond_rain)                          &
                             + dq_auto(ic)
    end do

    if ( l_diags ) then
      if ( moist_proc_diags % diags_cond(i_cond_cl)%pt                         &
           % dq_aut % flag ) then
        i_super = moist_proc_diags % diags_cond(i_cond_cl)%pt                  &
                  % dq_aut % i_super
        do ic2 = 1, nc(i_cond_cl)
          ic = index_ic(ic2,i_cond_cl)
          diags_super(ic,i_super) = diags_super(ic,i_super)                    &
                                  - dq_auto(ic)
        end do
      end if
      if ( moist_proc_diags % diags_cond(i_cond_rain)%pt                       &
           % dq_aut % flag ) then
        i_super = moist_proc_diags % diags_cond(i_cond_rain)%pt                &
                  % dq_aut % i_super
        do ic2 = 1, nc(i_cond_cl)
          ic = index_ic(ic2,i_cond_cl)
          diags_super(ic,i_super) = diags_super(ic,i_super)                    &
                                  + dq_auto(ic)
        end do
      end if
    end if  ! ( l_diags )

  end if  ! Test whether to use compressed vs full-field calculation

end if  ! ( nc(i_cond_cl) > 0 )


return
end subroutine microphysics_2

end module microphysics_2_mod
