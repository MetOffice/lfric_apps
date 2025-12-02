! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_conv_incs_mod

use um_types, only: real_umphys

implicit none

! Routine below either saves values before convection, or subtracts values
! before to calculate increments, depending on the input switch i_call;
! the 2 allowed values of this switch are stored here:
integer, parameter :: i_call_save_before_conv = 1
integer, parameter :: i_call_diff_to_get_incs = 2

contains

! Routine to calculate convection increments using the latest values minus
! values before convection saved in the increment arrays
subroutine calc_conv_incs( i_call, l_conv_inc_w, z_theta, z_rho,               &
                           u_p, v_p, ustar_p, vstar_p, w, w_work,              &
                           u_th_n, v_th_n, u_th_np1, v_th_np1,                 &
                           theta_star, q_star, qcl_star, qcf_star,             &
                           qcf2_star, qrain_star, qgraup_star,                 &
                           cf_liquid_star, cf_frozen_star, bulk_cf_star,       &
                           dubydt_pout, dvbydt_pout, r_w,                      &
                           theta_inc, q_inc, qcl_inc, qcf_inc,                 &
                           qcf2_inc, qrain_inc, qgraup_inc,                    &
                           cf_liquid_inc, cf_frozen_inc, bulk_cf_inc )

use atm_fields_bounds_mod, only: tdims, pdims, pdims_s, wdims, wdims_s
use mphys_inputs_mod, only: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup
use cloud_inputs_mod, only: i_cld_vn
use pc2_constants_mod, only: i_cld_pc2
use timestep_mod, only: recip_timestep

implicit none

! Integer switch indicating whether this is the call to save fields before
! convection, or subtract the values before from after to get increments
integer, intent(in) :: i_call

! Switch for whether convection should increment w
logical, intent(in) :: l_conv_inc_w

! Model-level heights above surface
real(kind=real_umphys), intent(in) ::                                          &
                    z_theta        ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    z_rho          ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     pdims%k_start:pdims%k_end )

! u and v winds at start-of-timestep, interpolated to p-grid (on rho-levels)
real(kind=real_umphys), intent(in) ::                                          &
                    u_p            ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     pdims%k_start:pdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    v_p            ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     pdims%k_start:pdims%k_end )
! u and v winds before convection, interpolated to p-grid (on rho-levels)
real(kind=real_umphys), intent(in) ::                                          &
                    ustar_p        ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     pdims%k_start:pdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    vstar_p        ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     pdims%k_start:pdims%k_end )
! Vertical velocity w at start-of-timestep
real(kind=real_umphys), intent(in) ::                                          &
                    w              ( wdims_s%i_start:wdims_s%i_end,            &
                                     wdims_s%j_start:wdims_s%j_end,            &
                                     wdims_s%k_start:wdims_s%k_end )
! Temporary work array storing latest w if convection not updating the real w
real(kind=real_umphys), intent(in out) ::                                      &
                           w_work     ( wdims%i_start:wdims%i_end,             &
                                        wdims%j_start:wdims%j_end,             &
                                        1:wdims%k_end )

! u and v winds interpolated to theta-levels, at start-of-timestep
real(kind=real_umphys), intent(out) ::                                         &
                     u_th_n        ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end-1)
real(kind=real_umphys), intent(out) ::                                         &
                     v_th_n        ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end-1)
! u and v winds interpolated to theta-levels, updated by convection
real(kind=real_umphys), intent(in out) ::                                      &
                        u_th_np1   ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end-1)
real(kind=real_umphys), intent(in out) ::                                      &
                        v_th_np1   ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end-1)

! Latest vales of temperature and moisture fields
real(kind=real_umphys), intent(in) ::                                          &
                    theta_star     ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    q_star         ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qcl_star       ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qcf_star       ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qcf2_star      ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qrain_star     ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qgraup_star    ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    cf_liquid_star ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    cf_frozen_star ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    bulk_cf_star   ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )

! Convective u,v tendencies on rho-levels, output to atmos_physics2
! to update u,v on their native staggered grids
real(kind=real_umphys), intent(out) ::                                         &
                     dubydt_pout   ( pdims_s%i_start:pdims_s%i_end,            &
                                     pdims_s%j_start:pdims_s%j_end,            &
                                     pdims_s%k_start:pdims_s%k_end )
real(kind=real_umphys), intent(out) ::                                         &
                     dvbydt_pout   ( pdims_s%i_start:pdims_s%i_end,            &
                                     pdims_s%j_start:pdims_s%j_end,            &
                                     pdims_s%k_start:pdims_s%k_end )
! Total vertical velocity increment so far this timestep
real(kind=real_umphys), intent(in out) ::                                      &
                        r_w        ( wdims%i_start:wdims%i_end,                &
                                     wdims%j_start:wdims%j_end,                &
                                     1:wdims%k_end )

! Convection increment arrays
real(kind=real_umphys), intent(in out) ::                                      &
                        theta_inc     ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        q_inc         ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        qcl_inc       ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        qcf_inc       ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        qcf2_inc      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        qrain_inc     ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        qgraup_inc    ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        cf_liquid_inc ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        cf_frozen_inc ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        bulk_cf_inc   ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )

! Interpolation weight
real(kind=real_umphys) :: interp

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(none) private( i, j, k, interp )                        &
!$OMP SHARED( i_call, tdims, wdims, pdims, z_theta, z_rho,                     &
!$OMP         u_th_n, u_p, v_th_n, v_p, u_th_np1, ustar_p, v_th_np1, vstar_p,  &
!$OMP         l_conv_inc_w, r_w, w, w_work, theta_inc, theta_star,             &
!$OMP         q_inc, q_star, qcl_inc, qcl_star, qcf_inc, qcf_star,             &
!$OMP         l_mcr_qcf2, qcf2_inc, qcf2_star,                                 &
!$OMP         l_mcr_qrain, qrain_inc, qrain_star,                              &
!$OMP         l_mcr_qgraup, qgraup_inc, qgraup_star,                           &
!$OMP         i_cld_vn, cf_liquid_inc, cf_liquid_star,                         &
!$OMP         cf_frozen_inc, cf_frozen_star, bulk_cf_inc, bulk_cf_star,        &
!$OMP         dubydt_pout, dvbydt_pout, recip_timestep )

! Which call to this routine are we in?
select case (i_call)

case (i_call_save_before_conv)
  ! 1st call: save fields before convection...

  ! Horizontal winds are on rho-levels, but CoMorph needs them
  ! to be colocated with the other fields.  So make copies
  ! interpolated to theta-levels
!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end-1
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        interp = ( z_theta(i,j,k) - z_rho(i,j,k) )                             &
               / ( z_rho(i,j,k+1) - z_rho(i,j,k) )

        u_th_n(i,j,k) = (1.0-interp) * u_p(i,j,k)                              &
                      +      interp  * u_p(i,j,k+1)
        v_th_n(i,j,k) = (1.0-interp) * v_p(i,j,k)                              &
                      +      interp  * v_p(i,j,k+1)

        u_th_np1(i,j,k) = (1.0-interp) * ustar_p(i,j,k)                        &
                        +      interp  * ustar_p(i,j,k+1)
        v_th_np1(i,j,k) = (1.0-interp) * vstar_p(i,j,k)                        &
                        +      interp  * vstar_p(i,j,k+1)
      end do
    end do
  end do
!$OMP end do NOWAIT

  if ( l_conv_inc_w ) then

    ! Convert increment r_w to field of w with increment so far added on
!$OMP do SCHEDULE(STATIC)
    do k = 1, wdims%k_end
      do j = wdims%j_start, wdims%j_end
        do i = wdims%i_start, wdims%i_end
          r_w(i,j,k) = r_w(i,j,k) + w(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT

  else

    ! Make a separate work array for w if we don't want convection
    ! to actually modify the UM's w-field.
!$OMP do SCHEDULE(STATIC)
    do k = 1, wdims%k_end
      do j = wdims%j_start, wdims%j_end
        do i = wdims%i_start, wdims%i_end
          w_work(i,j,k) = r_w(i,j,k) + w(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT

  end if

  ! Save values of temperature and moisture fields before convection
  ! in the increment arrays...

!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        theta_inc(i,j,k) = theta_star(i,j,k)
        q_inc(i,j,k)     = q_star(i,j,k)
        qcl_inc(i,j,k)   = qcl_star(i,j,k)
        qcf_inc(i,j,k)   = qcf_star(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT

  if ( l_mcr_qcf2 ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qcf2_inc(i,j,k) = qcf2_star(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if
  if ( l_mcr_qrain ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qrain_inc(i,j,k) = qrain_star(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if
  if ( l_mcr_qgraup ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qgraup_inc(i,j,k) = qgraup_star(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if

  if ( i_cld_vn == i_cld_pc2 ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          cf_liquid_inc(i,j,k) = cf_liquid_star(i,j,k)
          cf_frozen_inc(i,j,k) = cf_frozen_star(i,j,k)
          bulk_cf_inc(i,j,k)   = bulk_cf_star(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if

case (i_call_diff_to_get_incs)
  ! 2nd call: subtract values before convection to get increments...

  ! Convert final u,v to convection u,v increments on theta-levels
!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end-1
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        interp = ( z_theta(i,j,k) - z_rho(i,j,k) )                             &
               / ( z_rho(i,j,k+1) - z_rho(i,j,k) )

        ! Re-interpolate from fields on rho-levels to recover
        ! the winds before convection increments were added on,
        ! and subtract those to get increments
        u_th_np1(i,j,k) = u_th_np1(i,j,k)                                      &
                        - ( (1.0-interp) * ustar_p(i,j,k)                      &
                          +      interp  * ustar_p(i,j,k+1) )
        v_th_np1(i,j,k) = v_th_np1(i,j,k)                                      &
                        - ( (1.0-interp) * vstar_p(i,j,k)                      &
                          +      interp  * vstar_p(i,j,k+1) )
      end do
    end do
  end do
!$OMP end do

  ! Interpolate increments onto rho-levels for output
!$OMP do SCHEDULE(STATIC)
  do k = 2, pdims%k_end-1
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        interp = ( z_rho(i,j,k)   - z_theta(i,j,k-1) )                         &
               / ( z_theta(i,j,k) - z_theta(i,j,k-1) )

        dubydt_pout(i,j,k) = (1.0-interp) * u_th_np1(i,j,k-1)                  &
                           +      interp  * u_th_np1(i,j,k)
        dvbydt_pout(i,j,k) = (1.0-interp) * v_th_np1(i,j,k-1)                  &
                           +      interp  * v_th_np1(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
  ! Increment bottom rho-level using bottom theta-level increment
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      dubydt_pout(i,j,1) = u_th_np1(i,j,1)
      dvbydt_pout(i,j,1) = v_th_np1(i,j,1)
    end do
  end do
!$OMP end do NOWAIT
  ! Interpolate top rho-level assuming increments go to zero
  ! at the model-top
  k = pdims%k_end
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      interp = ( z_rho(i,j,k)   - z_theta(i,j,k-1) )                           &
             / ( z_theta(i,j,k) - z_theta(i,j,k-1) )

      dubydt_pout(i,j,k) = (1.0-interp) * u_th_np1(i,j,k-1)
      dvbydt_pout(i,j,k) = (1.0-interp) * v_th_np1(i,j,k-1)
    end do
  end do
!$OMP end do

  ! Convert u,v increments to tendencies, as this is what atmos_physics2
  ! expects
!$OMP do SCHEDULE(STATIC)
  do k = 1, pdims%k_end
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        dubydt_pout(i,j,k) = dubydt_pout(i,j,k) * recip_timestep
        dvbydt_pout(i,j,k) = dvbydt_pout(i,j,k) * recip_timestep
      end do
    end do
  end do
!$OMP end do NOWAIT

  if ( l_conv_inc_w ) then

    ! Vertical winds; convert updated final w back into increment
!$OMP do SCHEDULE(STATIC)
    do k = 1, wdims%k_end
      do j = wdims%j_start, wdims%j_end
        do i = wdims%i_start, wdims%i_end
          r_w(i,j,k) = r_w(i,j,k) - w(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT

    ! If not incrementing w, just discard the w updated by CoMorph.
  end if

  ! Difference latest temperature and moisture fields with saved values from
  ! before convection, to compute increments
!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        theta_inc(i,j,k) = theta_star(i,j,k) - theta_inc(i,j,k)
        q_inc(i,j,k)     = q_star(i,j,k)     - q_inc(i,j,k)
        qcl_inc(i,j,k)   = qcl_star(i,j,k)   - qcl_inc(i,j,k)
        qcf_inc(i,j,k)   = qcf_star(i,j,k)   - qcf_inc(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT

  if ( l_mcr_qcf2 ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qcf2_inc(i,j,k) = qcf2_star(i,j,k) - qcf2_inc(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if
  if ( l_mcr_qrain ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qrain_inc(i,j,k) = qrain_star(i,j,k) - qrain_inc(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if
  if ( l_mcr_qgraup ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qgraup_inc(i,j,k) = qgraup_star(i,j,k) - qgraup_inc(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if

  if ( i_cld_vn == i_cld_pc2 ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          cf_liquid_inc(i,j,k) = cf_liquid_star(i,j,k) - cf_liquid_inc(i,j,k)
          cf_frozen_inc(i,j,k) = cf_frozen_star(i,j,k) - cf_frozen_inc(i,j,k)
          bulk_cf_inc(i,j,k)   = bulk_cf_star(i,j,k)   - bulk_cf_inc(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if

end select  ! case(i_call)

!$OMP end PARALLEL

return
end subroutine calc_conv_incs

end module calc_conv_incs_mod
