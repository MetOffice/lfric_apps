! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module fall_out_mod

implicit none

contains


! Routine to compute fall-out of hydrometeor from the parcel.
! Yields a simple implicit solution, where the fall-out is
! proportional to the input fall-speed and the mixing ratio
! of hydrometeor remaining after fall-out has been removed.
subroutine fall_out( n_points, wf_cond, dt_over_lz,                            &
                     q_cond_in, q_cond_out )

use comorph_constants_mod, only: real_cvprec, one

implicit none

! Number of points
integer, intent(in) :: n_points

! Fall-speed of the hydrometeor particles
real(kind=real_cvprec), intent(in) :: wf_cond(n_points)

! Factor delta_t / vert_len which scales fall-out increment
real(kind=real_cvprec), intent(in) :: dt_over_lz(n_points)

! Mixing ratio of hydrometeor before fall-out
real(kind=real_cvprec), intent(in) :: q_cond_in(n_points)
! Mixing ratio of hydrometeor after fall-out
real(kind=real_cvprec), intent(out) :: q_cond_out(n_points)

! Loop counter
integer :: ic

! Implicit time-discretization for fall-out:
! q_cond_in is the value of q_cond before fall-out.
! q_cond_out is the final value of q_cond after fall-out has
! been removed.
!
! Backward-time discretization; fall-out tendency depends on
! value after the increment has been subtracted:
! q_cond_fall = q_cond_out wf/lz dt
! q_cond_out = q_cond_in - q_cond_fall
! q_cond_out = q_cond_in - q_cond_out wf/lz dt
! q_cond_out ( 1 + wf/lz dt ) = q_cond_in
!
! therefore
!
! q_cond_out = q_cond_in / ( 1 + wf/lz dt )

do ic = 1, n_points
  ! Calculate q_cond after fall-out using implicit formula
  q_cond_out(ic) = q_cond_in(ic)                                               &
                 / ( one + wf_cond(ic) * dt_over_lz(ic) )
end do

return
end subroutine fall_out



! Alternative version which also computes the fall-out flux
! in kg m-2 s-1.  It also optionally compresses the
! calculation, and udpates the input value of q_cond
! instead of using a separate array for the output.
subroutine fall_out_flux( n_points, nc, index_ic,                              &
                          wf_cond, delta_t, vert_len,                          &
                          dt_over_rhod_lz,                                     &
                          q_cond, flux_cond,                                   &
                          l_diags, i_cond, moist_proc_diags,                   &
                          n_points_diag, n_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, zero, one, indi_thresh
use moist_proc_diags_type_mod, only: moist_proc_diags_type

implicit none

! Number of points
integer, intent(in) :: n_points

! List of points for compressed calculation
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Fall-speed of the hydrometeor particles
real(kind=real_cvprec), intent(in) :: wf_cond(n_points)

! Time interval
real(kind=real_cvprec), intent(in) :: delta_t(n_points)
! Vertical length of the parcel / level
real(kind=real_cvprec), intent(in) :: vert_len(n_points)
! Factor delta_t / ( rho_dry vert_len ) used to convert
! fall-out increment into fall-out flux for output
real(kind=real_cvprec), intent(in) :: dt_over_rhod_lz(n_points)

! Mixing ratio of hydrometeor in the parcel
! (gets updated with negative increment from fall-out)
real(kind=real_cvprec), intent(in out) :: q_cond(n_points)

! Fall-flux of hydrometeor out of the parcel / kg m-2 s-1
real(kind=real_cvprec), intent(out) :: flux_cond(n_points)

! Master switch for whether or not to calculate any diagnostics
logical, intent(in) :: l_diags
! Super-array address of current species,
! needed to address the correct diagnostic array field
integer, intent(in) :: i_cond
! Structure containing diagnostic flags etc.
type(moist_proc_diags_type), intent(in) :: moist_proc_diags
! Number of points in the diagnostics super-array
integer, intent(in) :: n_points_diag
! Total number of diagnostics in the super-array
integer, intent(in) :: n_diags
! Super-array to store all the diagnostics
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                         ( n_points_diag, n_diags )

! Loop counters
integer :: ic, ic2, i_super


! If calculation required at majority of points
if ( real(nc,real_cvprec) > indi_thresh * real(n_points,real_cvprec) ) then
  ! Full-field calculation

  do ic = 1, n_points
    ! Calculate increment due to fall-out
    flux_cond(ic) = q_cond(ic) * ( one                                         &
     - one / ( one + wf_cond(ic) * delta_t(ic) / vert_len(ic) ) )
      ! (flux_cond now stores the mixing-ratio increment)
  end do

  ! Apply increment to q_cond
  do ic = 1, n_points
    q_cond(ic) = q_cond(ic) - flux_cond(ic)
  end do

  if ( l_diags ) then
    ! Update sedimentation increment diagnostic, if requested
    if ( moist_proc_diags % diags_cond(i_cond)%pt                              &
         % dq_fall % flag ) then
      i_super = moist_proc_diags % diags_cond(i_cond)%pt                       &
                % dq_fall % i_super
      do ic = 1, n_points
        diags_super(ic,i_super) = diags_super(ic,i_super)                      &
                                - flux_cond(ic)
      end do
    end if
  end if

  do ic = 1, n_points
    ! Convert increment to fall-flux
    flux_cond(ic) = flux_cond(ic) / dt_over_rhod_lz(ic)
  end do

else
  ! Compressed version of the same calculation

  ! First fall-flux needs to be set to zero at points that
  ! aren't in the compression list
  do ic = 1, n_points
    flux_cond(ic) = zero
  end do

  do ic2 = 1, nc
    ic = index_ic(ic2)
    flux_cond(ic) = q_cond(ic) * ( one                                         &
     - one / ( one + wf_cond(ic) * delta_t(ic) / vert_len(ic) ) )
    q_cond(ic) = q_cond(ic) - flux_cond(ic)
  end do

  if ( l_diags ) then
    if ( moist_proc_diags % diags_cond(i_cond)%pt                              &
         % dq_fall % flag ) then
      i_super = moist_proc_diags % diags_cond(i_cond)%pt                       &
                % dq_fall % i_super
      do ic2 = 1, nc
        ic = index_ic(ic2)
        diags_super(ic,i_super) = diags_super(ic,i_super)                      &
                                - flux_cond(ic)
      end do
    end if
  end if

  do ic2 = 1, nc
    ic = index_ic(ic2)
    flux_cond(ic) = flux_cond(ic) / dt_over_rhod_lz(ic)
  end do

end if


return
end subroutine fall_out_flux


end module fall_out_mod
