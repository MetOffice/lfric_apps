! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module proc_incr_mod

implicit none

contains

! Subroutine computes the increment to a hydrometeor species
! mixing ratio due to some process (either condensation
! / evaporation, or melting) after we've already computed linear
! coefficients for the process rate and implicitly solved for
! T and q using those coefficients

! Note: the input process rate coefficients (coefs) are assumed
! to have already been scaled by the timestep, so that this
! routine yeilds an increment, not a tendency.

subroutine proc_incr( n_points, nc, index_ic, l_full_do,                       &
                      imp_temp, imp_q_vap,                                     &
                      coefs, dq_proc )

use comorph_constants_mod, only: real_cvprec
use phase_change_coefs_mod, only: n_coefs, i_q, i_t, i_0

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where calculations actually needed
integer, intent(in) :: nc
! Indices of those points
integer, intent(in) :: index_ic(nc)
! Flag for whether to use full do-loops instead of indirect indexing
logical, intent(in) :: l_full_do

! Estimated temperature and water vapour mixing ratio
! after implicit phase changes (minus reference values)
real(kind=real_cvprec), intent(in) :: imp_temp(n_points)
real(kind=real_cvprec), intent(in) :: imp_q_vap(n_points)

! Coefficients in the process rate formula;
! these are for either the condensation/evaporation rate
! or the melting rate
real(kind=real_cvprec), intent(in) :: coefs( n_points, n_coefs )

! Amount of water species mixing ratio transferred, either from
! vapour to liquid or ice (for condensation), or from ice to
! liquid (for melting)
real(kind=real_cvprec), intent(in out) :: dq_proc(n_points)
! (intent inout to preserve values at points where no
!  calculation is done in compressed case)

! Loop counter
integer :: ic, ic2

! If calculation required at majority of points
if ( l_full_do ) then
  ! Just do calculation at all points

  do ic = 1, n_points
    dq_proc(ic) =                                                              &
            coefs(ic,i_q) * imp_q_vap(ic)                                      &
          + coefs(ic,i_t) * imp_temp(ic)                                       &
          + coefs(ic,i_0)
  end do

  ! If calculations only needed at a minority of points
else
  ! Do indirect indexing version of the same calculation

  do ic2 = 1, nc
    ic = index_ic(ic2)
    dq_proc(ic) =                                                              &
            coefs(ic,i_q) * imp_q_vap(ic)                                      &
          + coefs(ic,i_t) * imp_temp(ic)                                       &
          + coefs(ic,i_0)
  end do

end if

return
end subroutine proc_incr


end module proc_incr_mod
