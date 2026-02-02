! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_qcf2_incs_mod

implicit none

! Indicators for what the subroutine below should do on each of the 3 separate
! places it is called
integer, parameter :: i_call_combine_in_qcf2 = 1
integer, parameter :: i_call_subtract_qcf = 2
integer, parameter :: i_call_repartition = 3

contains

! Routine handles the 2 UM ice mass variables on input / output to comorph
! when comorph is using just 1 ice variable and the UM is using 2 separate
! fields for "ice" and "snow".
subroutine calc_qcf2_incs( i_call,                                             &
                           m_cf, m_cf2, qcf_star, qcf2_star, qcf_inc, qcf2_inc)

use um_types, only: real_umphys
use atm_fields_bounds_mod, only: tdims, tdims_s

implicit none

! Indicator for which call to this routine this is
integer, intent(in) :: i_call

! Start-of-timestep mixing-ratios of the 2 ice categories
real(kind=real_umphys), intent(in) :: m_cf                                     &
                                      ( tdims_s%i_start:tdims_s%i_end,         &
                                        tdims_s%j_start:tdims_s%j_end,         &
                                        tdims_s%k_start:tdims_s%k_end )
real(kind=real_umphys), intent(in out) :: m_cf2                                &
                                      ( tdims_s%i_start:tdims_s%i_end,         &
                                        tdims_s%j_start:tdims_s%j_end,         &
                                        tdims_s%k_start:tdims_s%k_end )

! Latest mixing-ratios of the 2 ice categories
real(kind=real_umphys), intent(in out) :: qcf_star                             &
                                          ( tdims%i_start:tdims%i_end,         &
                                            tdims%j_start:tdims%j_end,         &
                                            1:tdims%k_end )
real(kind=real_umphys), intent(in out) :: qcf2_star                            &
                                          ( tdims%i_start:tdims%i_end,         &
                                            tdims%j_start:tdims%j_end,         &
                                            1:tdims%k_end )

! Increments to the 2 ice categories
real(kind=real_umphys), intent(in out) :: qcf_inc                              &
                                          ( tdims%i_start:tdims%i_end,         &
                                            tdims%j_start:tdims%j_end,         &
                                            1:tdims%k_end )
real(kind=real_umphys), intent(in out) :: qcf2_inc                             &
                                          ( tdims%i_start:tdims%i_end,         &
                                            tdims%j_start:tdims%j_end,         &
                                            1:tdims%k_end )

! Part of convection ice increment passed to aggregates.
real(kind=real_umphys) :: dqcf

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k, dqcf )                          &
!$OMP SHARED( tdims, i_call,                                                   &
!$OMP         m_cf, m_cf2, qcf_star, qcf2_star, qcf_inc, qcf2_inc )

select case ( i_call )

case ( i_call_combine_in_qcf2 )
  ! Put all the ice-cloud mass in the qcf2 field to pass into comorph

!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        m_cf2(i,j,k)     = m_cf2(i,j,k)     + m_cf(i,j,k)
        qcf2_star(i,j,k) = qcf2_star(i,j,k) + qcf_star(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT

case ( i_call_subtract_qcf )
  ! Subtract qcf off from qcf2 again after the call to comorph

!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        m_cf2(i,j,k)     = m_cf2(i,j,k)     - m_cf(i,j,k)
        qcf2_star(i,j,k) = qcf2_star(i,j,k) - qcf_star(i,j,k)
        if ( qcf2_star(i,j,k) < 0.0 ) then
          ! If comorph produced a negative increment to qcf+qcf2,
          ! its possible for the above subtraction to result in
          ! negative qcf2.  If this happens, subtract from qcf
          qcf_star(i,j,k) = qcf_star(i,j,k) + qcf2_star(i,j,k)
          qcf2_star(i,j,k) = 0.0
        end if
      end do
    end do
  end do
!$OMP END DO NOWAIT

case ( i_call_repartition )
  ! Calculations so far have put the convection increment to
  ! qcf_total=qcf+qcf2 entirely in the qcf2 field (crystals),
  ! and left the mass of aggregates alone.
  ! Now repartition the increment between crystals and aggregates
  ! so that we detrain only aggregates:

!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        ! Compute the convection ice mass increment to add to qcf,
        ! with appropriate limits to avoid making negative qcf or qcf2
        dqcf = max( min( qcf2_inc(i,j,k),                                      &
                         qcf2_star(i,j,k) ),                                   &
                    -qcf_star(i,j,k) )
        ! Adjust latest ice masses and increments consistently
        qcf2_star(i,j,k) = qcf2_star(i,j,k) - dqcf
        qcf2_inc(i,j,k)  = qcf2_inc(i,j,k)  - dqcf
        qcf_star(i,j,k)  = qcf_star(i,j,k)  + dqcf
        qcf_inc(i,j,k)   = qcf_inc(i,j,k)   + dqcf
      end do
    end do
  end do
!$OMP END DO NOWAIT

end select

!$OMP END PARALLEL


return
end subroutine calc_qcf2_incs

end module calc_qcf2_incs_mod
