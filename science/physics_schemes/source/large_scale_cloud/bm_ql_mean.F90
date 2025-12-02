! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud

module bm_ql_mean_mod

implicit none

character(len=*), parameter, private :: ModuleName = 'BM_QL_MEAN_MOD'

contains

! Subroutine to calculate the mass-weighted vertical mean q + qcl
! up to the current height at each model-level.  This is used in the
! construction of the height-varying minimum-allowed unimodal variance
! in the bimodal cloud-scheme
subroutine bm_ql_mean( nlevels, qt_in, p_theta_levels, ql_mean )

use um_types,              only: real_umphys
use atm_fields_bounds_mod, only: tdims, pdims
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

implicit none

! Number of model-levels
integer, intent(in) :: nlevels

! Input profile of qt = q + qcl
real(kind=real_umphys), intent(in) :: qt_in                                    &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )

! Pressure profile on theta-levels
real(kind=real_umphys), intent(in) :: p_theta_levels                           &
                                      ( pdims%i_start:pdims%i_end,             &
                                        pdims%j_start:pdims%j_end, nlevels )

! Output mass-weighted vertical mean of qt_in up to current level
real(kind=real_umphys), intent(out) :: ql_mean                                 &
                                       ( tdims%i_start:tdims%i_end,            &
                                         tdims%j_start:tdims%j_end, nlevels )

! Loop counters
integer :: i, j, k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='BM_QL_MEAN'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!$OMP PARALLEL DEFAULT(none) private(i,j,k)                                    &
!$OMP SHARED( nlevels, tdims, qt_in, ql_mean, p_theta_levels )

! Initialise to zero in first model-level
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    ql_mean(i,j,1) = 0.0
  end do
end do
!$OMP end do

! Integrate qt dp up to current model-level
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  ! j-loop is outermost as the k-loop can't be parallelised
  do k = 2, nlevels
    do i = tdims%i_start, tdims%i_end
      ql_mean(i,j,k) = ql_mean(i,j,k-1)                                        &
                     + 0.5*( qt_in(i,j,k-1) + qt_in(i,j,k) )                   &
                     * ( p_theta_levels(i,j,k-1) - p_theta_levels(i,j,k) )
    end do
  end do
end do
!$OMP end do

! Normalise by pressure difference to get (hydrostatic) mass-weighted-mean
!$OMP do SCHEDULE(STATIC)
do k = 2, nlevels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ql_mean(i,j,k) = ql_mean(i,j,k)                                          &
                     / ( p_theta_levels(i,j,1) - p_theta_levels(i,j,k) )
    end do
  end do
end do
!$OMP end do

! Set value in first model-level
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    ql_mean(i,j,1) = qt_in(i,j,1)
  end do
end do
!$OMP end do

!$OMP end PARALLEL


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine bm_ql_mean

end module bm_ql_mean_mod
