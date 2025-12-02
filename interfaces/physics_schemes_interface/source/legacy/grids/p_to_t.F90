! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine P_TO_T
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
module p_to_t_mod

implicit none

character(len=*), parameter, private :: ModuleName='P_TO_T_MOD'

contains

subroutine p_to_t                                                              &
  ( row_length, rows, halo_i, halo_j, halo_i_data, halo_j_data, levels         &
  , r_theta_levels, r_rho_levels, field_in, field_out )


! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v10.3 programming standards.
!   This routine is thread safe and can be called from an
!   OMP parallel region

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
!$ use compute_chunk_size_mod, only: compute_chunk_size  ! Note OpenMP sentinel

implicit none

! Arguments with intent in. ie: input variables.
integer, intent(in) :: row_length
integer, intent(in) :: rows
integer, intent(in) :: levels
integer, intent(in) :: halo_i
integer, intent(in) :: halo_j
integer, intent(in) :: halo_i_data
integer, intent(in) :: halo_j_data

real, intent(in)  :: r_theta_levels ( 1-halo_i:row_length+halo_i               &
                                    , 1-halo_j:rows+halo_j, 0:levels+1 )

real, intent(in)  :: r_rho_levels   ( 1-halo_i:row_length+halo_i               &
                                    , 1-halo_j:rows+halo_j, levels+1 )

real, intent(in)  :: field_in ( 1-halo_i_data:row_length+halo_i_data           &
                              , 1-halo_j_data:rows+halo_j_data, levels+1 )


! Arguments with intent out. ie: output variables.
real, intent(out) :: field_out ( 1-halo_i_data:row_length+halo_i_data          &
                               , 1-halo_j_data:rows+halo_j_data, levels )


! Local variables.

integer :: i, j, k
integer :: ompt_start, ompt_end

real :: weight_1
real :: weight_2
real :: weight_3

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName = 'P_TO_T'

! ----------------------------------------------------------------------
! Interpolate field_in (on P grid) to T grid (field_out).
! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This procedure is now thread safe. i.e. it can be called from an
! OMP parallel region. Compute_chunk_size will determine the start and end
! iterations for each thread as if it were a static schedule.
! If called from a serial region, start and end are 1 and levels

ompt_start = 1
ompt_end   = levels
! Only call compute_chunk_size if compiling with OMP
! The procedure call is protected by the optional compile
! sentinel
!$ call compute_chunk_size(1,levels,ompt_start,ompt_end)

do k=ompt_start, ompt_end
  do j = 1-halo_j_data, rows+halo_j_data
    do i = 1-halo_i_data, row_length+halo_i_data

      weight_1 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
      weight_2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
      weight_3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
      field_out (i,j, k) =                                                     &
               weight_2/weight_1 * field_in (i,j,k+1)                          &
             + weight_3/weight_1 * field_in (i,j,k)
    end do
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine p_to_t

end module p_to_t_mod
