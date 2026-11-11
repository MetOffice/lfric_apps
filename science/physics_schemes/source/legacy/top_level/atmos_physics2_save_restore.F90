! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!  Routines to save convection and other variables in routine atmos_physics2
!  on 1st ENDGame cycle when l_quick_ap2=.true. for conv_diag, bl_ctl etc
!  and restore them for 2nd ENDGame cycle.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level

module atmos_physics2_save_restore_mod

implicit none

character(len=*), parameter, private :: ModuleName =                           &
                                       'ATMOS_PHYSICS2_SAVE_RESTORE_MOD'

contains

subroutine ap2_init_conv_diag( rows, row_length, ntml, ntpar, nlcl, cumulus,   &
    l_shallow, l_mid_level, delthvu, ql_ad, zhpar, dzh, qcl_inv_top, zlcl,     &
    zlcl_uv, conv_type, no_cumulus, w_max, w_copy, L_cape_opt_345)

use nlsizes_namelist_mod, only: model_levels
use cv_run_mod,  only: cape_bottom, cape_top
use bl_option_mod, only: zero
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use um_types, only: r_bl

implicit none

! Variables passed in through the argument list...
integer, intent(in) ::                                                         &
  row_length,                                                                  &
  rows

integer, intent(out) ::                                                        &
  ntml      (row_length, rows),                                                &
  ntpar     (row_length, rows),                                                &
  nlcl      (row_length, rows),                                                &
  conv_type (row_length, rows)

real(r_bl), intent(in) ::                                                      &
  w_copy        (row_length, rows, 0:model_levels)

real(r_bl), intent(out) ::                                                     &
  zhpar         (row_length, rows),                                            &
  zlcl          (row_length, rows),                                            &
  zlcl_uv       (row_length, rows),                                            &
  delthvu       (row_length, rows),                                            &
  ql_ad         (row_length, rows),                                            &
  w_max         (row_length, rows),                                            &
  dzh           (row_length, rows),                                            &
  qcl_inv_top   (row_length, rows)

logical, intent(in) ::                                                         &
  L_cape_opt_345

logical, intent(out) ::                                                        &
  cumulus      (row_length, rows),                                             &
  no_cumulus   (row_length, rows),                                             &
  l_shallow    (row_length, rows),                                             &
  l_mid_level  (row_length, rows)

! Local variables...
integer :: i, j, k  ! loop indices
character(len=*), parameter ::  RoutineName = 'AP2_INIT_CONV_DIAG'
! Dr Hook
!==============================
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise conv_diag output arrays

!$OMP PARALLEL do DEFAULT(none) private(i,j) SCHEDULE(STATIC)                  &
!$OMP SHARED(rows, row_length, ntml, ntpar, nlcl, conv_type, cumulus,          &
!$OMP no_cumulus, l_shallow, l_mid_level, delthvu, ql_ad, zhpar, dzh,          &
!$OMP qcl_inv_top, zlcl, zlcl_uv, w_max)
do j = 1, rows
  do i = 1, row_length
    ntml(i,j)       = 1
    ntpar(i,j)      = 1
    nlcl(i,j)       = 1
    cumulus(i,j)    = .false.
    no_cumulus(i,j) = .false.
    l_shallow(i,j)  = .false.
    l_mid_level(i,j)= .false.
    conv_type(i,j)  = 0
    delthvu(i,j)    = zero
    ql_ad(i,j)      = zero
    zhpar(i,j)      = zero
    dzh(i,j)        = zero
    qcl_inv_top(i,j)= zero
    zlcl(i,j)       = zero
    zlcl_uv(i,j)    = zero
    ! Initialise the w_max array.
    w_max(i,j)      = zero
  end do
end do
!$OMP end PARALLEL do

if (L_cape_opt_345) then
!$OMP PARALLEL DEFAULT(none) private(i,j,k)                                    &
!$OMP SHARED(rows, row_length, cape_bottom, cape_top,                          &
!$OMP        w_copy, w_max)
  !   Find w_max for each column. The w_max array is initialised just
  !   before the start of the OpenMP parallel region.
  do k =  cape_bottom, cape_top
!$OMP do SCHEDULE(STATIC)
    do j = 1, rows
      do i = 1, row_length
        w_max(i,j) = max(w_max(i,j), w_copy(i,j,k))
      end do
    end do
!$OMP end do NOWAIT
  end do
!$OMP end PARALLEL
end if  ! L_cape_opt_345

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine ap2_init_conv_diag
end module atmos_physics2_save_restore_mod
