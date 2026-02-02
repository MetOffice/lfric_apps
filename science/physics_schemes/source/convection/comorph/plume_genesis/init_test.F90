! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module init_test_mod

implicit none

contains


! Subroutine to calculate the mask of points where convection
! initiation mass sources might occur on a given level.
! Later calculations will further narrow down the subset
! of points where convective initiation actually occurs.

! Note: this calculation is done on the full fields, therefore
! it is done at the precision of the host-model (real_hmprec),
! not that of the rest of the convection scheme (real_cvprec)

subroutine init_test( grid, fields, virt_temp, l_init_poss )

use comorph_constants_mod, only: real_hmprec,                                  &
                     nx_full, ny_full, k_bot_conv, k_top_conv, k_top_init,     &
                     l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup
use grid_type_mod, only: grid_type
use fields_type_mod, only: fields_type
use test_unstable_mod, only: test_unstable


implicit none

! Structure containing pointers to height and pressure
type(grid_type), intent(in) :: grid

! Structure containing pointers to primary fields
type(fields_type), intent(in) :: fields

! Virtual temperature profile
real(kind=real_hmprec), intent(in) :: virt_temp                                &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Output mask for points where convective mass sources
! might initiate
logical, intent(out) :: l_init_poss                                            &
                     ( nx_full, ny_full, k_bot_conv:k_top_init )

! Lower and upper bounds on the height and pressure arrays
integer :: lb_z(3), ub_z(3)
integer :: lb_p(3), ub_p(3)

! Lower and upper bounds of mixing ratio arrays
integer :: lb_v(3), ub_v(3)
integer :: lb_l(3), ub_l(3)
integer :: lb_r(3), ub_r(3)
integer :: lb_f(3), ub_f(3)
integer :: lb_s(3), ub_s(3)
integer :: lb_g(3), ub_g(3)

! Loop counters
integer :: i, j, k


! Find bounds of required 3-D arrays
lb_z = lbound(grid % height_full)
ub_z = ubound(grid % height_full)
lb_p = lbound(grid % pressure_full)
ub_p = ubound(grid % pressure_full)
lb_v = lbound(fields % q_vap)
ub_v = ubound(fields % q_vap)
lb_l = lbound(fields % q_cl)
ub_l = ubound(fields % q_cl)
lb_r = lbound(fields % q_rain)
ub_r = ubound(fields % q_rain)
lb_f = lbound(fields % q_cf)
ub_f = ubound(fields % q_cf)
lb_s = lbound(fields % q_snow)
ub_s = ubound(fields % q_snow)
lb_g = lbound(fields % q_graup)
ub_g = ubound(fields % q_graup)

! Loop over levels
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                               &
!$OMP SHARED( nx_full, ny_full, k_bot_conv, k_top_conv, k_top_init,            &
!$OMP         grid, fields, virt_temp, l_init_poss, l_cv_snow,                 &
!$OMP         lb_z, ub_z, lb_p, ub_p, lb_v, ub_v, lb_l, ub_l, lb_r, ub_r,      &
!$OMP         lb_f, ub_f, lb_s, ub_s, lb_g, ub_g )                             &
!$OMP PRIVATE( i, j, k )
do k = k_bot_conv, k_top_init

  ! Initialise mask of points to false
  do j = 1, ny_full
    do i = 1, nx_full
      l_init_poss(i,j,k) = .false.
    end do
  end do

  ! Flag points as possible initiation sources if they are statically unstable,
  ! based on a dry-adiabat up or down vs the start-of-timestep Tv profile
  ! (calculations done inside test_unstable)...

  ! If not at the model-top, try lifting up one level
  if ( k < k_top_conv ) then

    call test_unstable( virt_temp(:,:,k), virt_temp(:,:,k+1),                  &
                        lb_z(1:2), ub_z(1:2), grid % height_full(:,:,k),       &
                                              grid % height_full(:,:,k+1),     &
                        lb_p(1:2), ub_p(1:2), grid % pressure_full(:,:,k),     &
                                              grid % pressure_full(:,:,k+1),   &
                        lb_v(1:2), ub_v(1:2), fields % q_vap(:,:,k),           &
                        lb_l(1:2), ub_l(1:2), fields % q_cl(:,:,k),            &
                        lb_r(1:2), ub_r(1:2), fields % q_rain(:,:,k),          &
                        lb_f(1:2), ub_f(1:2), fields % q_cf(:,:,k),            &
                        lb_s(1:2), ub_s(1:2), fields % q_snow(:,:,k),          &
                        lb_g(1:2), ub_g(1:2), fields % q_graup(:,:,k),         &
                        l_init_poss(:,:,k) )

  end if  ! ( k < k_top_conv )

  ! If not at the model-bottom, try subsiding down one level
  if ( k > k_bot_conv ) then

    call test_unstable( virt_temp(:,:,k), virt_temp(:,:,k-1),                  &
                        lb_z(1:2), ub_z(1:2), grid % height_full(:,:,k),       &
                                              grid % height_full(:,:,k-1),     &
                        lb_p(1:2), ub_p(1:2), grid % pressure_full(:,:,k),     &
                                              grid % pressure_full(:,:,k-1),   &
                        lb_v(1:2), ub_v(1:2), fields % q_vap(:,:,k),           &
                        lb_l(1:2), ub_l(1:2), fields % q_cl(:,:,k),            &
                        lb_r(1:2), ub_r(1:2), fields % q_rain(:,:,k),          &
                        lb_f(1:2), ub_f(1:2), fields % q_cf(:,:,k),            &
                        lb_s(1:2), ub_s(1:2), fields % q_snow(:,:,k),          &
                        lb_g(1:2), ub_g(1:2), fields % q_graup(:,:,k),         &
                        l_init_poss(:,:,k) )

  end if  ! ( k > k_bot_conv )

  ! Flag points as potential initiation sources if they
  ! contain condensed water of any sort
  call moist_test( lb_l(1:2), ub_l(1:2), fields % q_cl(:,:,k),                 &
                   l_init_poss(:,:,k) )
  if ( l_cv_rain ) then
    call moist_test( lb_r(1:2), ub_r(1:2), fields % q_rain(:,:,k),             &
                     l_init_poss(:,:,k) )
  end if
  if ( l_cv_cf ) then
    call moist_test( lb_f(1:2), ub_f(1:2), fields % q_cf(:,:,k),               &
                     l_init_poss(:,:,k) )
  end if
  if ( l_cv_snow ) then
    call moist_test( lb_s(1:2), ub_s(1:2), fields % q_snow(:,:,k),             &
                     l_init_poss(:,:,k) )
  end if
  if ( l_cv_graup ) then
    call moist_test( lb_g(1:2), ub_g(1:2), fields % q_graup(:,:,k),            &
                     l_init_poss(:,:,k) )
  end if

end do  ! k = k_bot_conv, k_top_init
!$OMP END PARALLEL DO


return
end subroutine init_test


! Subroutine to test whether a condensed water field is nonzero,
! and set l_init to true at points where it is.
! The main point of having a separate subroutine for this trivial
! task is so that the condensed water array slice qc (which is
! a pointer component of the fields structure in init_test above)
! can be declared as a normal array, allowing faster
! vectorisation.
subroutine moist_test( lb, ub, qc, l_init_poss )

use comorph_constants_mod, only: real_hmprec, nx_full, ny_full

implicit none

! Array lower and upper bounds
integer, intent(in) :: lb(2), ub(2)
! 2-D array slice of condensed water
real(kind=real_hmprec), intent(in) :: qc( lb(1):ub(1), lb(2):ub(2) )
! Flag for whether each point contains any condensed water
logical, intent(in out) :: l_init_poss( nx_full, ny_full )
real(kind=real_hmprec), parameter :: zero = 0.0_real_hmprec
! Loop counters
integer :: i, j

do j = 1, ny_full
  do i = 1, nx_full
    l_init_poss(i,j) = l_init_poss(i,j) .or. qc(i,j) > zero
  end do
end do

return
end subroutine moist_test


end module init_test_mod
