! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To solve the tri-diagonal equations for the prognostic
!           variables in the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_implic_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_IMPLIC_MOD'
contains

subroutine mym_implic(levels, kst, ken, aa, bb, cc, qq)

use atm_fields_bounds_mod, only: pdims
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

integer, intent(in) ::                                                         &
   levels,                                                                     &
                  ! number of levels of variables to be solved
   kst,                                                                        &
                  ! index of start level to be solved
   ken
                  ! index of emd level to be solved

real(kind=r_bl), intent(in out) ::                                      &
   aa(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, levels),            &
                  ! coefficients of fields on level K-1
                  ! in the tri-diagonal equation
   bb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, levels),            &
                  ! coefficients on fields level K
                  ! in the tri-diagonal equation
   cc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, levels),            &
                  ! coefficients on fields level K+1
                  ! in the tri-diagonal equation
   qq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, levels)
                  ! right hand side of the tri-diagonal equation

! Local variables
integer ::                                                                     &
   i, j, k
                  ! Loop indexes

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_IMPLIC'

! Solve from top to bottom
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = ken, kst + 1, -1
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ! normalising so that bb = 1.0
      aa(i, j, k) = aa(i, j, k) / bb(i, j, k)
      qq(i, j, k) = qq(i, j, k) / bb(i, j, k)

      bb(i, j, k - 1) = bb(i, j, k - 1) - cc(i, j, k - 1) * aa(i, j, k)
      qq(i, j, k - 1) = qq(i, j, k - 1) - cc(i, j, k - 1) * qq(i, j, k)
    end do
  end do
end do

do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    qq(i, j, kst) = qq(i, j, kst) / bb(i, j, kst)
  end do
end do

! Solve from bottom to top
do k = kst + 1, ken
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      qq(i, j, k) = qq(i, j, k) - aa(i, j, k) * qq(i, j, k - 1)
    end do
  end do
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mym_implic
end module mym_implic_mod
