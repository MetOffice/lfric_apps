! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the average w over a layer
!
module mean_w_layer_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName='MEAN_W_LAYER_MOD'

contains

! Subroutine Interface:
subroutine mean_w_layer(nunstable,row_length,rows,model_levels,                &
               k_start, index_i, index_j,                                      &
               depth, z_full_c, z_half_c, w, dmass_theta,                      &
               w_avg)

! ------------------------------------------------------------------------------
! Description:
!   Calculate the average w over a layer
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------
use parkind1, only: jpim, jprb       !DrHook
use yomhook,  only: lhook, dr_hook   !DrHook
use bl_option_mod, only: zero

implicit none

! Subroutine arguments
integer, intent(in) ::                                                         &
  nunstable            & ! Number of parcel ascents
 ,row_length           & ! Local number of points on a row
 ,rows                 & ! Local number of rows in a theta field
 ,model_levels           ! Number of model levels


integer, intent(in) ::                                                         &
  k_start(nunstable)   & ! Level above which require average
 ,index_i(nunstable)   & ! column number of unstable points
 ,index_j(nunstable)     ! row number of unstable points


real(kind=r_bl), intent(in)    ::                                              &
  depth                  ! depth of layer (m)

real(kind=r_bl), intent(in)    ::                                              &
  z_full_c(nunstable, model_levels)      & ! Height theta lev (m)
 ,z_half_c(nunstable, model_levels)      & ! Height uv lev    (m)
 ,w(row_length,rows,0:model_levels)      & ! vertical velocity (m/s)
 ,dmass_theta(nunstable,model_levels)      ! r**2rho*dr on theta levels

real(kind=r_bl), intent(out)   ::                                              &
  w_avg(nunstable)       ! Average value of w in required layer (m/s)


! Local variables

integer :: ii, k,i, j            !Loop counters

real(kind=r_bl) ::                                                             &
  mass                   ! mass of required layer

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MEAN_W_LAYER'

!-------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------

!$OMP  PARALLEL do DEFAULT(none) private(ii,k,i,j,mass)                        &
!$OMP  SHARED(nunstable,w_avg,index_i,index_j,model_levels,k_start,z_full_c,   &
!$OMP         z_half_c,depth,dmass_theta,w)                                    &
!$OMP  SCHEDULE(STATIC)
! The first index of the loop is "ii" so that it can be
! parallelised using OpenMP.
do ii=1,nunstable

  w_avg(ii) = zero
  mass      = zero

  i = index_i(ii)
  j = index_j(ii)
  do k=1,model_levels-1
    if (k >= k_start(ii) .and.                                                 &
        z_full_c(ii,k) <= (z_half_c(ii,k_start(ii)+1)+depth)) then

      mass  = mass + dmass_theta(ii,k)
      w_avg(ii) = w_avg(ii) + w(i,j,k)*dmass_theta(ii,k)

    end if
  end do

  if (mass  >  zero ) then
    w_avg(ii) = w_avg(ii)/mass
  end if
end do
!$OMP end PARALLEL do

!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine mean_w_layer

end module mean_w_layer_mod
