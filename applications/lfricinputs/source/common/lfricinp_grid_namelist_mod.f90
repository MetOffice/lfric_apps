! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_grid_namelist_mod
! Description:
! Namelist used to define a global regular grid

! lfric modules
use constants_mod,    only: imdi, rmdi

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64, real64

implicit none

! Grid namelist - names matching namelist created by weight generation script
real(kind=real64) :: lambda_origin_targ = rmdi ! Note this is P grid origin
real(kind=real64) :: phi_origin_targ = rmdi    ! Note this is P grid origin
real(kind=real64) :: phi_pole =  rmdi          ! Latitude of north pole
real(kind=real64) :: lambda_pole = rmdi        ! Longitude of north polexs
real(kind=real64) :: delta_lambda_targ = rmdi  ! Grid spacing x direction
real(kind=real64) :: delta_phi_targ = rmdi     ! Grid spacing y direction
integer(kind=int64) :: points_lambda_targ = imdi ! Num points x direction
integer(kind=int64) :: points_phi_targ = imdi    ! Num points y direction
integer(kind=int64) :: igrid_targ = imdi ! Grid staggering
logical :: rotated = .false. ! Does grid have a rotated pole?

namelist /grid/ lambda_origin_targ,     phi_origin_targ,    &
                lambda_pole,            phi_pole,           &
                delta_lambda_targ,      delta_phi_targ,     &
                points_lambda_targ,     points_phi_targ,    &
                igrid_targ,             rotated


end module lfricinp_grid_namelist_mod
