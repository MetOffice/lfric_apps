! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_grid_namelist_mod
! Description:
! Namelist used to define a global regular grid

! lfric modules
USE constants_mod,    ONLY: imdi, rmdi

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64, real64

IMPLICIT NONE

! Grid namelist - names matching namelist created by weight generation script
REAL(KIND=real64) :: lambda_origin_targ = rmdi ! Note this is P grid origin
REAL(KIND=real64) :: phi_origin_targ = rmdi    ! Note this is P grid origin
REAL(KIND=real64) :: phi_pole =  rmdi          ! Latitude of north pole
REAL(KIND=real64) :: lambda_pole = rmdi        ! Longitude of north polexs
REAL(KIND=real64) :: delta_lambda_targ = rmdi  ! Grid spacing x direction
REAL(KIND=real64) :: delta_phi_targ = rmdi     ! Grid spacing y direction
INTEGER(KIND=int64) :: points_lambda_targ = imdi ! Num points x direction
INTEGER(KIND=int64) :: points_phi_targ = imdi    ! Num points y direction
INTEGER(KIND=int64) :: igrid_targ = imdi ! Grid staggering
LOGICAL :: rotated = .FALSE. ! Does grid have a rotated pole?

NAMELIST /grid/ lambda_origin_targ,     phi_origin_targ,    &
                lambda_pole,            phi_pole,           &
                delta_lambda_targ,      delta_phi_targ,     &
                points_lambda_targ,     points_phi_targ,    &
                igrid_targ,             rotated


END MODULE lfricinp_grid_namelist_mod
