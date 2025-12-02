! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level

! Purpose: STASH parameter definitions

module version_mod

implicit none

! Max. no. of STASH sections per internal model
integer,parameter :: nsectp=99
! Max. no. of STASH items per section
integer,parameter :: nitemp=999
! Max. no. of STASH list records (prognostic + diagnostic)
integer,parameter :: nrecdp=4500
! Max. no. of output times tables in STASHC
integer,parameter :: ntimep=100
! Max. no. of time profiles in STASHC
integer,parameter :: nproftp=100
! Max. no. of domain profiles/levels lists in STASHC (used for both)
integer,parameter :: nprofdp=100
! Max. total no. of time series in STASHC
integer,parameter :: NTimSerP=1500
! Max. no. time series per domain profile
integer,parameter :: tsdp=250
! Max. no. of useage profiles in STASHC
integer,parameter :: nprofup=40
! Max. no. of packages in STASHC
integer,parameter :: nprofpp=30
! Max. no. of levels in a levels list
integer,parameter :: nlevp=50
! Max. no. of pseudo levels in a  pseudo levels list
integer,parameter :: npslevp=5000
! Max. no. of pseudo levels lists in STASHC
integer,parameter :: npslistp=150
! Max. no. non-blank records in PPXREF file
integer,parameter :: ndiagp=6000
integer,parameter :: ndiagpm=nrecdp
integer,parameter :: nelemp=39
integer,parameter :: nlevp_S=nlevp*6+1
integer,parameter :: nlevlstsp=nprofdp

end module version_mod
