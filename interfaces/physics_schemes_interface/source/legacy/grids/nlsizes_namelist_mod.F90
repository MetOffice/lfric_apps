! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Sizes for many of the UM's main, dynamic data arrays
!
module nlsizes_namelist_mod

! Description:
!  Module-based interface to the nlsizes namelist and associated declarations.
!  Contains the sizes needed for the dynamic allocation of the main data arrays
!  within the model.
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v10.3 programming standards.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Grids
!

use missing_data_mod,   only: imdi

use jules_soil_mod,         only: sm_levels

implicit none

private :: imdi

character (len=*), private, parameter :: ModuleName  = "NLSIZES_NAMELIST_MOD"

save

! Main sizes of fields for each submodel
integer :: global_row_length = imdi ! Points per global row
integer :: global_rows       = imdi ! Number of global (theta) rows
integer :: model_levels      = imdi ! Number of model levels
integer :: bl_levels         = imdi ! Number of boundary-layer-levels

integer :: row_length           ! Number of points per local row
integer :: rows                 ! Number of local (theta) rows

integer :: ntiles               ! Number of land surface tiles

integer :: n_cca_lev            ! Number of CCA levels

end module nlsizes_namelist_mod
