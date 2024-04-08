! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfric2um_initialise_lfric2um_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: lfric2um_initialise_lfric2um

CONTAINS

!> @brief Initialises lfric2um specific infrastructure
!> This includes reading namelists and stashmaster file,
!> and reading in gridding weights
SUBROUTINE lfric2um_initialise_lfric2um()
USE lfric2um_namelists_mod, ONLY: lfric2um_config
USE lfricinp_stashmaster_mod, ONLY: lfricinp_read_stashmaster
USE lfricinp_stash_to_lfric_map_mod, ONLY: lfricinp_init_stash_to_lfric_map
USE lfric2um_regrid_weights_mod, ONLY: lfric2um_regrid_weightsfile_ctl
IMPLICIT NONE

! Read namelists
CALL lfric2um_config%load_namelists()

! Read in STASHmaster file
CALL lfricinp_read_stashmaster(lfric2um_config%stashmaster_file)

! Read in weights files
CALL lfric2um_regrid_weightsfile_ctl()

END SUBROUTINE lfric2um_initialise_lfric2um

END MODULE lfric2um_initialise_lfric2um_mod
