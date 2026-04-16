! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfric2um_initialise_lfric2um_mod

implicit none

private

public :: lfric2um_initialise_lfric2um

contains

!> @brief Initialises lfric2um specific infrastructure
!> This includes reading namelists and stashmaster file,
!> and reading in gridding weights
subroutine lfric2um_initialise_lfric2um()

use lfricinp_stashmaster_mod,        only: lfricinp_read_stashmaster
use lfricinp_stash_to_lfric_map_mod, only: lfricinp_init_stash_to_lfric_map
use lfricinp_um_grid_mod,            only: lfricinp_set_grid_from_namelist

use lfric2um_namelists_mod,          only: lfric2um_config
use lfric2um_regrid_weights_mod,     only: lfric2um_regrid_weightsfile_ctl

implicit none

! Read namelists
call lfric2um_config%load_namelists()

! Read in STASHmaster file
call lfricinp_read_stashmaster(lfric2um_config%stashmaster_file)

! Set um_grid from lfric2um namelist
call lfricinp_set_grid_from_namelist(lfric2um_config%num_snow_layers, &
                                     lfric2um_config%num_surface_types)

! Read in weights files
call lfric2um_regrid_weightsfile_ctl()

end subroutine lfric2um_initialise_lfric2um

end module lfric2um_initialise_lfric2um_mod
