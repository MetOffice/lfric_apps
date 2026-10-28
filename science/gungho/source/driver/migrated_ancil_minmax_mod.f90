!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Temporary development aid: logs the min/max of the ancillary fields
!>        that have been migrated to the file-owned time-varying read.
!> @details Self-contained, so that the same module can be added unchanged to a
!>          build without the migration and the logs compared. Fields are found
!>          by name in the depository, which holds them on both the legacy and
!>          the file-owned paths. Each line has the form
!>          "Min/max migrated_ancil init <field> = <min> <max>" or
!>          "Min/max migrated_ancil step <n> <field> = <min> <max>".
module migrated_ancil_minmax_mod

  use constants_mod,           only : i_def, str_def
  use driver_modeldb_mod,      only : modeldb_type
  use field_collection_mod,    only : field_collection_type
  use field_mod,               only : field_type
  use log_mod,                 only : LOG_LEVEL_INFO
  use sci_field_minmax_alg_mod, only : log_field_minmax
  use aerosol_config_mod,      only : glomap_mode,               &
                                      glomap_mode_climatology,   &
                                      glomap_mode_dust_and_clim, &
                                      easyaerosol_cdnc,          &
                                      easyaerosol_sw,            &
                                      easyaerosol_lw
  use jules_radiation_config_mod, only : l_albedo_obs, l_sea_alb_var_chl

  implicit none

  private
  public :: log_migrated_ancil_minmax

  !> Log every this many timesteps (the initial log is always written)
  integer(i_def), parameter :: log_every_n_steps = 1

  character(len=str_def), parameter :: aerosol_names(17) =                  &
    [ character(len=str_def) ::                                             &
      "acc_sol_bc", "acc_sol_om", "acc_sol_su", "acc_sol_ss", "n_acc_sol",  &
      "ait_sol_bc", "ait_sol_om", "ait_sol_su", "n_ait_sol",                &
      "ait_ins_bc", "ait_ins_om", "n_ait_ins",                              &
      "cor_sol_bc", "cor_sol_om", "cor_sol_su", "cor_sol_ss", "n_cor_sol" ]

  character(len=str_def), parameter :: albedo_names(2) =                    &
    [ character(len=str_def) :: "albedo_obs_vis", "albedo_obs_nir" ]

  character(len=str_def), parameter :: cdnc_names(1) =                      &
    [ character(len=str_def) :: "cloud_drop_no_conc" ]

  character(len=str_def), parameter :: easy_sw_names(3) =                   &
    [ character(len=str_def) ::                                             &
      "easy_asymmetry_sw", "easy_absorption_sw", "easy_extinction_sw" ]

  character(len=str_def), parameter :: easy_lw_names(3) =                   &
    [ character(len=str_def) ::                                             &
      "easy_asymmetry_lw", "easy_absorption_lw", "easy_extinction_lw" ]

  character(len=str_def), parameter :: chloro_names(1) =                    &
    [ character(len=str_def) :: "chloro_sea" ]

  ! Always logged. The plant functional type fields are included whether or
  ! not their migration is switched on.
  character(len=str_def), parameter :: always_names(3) =                    &
    [ character(len=str_def) :: "ozone", "canopy_height", "leaf_area_index" ]

contains

  !> @brief Logs the min/max of each migrated ancillary field that is in use
  !> @param[in,out] modeldb  The model database holding the depository
  !> @param[in]     at_init  True for the initial log, false for a timestep log
  subroutine log_migrated_ancil_minmax( modeldb, at_init )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb
    logical,            intent(in)    :: at_init

    type(field_collection_type), pointer :: depository
    character(len=str_def)               :: prefix

    if ( at_init ) then
      prefix = 'migrated_ancil init '
    else
      if ( mod( modeldb%clock%get_step(), log_every_n_steps ) /= 0 ) return
      write( prefix, '(A,I0,A)' ) 'migrated_ancil step ', &
                                  modeldb%clock%get_step(), ' '
    end if

    depository => modeldb%fields%get_field_collection( "depository" )

    if ( glomap_mode == glomap_mode_climatology .or. &
         glomap_mode == glomap_mode_dust_and_clim ) then
      call log_fields( depository, prefix, aerosol_names )
    end if
    if ( l_albedo_obs )      call log_fields( depository, prefix, albedo_names )
    if ( easyaerosol_cdnc )  call log_fields( depository, prefix, cdnc_names )
    if ( easyaerosol_sw )    call log_fields( depository, prefix, easy_sw_names )
    if ( easyaerosol_lw )    call log_fields( depository, prefix, easy_lw_names )
    if ( l_sea_alb_var_chl ) call log_fields( depository, prefix, chloro_names )
    call log_fields( depository, prefix, always_names )

    nullify( depository )

  end subroutine log_migrated_ancil_minmax

  !> @brief Logs the min/max of the named fields that exist in the depository
  subroutine log_fields( depository, prefix, names )

    implicit none

    type(field_collection_type), intent(inout) :: depository
    character(len=*),            intent(in)    :: prefix
    character(len=*),            intent(in)    :: names(:)

    type(field_type), pointer :: fld => null()
    integer(i_def)            :: i

    do i = 1, size( names )
      if ( depository%field_exists( trim(names(i)) ) ) then
        call depository%get_field( trim(names(i)), fld )
        call log_field_minmax( LOG_LEVEL_INFO, &
                               trim(prefix)//trim(names(i)), fld )
      end if
    end do

    nullify( fld )

  end subroutine log_fields

end module migrated_ancil_minmax_mod
