!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Sets up I/O configuration from within GungHo.
!>
!> @details Collects configuration information relevant for the I/O subsystem
!>          and formats it so that it can be passed to the infrastructure
!>
module gungho_setup_io_mod

  use constants_mod,             only: r_def, i_def, str_def, &
                                       str_max_filename, r_second
  use driver_modeldb_mod,        only: modeldb_type
  use file_mod,                  only: FILE_MODE_READ, &
                                       FILE_MODE_WRITE
  use lfric_xios_file_mod,       only: lfric_xios_file_type, &
                                       OPERATION_TIMESERIES, &
                                       CONVENTION_CF
  use lfric_xios_write_mod,      only: create_checkpoint_list
  use linked_list_mod,           only: linked_list_type
  use log_mod,                   only: log_event, log_level_error, &
                                       log_scratch_space, log_level_info
  use model_clock_mod,           only: model_clock_type
  ! Configuration modules
  use files_config_mod,          only: ancil_directory,           &
                                       checkpoint_stem_name,      &
                                       land_area_ancil_path,      &
                                       orography_mean_ancil_path, &
                                       orography_subgrid_ancil_path,&
                                       aerosols_ancil_path,       &
                                       albedo_nir_ancil_path,     &
                                       albedo_vis_ancil_path,     &
                                       emiss_bc_biofuel_ancil_path,&
                                       emiss_bc_fossil_ancil_path, &
                                       emiss_bc_biomass_ancil_path,&
                                       emiss_bc_biomass_hi_ancil_path,&
                                       emiss_bc_biomass_lo_ancil_path,&
                                       gas_mmr_ancil_path,         &
                                       emiss_dms_land_ancil_path,  &
                                       dms_conc_ocean_ancil_path,  &
                                       emiss_monoterp_ancil_path,  &
                                       emiss_om_biofuel_ancil_path,&
                                       emiss_om_fossil_ancil_path, &
                                       emiss_om_biomass_ancil_path,&
                                       emiss_om_biomass_hi_ancil_path,&
                                       emiss_om_biomass_lo_ancil_path,&
                                       emiss_so2_low_ancil_path,  &
                                       emiss_so2_high_ancil_path, &
                                       emiss_so2_nat_ancil_path,  &
                                       emiss_c2h6_ancil_path,     &
                                       emiss_c3h8_ancil_path,     &
                                       emiss_c5h8_ancil_path,     &
                                       emiss_ch4_ancil_path,      &
                                       emiss_co_ancil_path,       &
                                       emiss_hcho_ancil_path,     &
                                       emiss_me2co_ancil_path,    &
                                       emiss_mecho_ancil_path,    &
                                       emiss_nh3_ancil_path,      &
                                       emiss_no_ancil_path,       &
                                       emiss_meoh_ancil_path,     &
                                       emiss_no_aircrft_ancil_path, &
                                       ea_ancil_directory,        &
                                       cloud_drop_no_conc_ancil_path, &
                                       easy_asymmetry_sw_ancil_path,&
                                       easy_asymmetry_lw_ancil_path,&
                                       easy_absorption_sw_ancil_path,&
                                       easy_absorption_lw_ancil_path,&
                                       easy_extinction_sw_ancil_path,&
                                       easy_extinction_lw_ancil_path,&
                                       hydtop_ancil_path,         &
                                       h2o2_limit_ancil_path,     &
                                       ho2_ancil_path,            &
                                       no3_ancil_path,            &
                                       o3_ancil_path,             &
                                       oh_ancil_path,             &
                                       ozone_ancil_path,          &
                                       plant_func_ancil_path,     &
                                       sea_ancil_path,            &
                                       sea_ice_ancil_path,        &
                                       snow_analysis_ancil_path,  &
                                       soil_ancil_path,           &
                                       soil_dust_ancil_path,      &
                                       emiss_murk_ancil_path,     &
                                       soil_rough_ancil_path,     &
                                       sst_ancil_path,            &
                                       surface_frac_ancil_path,   &
                                       urban_ancil_path,          &
                                       start_dump_filename,       &
                                       start_dump_directory,      &
                                       iau_path,                  &
#ifdef UM_PHYSICS
                                       iau_pert_path,             &
#endif
                                       iau_sst_path,              &
                                       iau_surf_path,             &
                                       lbc_filename,              &
                                       lbc_directory,             &
                                       ls_filename,               &
                                       ls_directory,              &
                                       coarse_ancil_directory,    &
                                       internal_flux_ancil_path
  use initialization_config_mod, only: init_option,               &
                                       init_option_fd_start_dump, &
                                       init_option_checkpoint_dump,&
                                       ancil_option,              &
                                       ancil_option_start_dump,   &
                                       ancil_option_fixed,        &
                                       ancil_option_idealised,    &
                                       ancil_option_updating,     &
                                       lbc_option,                &
                                       lbc_option_gungho_file,    &
                                       lbc_option_um2lfric_file,  &
                                       ls_option,                 &
                                       ls_option_file,            &
                                       sst_source,                &
                                       sst_source_start_dump,     &
                                       coarse_aerosol_ancil,      &
                                       coarse_orography_ancil,    &
                                       coarse_ozone_ancil,        &
                                       snow_source,               &
                                       snow_source_surf
  use io_config_mod,             only: use_xios_io,               &
                                       diagnostic_frequency,      &
                                       checkpoint_write,          &
                                       checkpoint_read,           &
                                       checkpoint_times,          &
                                       end_of_run_checkpoint,     &
                                       write_dump, write_diag,    &
                                       diag_active_files,         &
                                       diag_always_on_sampling
  use orography_config_mod,      only: orog_init_option,          &
                                       orog_init_option_ancil,    &
                                       orog_init_option_start_dump
  use section_choice_config_mod, only: iau,                       &
                                       iau_sst,                   &
                                       iau_surf
  use time_config_mod,           only: timestep_start,            &
                                       timestep_end
  use derived_config_mod,        only: l_couple_sea_ice
  use lfric_string_mod,          only: split_string
#ifdef UM_PHYSICS
  use jules_surface_config_mod,  only: l_vary_z0m_soil, l_urban2t
  use specified_surface_config_mod, only: internal_flux_method, &
                                       internal_flux_method_non_uniform, &
                                       surf_temp_forcing, &
                                       surf_temp_forcing_int_flux
  use jules_radiation_config_mod, only: l_sea_alb_var_chl, l_albedo_obs
  use aerosol_config_mod,        only: glomap_mode,               &
                                       glomap_mode_climatology,   &
                                       glomap_mode_dust_and_clim, &
                                       glomap_mode_ukca,          &
                                       emissions, emissions_GC3,  &
                                       emissions_GC5,             &
                                       easyaerosol_cdnc,          &
                                       easyaerosol_sw,            &
                                       easyaerosol_lw,            &
                                       murk_prognostic
  use chemistry_config_mod,      only: chem_scheme,               &
                                       chem_scheme_strattrop,     &
                                       chem_scheme_strat_test,    &
                                       chem_scheme_offline_ox
  use iau_config_mod,            only: iau_use_pertinc
#endif

  implicit none

  private
  public :: init_gungho_files

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Adds details of all files-of-interest to a list.
  !>
  !> @param[out]   files_list Array of lfric_xios_file_type objects.
  !> @param[inout] modeldb    The modeldb holding the model state [unused here
  !>                          but required by the procedure interface]
  !>
  subroutine init_gungho_files( files_list, modeldb )

    implicit none

    type(linked_list_type),                intent(out) :: files_list
    type(modeldb_type), optional,          intent(inout)  :: modeldb

    character(len=str_max_filename) :: checkpoint_write_fname,  &
                                       checkpoint_context_name, &
                                       checkpoint_read_fname,   &
                                       dump_fname,              &
                                       ancil_fname,             &
                                       iau_fname,               &
                                       iau_surf_fname,          &
                                       lbc_fname,               &
                                       ls_fname
    character(:), allocatable :: split_checkpoint_stem_name(:)
    real(r_second), allocatable :: all_checkpoint_times(:)
#ifdef UM_PHYSICS
    character(len=str_max_filename) :: aerosol_ancil_directory
    character(len=str_max_filename) :: ozone_ancil_directory
#endif
    integer(i_def)                  :: ts_start, ts_end
    integer(i_def)                  :: rc
    integer(i_def)                  :: i
    integer(i_def)                  :: time_point

    ! Only proceed if XIOS is being used for I/O
    if (.not. use_xios_io) return

    ! Get time configuration in integer form
    read(timestep_start,*,iostat=rc)  ts_start
    read(timestep_end,*,iostat=rc)    ts_end

    ! Setup initial output file - no initial output for continuation runs
    if ( .not. checkpoint_read ) then
      call files_list%insert_item( lfric_xios_file_type( "lfric_initial",         &
                                                         xios_id="lfric_initial", &
                                                         io_mode=FILE_MODE_WRITE ) )
    end if

    ! Setup diagnostic output files
    if (write_diag) then
      do i=1, size(diag_active_files)
        call files_list%insert_item(                                                &
              lfric_xios_file_type(                                                 &
                diag_active_files(i),                                               &
                xios_id=diag_active_files(i),                                       &
                io_mode=FILE_MODE_WRITE,                                            &
                freq=diagnostic_frequency,                                          &
                is_diag=.true.,                                                     &
                diag_always_on_sampling=diag_always_on_sampling                     &
              )                                                                     &
        )
      end do
    end if

    ! Setup dump-writing context information
    if ( write_dump ) then
      ! Create dump filename from base name and end timestep
      write(dump_fname,'(A)') &
         trim(start_dump_directory)//'/'//trim(start_dump_filename),"_", &
         timestep_end
      ! In the future, Gungho/lfric_atm should define the fields that are present
      ! in the dump file here.
      call files_list%insert_item( lfric_xios_file_type( dump_fname,              &
                                                         xios_id="lfric_fd_dump", &
                                                         io_mode=FILE_MODE_WRITE, &
                                                         freq=ts_end ) )
    end if

    ! Setup dump-reading context information
    if( ((init_option == init_option_fd_start_dump .or. &
         ancil_option == ancil_option_start_dump) .and. &
         .not. checkpoint_read) .or. &
         orog_init_option == orog_init_option_start_dump) then
      ! Create dump filename from stem
      write(dump_fname,'(A)') trim(start_dump_directory)//'/'// &
                              trim(start_dump_filename)

      ! In the future, Gungho/lfric_atm should define the fields that are present
      ! in the dump file here.
      call files_list%insert_item( lfric_xios_file_type( dump_fname,                   &
                                                         xios_id="read_lfric_fd_dump", &
                                                         io_mode=FILE_MODE_READ ) )
    end if

#ifdef UM_PHYSICS
    ! Setup ancillary files
    if( ancil_option == ancil_option_fixed .or. &
        ancil_option == ancil_option_updating ) then

      ! Only read static ancils on cold start from UM
      if (init_option == init_option_fd_start_dump .and. &
           .not. checkpoint_read) then

        ! Set orography ancil filename from namelist
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(orography_subgrid_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                         xios_id="orography_subgrid_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        ! Set land area ancil filename from namelist
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(land_area_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                         xios_id="land_area_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        ! Set surface fraction ancil filename from namelist
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(surface_frac_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                         xios_id="surface_frac_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        if (l_urban2t) then
          ! Set urban ancil filename from namelist
          write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                   trim(urban_ancil_path)
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                           xios_id="urban_ancil", &
                                                           io_mode=FILE_MODE_READ ) )
        end if

        ! Set soil ancil filename from namelist
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(soil_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                         xios_id="soil_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        if (l_vary_z0m_soil) then
          ! Set soil roughness ancil filename from namelist
          write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                   trim(soil_rough_ancil_path)
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                           xios_id="soil_rough_ancil", &
                                                           io_mode=FILE_MODE_READ ) )
        end if

        ! Set topmodel hydrology filename from namelist
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(hydtop_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                         xios_id="hydtop_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        if ( ( glomap_mode == glomap_mode_dust_and_clim ).or. &
             ( glomap_mode == glomap_mode_ukca   ) ) then

          ! Set aerosol emission ancil filenames from namelist
          write(ancil_fname,'(A)') trim(ancil_directory)//'/'//                &
                                   trim(soil_dust_ancil_path)
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="soil_dust_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
        end if

      end if ! static ancils on cold start

      ! Only read updating ancils on new run, if updating or using surf
      if (ancil_option == ancil_option_updating .or. &
           .not. checkpoint_read) then
        ! need to add .or. use_surf_analysis here

        if (snow_source == snow_source_surf ) then
          if (snow_analysis_ancil_path(1:1) == '/') then
            write(ancil_fname,'(A)') trim(snow_analysis_ancil_path)
          else
            write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                     trim(snow_analysis_ancil_path)
          end if
          call files_list%insert_item( lfric_xios_file_type( ancil_fname, &
                                          xios_id="snow_analysis_ancil",  &
                                          io_mode=FILE_MODE_READ ) )
        end if

        ! Set sea surface temperature ancil filename from namelist
        if (sst_source /= sst_source_start_dump) then
          if (sst_ancil_path(1:1) == '/') then
            write(ancil_fname,'(A)') trim(sst_ancil_path)
          else
            write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                    trim(sst_ancil_path)
          end if
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="sst_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
        end if

        ! Set sea ice ancil filename from namelist
        if (.not. l_couple_sea_ice) then
          if (sea_ice_ancil_path(1:1) == '/') then
            write(ancil_fname,'(A)') trim(sea_ice_ancil_path)
          else
            write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                    trim(sea_ice_ancil_path)
          end if
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="sea_ice_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
        end if

      end if

      ! Only read updating ancils on new run or if updating
      if (ancil_option == ancil_option_updating .or. &
           .not. checkpoint_read) then

        ! Set plant functional type ancil filename from namelist
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(plant_func_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                         xios_id="plant_func_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        ! Set sea chlorophyll ancil filename from namelist
        if ( l_sea_alb_var_chl ) then
          write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                   trim(sea_ancil_path)
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                           xios_id="sea_ancil", &
                                                           io_mode=FILE_MODE_READ ) )
        end if

        if ( l_albedo_obs ) then
          ! Set albedo_vis ancil filename from namelist
          write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                   trim(albedo_vis_ancil_path)
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="albedo_vis_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

          ! Set albedo_nir ancil filename from namelist
          write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                   trim(albedo_nir_ancil_path)
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="albedo_nir_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
        end if

        ! Set ozone filename from namelist
        if ( coarse_ozone_ancil ) then
          ozone_ancil_directory = coarse_ancil_directory
        else
          ozone_ancil_directory = ancil_directory
        end if
        write(ancil_fname,'(A)') trim(ozone_ancil_directory)//'/'// &
                                 trim(ozone_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                         xios_id="ozone_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        if (murk_prognostic) then
          ! Set aerosol emission ancil filenames from namelist
          if (emiss_murk_ancil_path(1:1) == '/') then
            write(ancil_fname,'(A)') trim(emiss_murk_ancil_path)
          else
            write(ancil_fname,'(A)') trim(ancil_directory)//'/'//              &
                                     trim(emiss_murk_ancil_path)
          end if
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="emiss_murk_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
        end if

        if ( ( glomap_mode == glomap_mode_dust_and_clim ) .or.    &
             ( glomap_mode == glomap_mode_climatology   ) ) then

          ! Set aerosol ancil filename from namelist
          if ( coarse_aerosol_ancil ) then
            aerosol_ancil_directory = coarse_ancil_directory
          else
            aerosol_ancil_directory = ancil_directory
          end if
          write(ancil_fname,'(A)') trim(aerosol_ancil_directory)//'/'// &
                                   trim(aerosols_ancil_path)
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="aerosols_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
        end if

      end if ! updating or a new run

    end if ! fixed or updating ancils

    ! Chemistry ancils
    if ( (chem_scheme == chem_scheme_strattrop    .or.        &
          chem_scheme == chem_scheme_strat_test)   .and.     &
         ancil_option == ancil_option_updating ) then

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_c2h6_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                &
                                                         xios_id="emiss_c2h6_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_c3h8_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                &
                                                         xios_id="emiss_c3h8_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_c5h8_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                &
                                                         xios_id="emiss_c5h8_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_ch4_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,               &
                                                         xios_id="emiss_ch4_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_co_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,              &
                                                         xios_id="emiss_co_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_hcho_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                &
                                                         xios_id="emiss_hcho_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_me2co_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                 &
                                                         xios_id="emiss_me2co_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_mecho_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                 &
                                                         xios_id="emiss_mecho_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_nh3_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,               &
                                                         xios_id="emiss_nh3_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_no_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,              &
                                                         xios_id="emiss_no_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_meoh_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                &
                                                         xios_id="emiss_meoh_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_no_aircrft_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                      &
                                                         xios_id="emiss_no_aircrft_ancil", &
                                                         io_mode=FILE_MODE_READ ) )


    endif !chemistry ancils

    if ( glomap_mode == glomap_mode_ukca   .and.            &
         ancil_option == ancil_option_updating ) then

      ! Set aerosol emission ancil filenames from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_bc_biofuel_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                      &
                                                         xios_id="emiss_bc_biofuel_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_bc_fossil_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                     &
                                                         xios_id="emiss_bc_fossil_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      if (emissions == emissions_GC3) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(emiss_bc_biomass_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,                      &
                                                           xios_id="emiss_bc_biomass_ancil", &
                                                           io_mode=FILE_MODE_READ ) )
      else if (emissions == emissions_GC5) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(emiss_bc_biomass_hi_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,                         &
                                                           xios_id="emiss_bc_biomass_hi_ancil", &
                                                           io_mode=FILE_MODE_READ ) )

        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(emiss_bc_biomass_lo_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,                         &
                                                           xios_id="emiss_bc_biomass_lo_ancil", &
                                                           io_mode=FILE_MODE_READ ) )
      end if

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_dms_land_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                    &
                                                         xios_id="emiss_dms_land_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(dms_conc_ocean_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                    &
                                                         xios_id="dms_conc_ocean_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_monoterp_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                    &
                                                         xios_id="emiss_monoterp_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_om_biofuel_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                      &
                                                         xios_id="emiss_om_biofuel_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_om_fossil_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                     &
                                                         xios_id="emiss_om_fossil_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      if (emissions == emissions_GC3) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(emiss_om_biomass_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,                      &
                                                           xios_id="emiss_om_biomass_ancil", &
                                                           io_mode=FILE_MODE_READ ) )
      else if (emissions == emissions_GC5) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(emiss_om_biomass_hi_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,                         &
                                                           xios_id="emiss_om_biomass_hi_ancil", &
                                                           io_mode=FILE_MODE_READ ) )

        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(emiss_om_biomass_lo_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,                         &
                                                           xios_id="emiss_om_biomass_lo_ancil", &
                                                           io_mode=FILE_MODE_READ ) )
      end if

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_so2_low_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,                   &
                                                         xios_id="emiss_so2_low_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

      if (emissions == emissions_GC3) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(emiss_so2_high_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,                    &
                                                           xios_id="emiss_so2_high_ancil", &
                                                           io_mode=FILE_MODE_READ ) )
      end if

      ! single time file only needs reading on cold start
      if (init_option == init_option_fd_start_dump .and. &
           .not. checkpoint_read) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(emiss_so2_nat_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,        &
                                                         xios_id="emiss_so2_nat_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
      end if

      ! Setup Offline oxidants ancillary files
      if ( chem_scheme == chem_scheme_offline_ox ) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(h2o2_limit_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,              &
                                                         xios_id="h2o2_limit_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(ho2_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,       &
                                                         xios_id="ho2_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(no3_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,       &
                                                         xios_id="no3_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(o3_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="o3_ancil", &
                                                         io_mode=FILE_MODE_READ ) )

        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(oh_ancil_path)
        call files_list%insert_item( lfric_xios_file_type( ancil_fname,      &
                                                         xios_id="oh_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
      end if  ! Offline oxidants scheme

    end if !ukca ancils


   ! Easy Aerosols
   if ( easyaerosol_cdnc .and. ancil_option == ancil_option_updating )  then
      write(ancil_fname,'(A)') trim(ea_ancil_directory)//'/'// &
                               trim(cloud_drop_no_conc_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,           &
           xios_id="cloud_drop_no_conc_ancil", io_mode=FILE_MODE_READ ) )
   endif
   if ( easyaerosol_sw .and. ancil_option == ancil_option_updating )  then
      write(ancil_fname,'(A)') trim(ea_ancil_directory)//'/'// &
                               trim(easy_asymmetry_sw_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,           &
           xios_id="easy_asymmetry_sw_ancil", io_mode=FILE_MODE_READ ) )
      write(ancil_fname,'(A)') trim(ea_ancil_directory)//'/'// &
                               trim(easy_absorption_sw_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,           &
           xios_id="easy_absorption_sw_ancil", io_mode=FILE_MODE_READ ) )
      write(ancil_fname,'(A)') trim(ea_ancil_directory)//'/'// &
                               trim(easy_extinction_sw_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,           &
           xios_id="easy_extinction_sw_ancil", io_mode=FILE_MODE_READ ) )
   endif
   if ( easyaerosol_lw .and. ancil_option == ancil_option_updating )  then
      write(ancil_fname,'(A)') trim(ea_ancil_directory)//'/'// &
                               trim(easy_asymmetry_lw_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,           &
           xios_id="easy_asymmetry_lw_ancil", io_mode=FILE_MODE_READ ) )
      write(ancil_fname,'(A)') trim(ea_ancil_directory)//'/'// &
                               trim(easy_absorption_lw_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,           &
           xios_id="easy_absorption_lw_ancil", io_mode=FILE_MODE_READ ) )
      write(ancil_fname,'(A)') trim(ea_ancil_directory)//'/'// &
                               trim(easy_extinction_lw_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,           &
           xios_id="easy_extinction_lw_ancil", io_mode=FILE_MODE_READ ) )
   endif

    ! Radiatively active gases
    if ( ancil_option == ancil_option_idealised ) then
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(gas_mmr_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname, &
                                                         xios_id="gas_mmr_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
    end if

    ! Prescribed forcing of surface temperature e.g. via internal flux
    if ( surf_temp_forcing == surf_temp_forcing_int_flux .and. &
         internal_flux_method == internal_flux_method_non_uniform ) then
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(internal_flux_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname, &
                                                         xios_id="int_flux_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
    end if
#endif

    ! Setup orography ancillary file
    if ( orog_init_option == orog_init_option_ancil ) then
        if ( coarse_orography_ancil ) then
          ! Set orography ancil filename from namelist
          write(ancil_fname,'(A)') trim(coarse_ancil_directory)//'/'// &
          trim(orography_mean_ancil_path)
          call files_list%insert_item( lfric_xios_file_type( ancil_fname,   &
                                          xios_id="coarse_orography_ancil", &
                                          io_mode=FILE_MODE_READ ) )
        end if
      ! Set orography ancil filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(orography_mean_ancil_path)
      call files_list%insert_item( lfric_xios_file_type( ancil_fname,               &
                                                         xios_id="orography_mean_ancil", &
                                                         io_mode=FILE_MODE_READ ) )
    end if

    ! Setup the IAU file
    if ( iau ) then
      write(iau_fname,'(A)') trim(iau_path)
      call files_list%insert_item( lfric_xios_file_type( iau_fname,     &
                                                         xios_id="iau", &
                                                         io_mode=FILE_MODE_READ ))
#ifdef UM_PHYSICS
      ! Setup the IAU pert increments file
      if ( iau_use_pertinc ) then
        call files_list%insert_item( lfric_xios_file_type( trim(iau_pert_path), &
                                                           xios_id="iau_pert",  &
                                                           io_mode=FILE_MODE_READ ))
      end if
#endif
    end if

    ! Setup the IAU surface inc file
    if ( iau_surf )  then
      write(iau_surf_fname,'(A)') trim(iau_surf_path)
      call files_list%insert_item( lfric_xios_file_type( iau_surf_fname,     &
                                                         xios_id="iau_surf", &
                                                         io_mode=FILE_MODE_READ ))
    end if

   ! Setup the IAU SST inc file
    if ( iau_sst )  then
      call files_list%insert_item( lfric_xios_file_type( trim(iau_sst_path), &
                                                         xios_id="iau_sst",  &
                                                         io_mode=FILE_MODE_READ ))
    end if

    ! Setup the lbc file
    if ( lbc_option == lbc_option_gungho_file .or. &
         lbc_option == lbc_option_um2lfric_file ) then
      write(lbc_fname,'(A)') trim(lbc_directory)//'/'// &
                             trim(lbc_filename)
      call files_list%insert_item( lfric_xios_file_type( lbc_fname,     &
                                                         xios_id="lbc", &
                                                         io_mode=FILE_MODE_READ ) )
    endif

    ! Setup the ls file
    if ( ls_option == ls_option_file ) then
      write(ls_fname,'(A)') trim(ls_directory)//'/'// &
                            trim(ls_filename)
      call files_list%insert_item( lfric_xios_file_type( ls_fname,     &
                                                         xios_id="ls", &
                                                         io_mode=FILE_MODE_READ ) )
    endif

    ! Setup checkpoint writing context information
    if (checkpoint_write) then
      call create_checkpoint_list( modeldb%clock, checkpoint_times, &
                                   end_of_run_checkpoint, all_checkpoint_times)
      if ( size(all_checkpoint_times) > 0 ) then
        do i = 1, size(all_checkpoint_times)

          if( all_checkpoint_times(i) < 0.0_r_def ) then
            write(log_scratch_space, '(A,F10.3)') &
              "Negative checkpoint time not allowed: ", all_checkpoint_times(i)
            call log_event( log_scratch_space, log_level_error )
          end if

          ! Convert checkpoint time in seconds to timestep number
          time_point = modeldb%clock%steps_from_seconds( all_checkpoint_times(i) )

          if( log10(real(time_point)) >= 10.0_r_def )then
            write(log_scratch_space, '(A,I0)') &
              "Timestep number too big to fit in checkpoint filename: ", time_point
            call log_event( log_scratch_space, log_level_error )
          end if

          ! Create checkpoint filename from stem and timepoint
          write(checkpoint_write_fname,'(A,A,I10.10)') &
                               trim(checkpoint_stem_name),"_", time_point

          split_checkpoint_stem_name = split_string( trim(checkpoint_stem_name), '/' )
          ! Check whether checkpoint stem name has a trailing slash (i.e. is a directory)
          if (split_checkpoint_stem_name(size(split_checkpoint_stem_name)) == "") then
            call log_event( &
              "Could not create checkpoint context name from filename: "// &
              trim(checkpoint_stem_name), &
              log_level_error )
          end if
          write(checkpoint_context_name, '(A,A,I10.10)') &
               trim(split_checkpoint_stem_name(size(split_checkpoint_stem_name))), &
               "_", time_point
          call log_event( "Adding checkpoint file for writing: "// &
                               trim(checkpoint_write_fname), log_level_info )

          call files_list%insert_item( lfric_xios_file_type( trim(checkpoint_write_fname),                &
                                                             xios_id=trim(checkpoint_context_name), &
                                                             io_mode=FILE_MODE_WRITE,               &
                                                             freq=time_point,                       &
                                                             field_group_id="checkpoint_fields",    &
                                                             file_convention=CONVENTION_CF) )

        end do
      else
        call log_event("Checkpointing has been turned on but no checkpoint " // &
                       "times have been specified.", log_level_error)
      end if
    end if

    ! Setup checkpoint reading context information
    if ( checkpoint_read ) then
      ! Create checkpoint filename from stem and (start - 1) timestep
      if( log10(real(ts_start)) >= 10.0_r_def )then
        call log_event( &
          "Number of timesteps too big to fit in checkpoint filename", &
          log_level_error )
      end if
      write(checkpoint_read_fname,'(A,A,I10.10)') &
                   trim(checkpoint_stem_name),"_", (ts_start - 1)
      call files_list%insert_item( lfric_xios_file_type( checkpoint_read_fname,           &
                                                         xios_id="lfric_checkpoint_read", &
                                                         io_mode=FILE_MODE_READ,          &
                                                         freq=ts_start - 1,               &
                                                         field_group_id="checkpoint_fields",    &
                                                         file_convention=CONVENTION_CF ) )
    end if

    ! Read checkpoint file as though it were start dump
    if ( init_option == init_option_checkpoint_dump .and. &
         .not. checkpoint_read ) then
      ! Create dump filename from stem
      write(dump_fname,'(A)') trim(start_dump_directory)//'/'// &
                              trim(start_dump_filename)
      call files_list%insert_item( lfric_xios_file_type( dump_fname,           &
                                                         xios_id="lfric_checkpoint_read", &
                                                         io_mode=FILE_MODE_READ,          &
                                                         field_group_id="checkpoint_fields",    &
                                                         file_convention=CONVENTION_CF ) )
    end if

  end subroutine init_gungho_files

end module gungho_setup_io_mod
