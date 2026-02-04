! ADD Licence?????
!================================================================
! Temporary code to check new and old confiuration objects return
! the same configuration values. This is to manage the transition
! of the codebase from using a namelist_collection_type to a
! config_type
!================================================================
module check_config_api_mod

  use constants_mod,           only: l_def, i_def, r_def,       &
                                     str_def, str_max_filename, &
                                     r_second, i_medium
  use config_mod,              only: config_type
  use log_mod,                 only: log_event, log_level_warning, &
                                     log_level_info
  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type

  implicit none

  private
  public :: check_config_api

contains

subroutine check_config_api( configuration, config )

  implicit none

  type(namelist_collection_type), intent(in) :: configuration
  type(config_type),              intent(in) :: config

  type(namelist_type), pointer :: aerosol_nml
  type(namelist_type), pointer :: base_mesh_nml
  type(namelist_type), pointer :: extrusion_nml
  type(namelist_type), pointer :: files_nml

  type(namelist_type), pointer :: planet_nml
  type(namelist_type), pointer :: io_nml
  type(namelist_type), pointer :: timestepping_nml
  type(namelist_type), pointer :: time_nml
  type(namelist_type), pointer :: mixed_solver_nml
  type(namelist_type), pointer :: microphysics_nml
  type(namelist_type), pointer :: multires_coupling_nml
  type(namelist_type), pointer :: multigrid_nml
  type(namelist_type), pointer :: formulation_nml
  type(namelist_type), pointer :: initialization_nml
  type(namelist_type), pointer :: boundaries_nml
  type(namelist_type), pointer :: finite_element_nml

! initial_temperature, bvf_square r_def
! gravity_wave_constants b_space i_def
  character(str_max_filename) :: start_dump_filename
  character(str_max_filename) :: checkpoint_stem_name
! lfric2lfric regrid_method i_def
  character(str_def) :: prime_mesh_name
  character(str_def) :: aerosol_mesh_name
  character(str_def) :: orography_mesh_name
  character(str_def) :: time_origin
  character(str_def) :: time_start
  character(str_max_filename) :: file_prefix

  logical(l_def) :: prepartitioned
  integer(i_def) :: geometry
  integer(i_def) :: moisture_formulation
  integer(i_def) :: lbc_option
  integer(i_def) :: ls_option
  integer(i_def) :: init_option
  integer(i_def) :: lbc_eos_height
  integer(i_def) :: model_eos_height

  real(r_def)    :: tau_r
  real(r_def)    :: atol
  real(r_def)    :: domain_height
  real(r_def)    :: planet_radius
  real(r_def)    :: scaled_radius
  real(r_second) :: dt
  integer(i_def) :: method
  integer(i_def) :: nlayers
  integer(i_def) :: element_order_h
  integer(i_def) :: element_order_v
  integer(i_medium) :: diag_freq
  logical(l_def) :: nodal
  logical(l_def) :: write_diag
  logical(l_def) :: use_xios_io
  logical(l_def) :: l_multigrid
  logical(l_def) :: use_multires_coupling
  logical(l_def) :: subroutine_timers
  logical(l_def) :: ozone_ancil
  logical(l_def) :: aero_ancil
  logical(l_def) :: murk_lbc
  logical(l_def) :: microphysics_casim
  logical(l_def) :: read_w2h_wind

!  character(str_def), allocatable :: chain_mesh_tags(:)
!  character(str_def), allocatable :: multires_coupling_mesh_tags(:)

  character(*), parameter :: message = 'Difference in config objects '

  call log_event('Validating Config/Configuration object data.', &
                 log_level_info)

  if (configuration%namelist_exists('base_mesh')) then
    base_mesh_nml => configuration%get_namelist('base_mesh')

    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    call base_mesh_nml%get_value( 'file_prefix', file_prefix )
    call base_mesh_nml%get_value( 'geometry', geometry )
    call base_mesh_nml%get_value( 'prepartitioned', prepartitioned )

    if (prime_mesh_name /= config%base_mesh%prime_mesh_name()) then
      call log_event(message//'prime_mesh_name', log_level_warning)
    end if

    if (file_prefix /= config%base_mesh%file_prefix()) then
      call log_event(message//'file_prefix', log_level_warning)
    end if

    if (geometry /= config%base_mesh%geometry()) then
      call log_event(message//'geometry', log_level_warning)
    end if

    if (prepartitioned .neqv. config%base_mesh%prepartitioned()) then
      call log_event(message//'prepartitioned', log_level_warning)
    end if
  end if

  if (configuration%namelist_exists('extrusion')) then
    extrusion_nml => configuration%get_namelist('extrusion')

    call extrusion_nml%get_value( 'method', method )
    call extrusion_nml%get_value( 'planet_radius', planet_radius )
    call extrusion_nml%get_value( 'domain_height', domain_height )
    call extrusion_nml%get_value( 'number_of_layers', nlayers )

    if (method /= config%extrusion%method()) then
      call log_event(message//'method', log_level_warning)
    end if

    if (planet_radius /= config%extrusion%planet_radius()) then
      call log_event(message//'planet_radius', log_level_warning)
    end if

    if (domain_height /= config%extrusion%domain_height()) then
      call log_event(message//'domain_height', log_level_warning)
    end if

    if (nlayers /= config%extrusion%number_of_layers()) then
      call log_event(message//'number_of_layers', log_level_warning)
    end if
  end if

  if (configuration%namelist_exists('io')) then
    io_nml => configuration%get_namelist('io')

    call io_nml%get_value( 'nodal_output_on_w3', nodal )
    call io_nml%get_value( 'write_diag', write_diag )
    call io_nml%get_value( 'use_xios_io', use_xios_io )
    call io_nml%get_value( 'subroutine_timers', subroutine_timers )
    call io_nml%get_value( 'diagnostic_frequency', diag_freq )

    if (diag_freq /= config%io%diagnostic_frequency()) then
      call log_event(message//'diagnostic_frequency', log_level_warning)
    end if

    if (nodal .neqv. config%io%nodal_output_on_w3()) then
      call log_event(message//'nodal_output_on_w3', log_level_warning)
    end if

    if (write_diag .neqv. config%io%write_diag()) then
      call log_event(message//'write_diag', log_level_warning)
    end if

    if (use_xios_io .neqv. config%io%use_xios_io()) then
      call log_event(message//'use_xios_io', log_level_warning)
    end if

    if (subroutine_timers .neqv. config%io%subroutine_timers()) then
      call log_event(message//'subroutine_timers', log_level_warning)
    end if
  end if

  if (configuration%namelist_exists('finite_element')) then
    finite_element_nml => configuration%get_namelist('finite_element')

    call finite_element_nml%get_value('element_order_h', element_order_h)
    call finite_element_nml%get_value('element_order_v', element_order_v)

    if (element_order_h /= config%finite_element%element_order_h()) then
      call log_event( message//'element_order_h', log_level_warning )
    end if

    if (element_order_v /= config%finite_element%element_order_v()) then
      call log_event( message//'element_order_v', log_level_warning )
    end if
  end if

  if (configuration%namelist_exists('formulation')) then
    formulation_nml => configuration%get_namelist('formulation')

    call formulation_nml%get_value('l_multigrid', l_multigrid)
    call formulation_nml%get_value('moisture_formulation', moisture_formulation)
    call formulation_nml%get_value('use_multires_coupling', use_multires_coupling)

    if (moisture_formulation /= config%formulation%moisture_formulation()) then
      call log_event( message//'moisture_formulation', log_level_warning )
    end if

    if (l_multigrid .neqv. config%formulation%l_multigrid()) then
      call log_event( message//'l_multigrid', log_level_warning )
    end if

    if (use_multires_coupling .neqv. config%formulation%use_multires_coupling()) then
      call log_event( message//'use_multires_coupling', log_level_warning )
    end if
  end if

  if (configuration%namelist_exists('planet')) then
    planet_nml => configuration%get_namelist('planet')

    call planet_nml%get_value( 'scaled_radius', scaled_radius )

    if (scaled_radius /= config%planet%scaled_radius()) then
      call log_event( message//'scaled_radius', log_level_warning )
    end if
  end if


  if (configuration%namelist_exists('timestepping')) then
    timestepping_nml => configuration%get_namelist('timestepping')

    call timestepping_nml%get_value( 'dt', dt )
    call timestepping_nml%get_value( 'tau_r', tau_r )

    if (dt /= config%timestepping%dt()) then
      call log_event( message//'dt', log_level_warning )
    end if

    if (tau_r /= config%timestepping%tau_r()) then
      call log_event( message//'tau_r', log_level_warning )
    end if

  end if


  if (configuration%namelist_exists('mixed_solver')) then
    mixed_solver_nml => configuration%get_namelist('mixed_solver')

    call mixed_solver_nml%get_value('mixed_solver_a_tol', atol)

    if (atol /= config%mixed_solver%mixed_solver_a_tol()) then
      call log_event( message//'mixed_solver_a_tol', log_level_warning )
    end if
  end if


  if (configuration%namelist_exists('initialization')) then
    initialization_nml => configuration%get_namelist('initialization')

    call initialization_nml%get_value('ls_option', ls_option)
    call initialization_nml%get_value('lbc_option', lbc_option)
    call initialization_nml%get_value('coarse_aerosol_ancil', aero_ancil)
    call initialization_nml%get_value('coarse_ozone_ancil', ozone_ancil)
    call initialization_nml%get_value('init_option', init_option)
    call initialization_nml%get_value('model_eos_height', model_eos_height)
    call initialization_nml%get_value('read_w2h_wind', read_w2h_wind)

    if (ls_option /= config%initialization%ls_option()) then
      call log_event( message//'ls_option', log_level_warning )
    end if

    if (lbc_option /= config%initialization%lbc_option()) then
      call log_event( message//'lbc_option', log_level_warning )
    end if

    if (init_option /= config%initialization%init_option()) then
      call log_event( message//'init_option', log_level_warning )
    end if

    if (aero_ancil .neqv. config%initialization%coarse_aerosol_ancil()) then
      call log_event( message//'aerosol_ancil', log_level_warning )
    end if

    if (ozone_ancil .neqv. config%initialization%coarse_ozone_ancil()) then
      call log_event( message//'ozone_ancil', log_level_warning )
    end if

    if (model_eos_height /= config%initialization%model_eos_height()) then
      call log_event( message//'model_eos_height', log_level_warning )
    end if

    if (read_w2h_wind .neqv. config%initialization%read_w2h_wind()) then
      call log_event( message//'read_w2h_wind', log_level_warning )
    end if

  end if


  if (configuration%namelist_exists('boundaries')) then
    boundaries_nml => configuration%get_namelist('boundaries')

    call boundaries_nml%get_value('lbc_eos_height', lbc_eos_height)

    if (lbc_eos_height /= config%boundaries%lbc_eos_height()) then
      call log_event( message//'lbc_eos_height', log_level_warning )
    end if
  end if

  if (configuration%namelist_exists('multires_coupling')) then
    multires_coupling_nml => configuration%get_namelist('multires_coupling')

    call multires_coupling_nml%get_value('aerosol_mesh_name', aerosol_mesh_name)
    call multires_coupling_nml%get_value('orography_mesh_name', orography_mesh_name)
!   call multires_coupling_nml%get_value('multires_coupling_mesh_tags', multires_coupling_mesh_tags)

    if (orography_mesh_name /= config%multires_coupling%orography_mesh_name()) then
      call log_event( message//'orography_mesh_name', log_level_warning )
    end if

!!$    if (multires_coupling_mesh_tags /= config%multires_coupling%multires_coupling_mesh_tags()) then
!!$      call log_event( message//'multires_coupling_mesh_tags', log_level_warning )
!!$    end if

    if (aerosol_mesh_name /= config%multires_coupling%aerosol_mesh_name()) then
      call log_event( message//'aerosol_mesh_name', log_level_warning )
    end if
  end if

  if (configuration%namelist_exists('aerosol')) then
    aerosol_nml => configuration%get_namelist('aerosol')

    call aerosol_nml%get_value('murk_lbc', murk_lbc)

    if (murk_lbc .neqv. config%aerosol%murk_lbc()) then
      call log_event( message//'murk_lbc', log_level_warning )
    end if
  end if

  if (configuration%namelist_exists('multigrid')) then
    multigrid_nml => configuration%get_namelist('multigrid')

!!$    call multigrid_nml%get_value('chain_mesh__tags', chain_mesh_tags)
!!$    chain_mesh_tags_2 = config%multigrid%chain_mesh_tags()
!!$    do i=1, size(chain_mesh_tags)
!!$      if (chain_mesh_tags /= config%multigrid%chain_mesh_tags()) then
!!$        call log_event( message//'chain_mesh_tags', log_level_warning )
!!$      end if
!!$    end do
   end if

  if (configuration%namelist_exists('microphysics')) then
    microphysics_nml => configuration%get_namelist('microphysics')

    call microphysics_nml%get_value('microphysics_casim', microphysics_casim)

    if (microphysics_casim .neqv. config%microphysics%microphysics_casim()) then
      call log_event( message//'microphysics_casim', log_level_warning )
    end if
  end if

  if (configuration%namelist_exists('time')) then
    time_nml => configuration%get_namelist('time')

    call time_nml%get_value('calendar_origin', time_origin)
    call time_nml%get_value('calendar_start', time_start)

    if (time_origin /= config%time%calendar_origin()) then
      call log_event( message//'calendar_origin', log_level_warning )
    end if

    if (time_start /= config%time%calendar_start())  then
      call log_event( message//'calendar_start', log_level_warning )
    end if
  end if

  if (configuration%namelist_exists('files')) then
    files_nml => configuration%get_namelist('files')

    call files_nml%get_value('start_dump_filename', start_dump_filename)
    call files_nml%get_value('checkpoint_stem_name', checkpoint_stem_name)

    if (start_dump_filename /= config%files%start_dump_filename()) then
      call log_event( message//'start_dump_filename', log_level_warning )
    end if

    if (checkpoint_stem_name /= config%files%checkpoint_stem_name())  then
      call log_event( message//'checkpoint_stem_name', log_level_warning )
    end if
  end if

end subroutine check_config_api

end module check_config_api_mod
