!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Function space discovery from XIOS metadata

module space_from_metadata_mod

  use log_mod,                         only: log_event,         &
                                             log_scratch_space, &
                                             log_level_trace,   &
                                             log_level_error
  use constants_mod,                   only: i_def, l_def, str_def
  use mesh_mod,                        only: mesh_type
  use fs_continuity_mod,               only: W3,                      &
                                             Wtheta,                  &
                                             W2H,                     &
                                             W2,                      &
                                             W0,                      &
                                             name_from_functionspace
  use function_space_collection_mod,   only: function_space_collection
  use function_space_mod,              only: function_space_type
  use lfric_xios_diag_mod,             only: get_field_order,      &
                                             get_field_grid_ref,   &
                                             get_field_domain_ref, &
                                             get_field_axis_ref
  use multidata_field_dimensions_mod,  only: get_ndata => &
                                                  get_multidata_field_dimension
  use extrusion_mod,                   only: TWOD
  use mesh_collection_mod,             only: mesh_collection
  use base_mesh_config_mod,            only: prime_mesh_name

  implicit none

  private

  ! field flavours
  character(str_def), parameter :: vanilla       = 'VanillaField'  ! scalar field on 3D mesh
  character(str_def), parameter :: vanilla_multi = 'VanillaMulti'  ! multidata field on 3D mesh
  character(str_def), parameter :: planar        = 'PlanarField'
  character(str_def), parameter :: tile          = 'TileField'
  character(str_def), parameter :: radiation     = 'RadiationField'

  ! grid names
  character(str_def), parameter :: full_level_face_grid                        &
    = 'full_level_face_grid'
  character(str_def), parameter :: var_full_face_grid                          &
    = 'var_full_face_grid'
  character(str_def), parameter :: var_face                                    &
    = 'var_face'
  character(str_def), parameter :: half_level_face_grid                        &
    = 'half_level_face_grid'
  character(str_def), parameter :: half_level_edge_grid                        &
    = 'half_level_edge_grid'
  character(str_def), parameter :: node_grid                                   &
    = 'node_grid'

  public :: space_from_metadata

contains

  ! if string is of shape "<prefix>_<suffix>",
  ! split it into prefix and suffix, updating the string
  function try_split(string, suffix, prefix) result(ok)
    character(*), intent(inout) :: string
    character(*), intent(inout) :: suffix
    character(*), intent(in) :: prefix
    logical(l_def) :: ok
    integer(i_def) :: m
    if (index(string, trim(prefix) // '_') == 1) then
      m = len_trim(prefix)
      suffix = string(m+2:) ! start one after the underscore
      string = prefix
      ok = .true.
    else
      ok = .false.
    end if
  end function try_split

  ! extract and remove axis reference from composite grid reference;
  ! for instance, "full_level_face_grid_cloud_subcols" yields
  !   grid_ref = "full_level_face_grid"
  !   axis_ref = "cloud_subcols"
  subroutine split_composite_grid_ref(grid_ref, axis_ref)
    implicit none
    character(*), intent(inout) :: grid_ref
    character(*), intent(out) :: axis_ref
    if (try_split(grid_ref, axis_ref, full_level_face_grid)) then
      ! noop
    else if (try_split(grid_ref, axis_ref, var_full_face_grid)) then
      ! noop
    else if (try_split(grid_ref, axis_ref, half_level_face_grid)) then
      ! noop
    else if (try_split(grid_ref, axis_ref, half_level_edge_grid)) then
      ! noop
    else if (try_split(grid_ref, axis_ref, node_grid)) then
      ! noop
    else if (try_split(grid_ref, axis_ref, var_face)) then
      ! noop
    else
      ! no axis_ref - keep grid_ref unchanged
      axis_ref = ''
    end if
  end subroutine split_composite_grid_ref

  ! derive function space enumerator from xios metadata
  function get_field_fsenum(xios_id, grid_ref, domain_ref)                    &
    result(fsenum)
    implicit none
    character(*), intent(in) :: xios_id
    character(*), intent(in) :: grid_ref
    character(*), intent(in) :: domain_ref
    integer(i_def) :: fsenum

    ! from RB's python metadata generator
    if (grid_ref == full_level_face_grid .or. &
        grid_ref == var_full_face_grid) then
      fsenum = Wtheta
    else if (grid_ref == half_level_face_grid                                 &
      .or. grid_ref == var_face                                               &
      .or. domain_ref == 'face') then
      fsenum = W3
    else if (grid_ref == half_level_edge_grid) then
      fsenum = W2H
    else if (grid_ref == node_grid) then
      fsenum = W0
    else if (domain_ref == "checkpoint_Wtheta") then
      fsenum = Wtheta
    else if (domain_ref == "checkpoint_W3") then
      fsenum = W3
    else if (domain_ref == "checkpoint_W2") then
      fsenum = W2
    else
      fsenum = 0 ! silence compiler warning
      write(log_scratch_space, *)                                             &
        'cannot derive function space enumerator for field: ' //              &
        trim(xios_id) //                                                      &
        ', grid_ref: ' // trim(grid_ref) //                                   &
        ', domain_ref: ' // trim(domain_ref)
      call log_event(log_scratch_space, log_level_error)
    end if
  end function get_field_fsenum

  ! derive field flavour (vanilla/planar/tile/radiation) from xios metadata
  function get_field_flavour(xios_id, grid_ref, domain_ref, axis_ref)         &
    result(flavour)
    implicit none
    character(*), intent(in) :: xios_id
    character(*), intent(in) :: grid_ref
    character(*), intent(in) :: domain_ref
    character(*), intent(in) :: axis_ref
    character(str_def) :: flavour

    if (grid_ref /= "") then
      if (domain_ref /= "") then
        write(log_scratch_space, *)                                           &
        'field ' // trim(xios_id) //                                          &
        'with grid_ref and domain_ref : ' //                                  &
        grid_ref // ' ' // domain_ref
        call log_event(log_scratch_space, log_level_error)
      end if
      if (axis_ref /= "") then
        flavour = vanilla_multi
      else
        flavour = vanilla
      end if
    else
      if (domain_ref /= "") then
        if (axis_ref /= "") then
          if (axis_ref == 'radiation_levels') then
            flavour = radiation
          else
            flavour = tile
          end if
        else
          ! only domain - must not happen except for face domain
          if (domain_ref == 'face') then
            flavour = planar
          else if (domain_ref == "checkpoint_Wtheta" .or. &
                   domain_ref == "checkpoint_W3" .or.     &
                   domain_ref == "checkpoint_W2"          ) then
            flavour = vanilla
          else
            write(log_scratch_space, *)                                       &
            'field ' // trim(xios_id) //                                      &
            ' with only unknown domain_ref: ' // domain_ref
            call log_event(log_scratch_space, log_level_error)
          end if
        end if
      else
        ! only axis - must not happen
        write(log_scratch_space, *)                                           &
        'field ' // trim(xios_id) //                                          &
        ' with only an axis_ref: ' // axis_ref
        call log_event(log_scratch_space, log_level_error)
      end if
    end if
  end function get_field_flavour

  ! derive multidata item name from xios metadata
  function get_field_tile_id(xios_id, axis_ref) result(tile_id)
    implicit none
    character(*), intent(in) :: xios_id
    character(*), intent(in) :: axis_ref
    character(str_def) :: tile_id
    tile_id = axis_ref
  end function get_field_tile_id

  !> @brief Space discovery from metadata for field given by XIOS id
  !> @param[in]           xios_id          XIOS id of field
  !> @param[in]           status           Field status (logging only)
  !> @param[in, optional] mesh_3d          Override 3D mesh
  !> @param[in, optional] mesh_2d          Override 2D mesh
  !> @param[in, optional] force_mesh       Override derived mesh
  !> @param[in, optional] force_rad_levels Override derived radiation levels
  !> @param[in, optional] force_ndata      Override derived ndata
  !> @param[in, optional] force_order_h    Override function space order in
  !!                                       horizontal
  !> @param[in, optional] force_order_v    Override function space order in
  !!                                       vertical
  !> @return                               Function space returned
  function space_from_metadata(xios_id, status,                            &
                               mesh_3d, mesh_2d, force_mesh,               &
                               force_rad_levels, force_order_h,            &
                               force_order_v, force_ndata) &
                               result(vector_space)

    implicit none

    character(*),                       intent(in)  :: xios_id
    character(*),                       intent(in)  :: status
    type(mesh_type), pointer, optional, intent(in)  :: mesh_3d
    type(mesh_type), pointer, optional, intent(in)  :: mesh_2d
    type(mesh_type), pointer, optional, intent(in)  :: force_mesh
    integer(kind=i_def), optional,      intent(in)  :: force_rad_levels
    integer(kind=i_def), optional,      intent(in)  :: force_order_h
    integer(kind=i_def), optional,      intent(in)  :: force_order_v
    integer(kind=i_def), optional,      intent(in)  :: force_ndata

    character(str_def)  :: grid_ref
    character(str_def)  :: domain_ref
    character(str_def)  :: axis_ref
    integer(kind=i_def) :: order_h
    integer(kind=i_def) :: order_v
    integer(kind=i_def) :: fsenum
    character(str_def)  :: flavour
    integer(kind=i_def) :: ndata

    type(mesh_type), pointer :: diag_mesh_3d
    type(mesh_type), pointer :: diag_mesh_2d

    type(mesh_type), pointer            :: this_mesh
    type(function_space_type), pointer  :: vector_space

#ifdef UNIT_TEST
    diag_mesh_3d => mesh_collection%get_mesh('test mesh: planar bi-periodic')
    diag_mesh_2d => mesh_collection%get_mesh('test mesh: planar bi-periodic')
#else
    ! if specific meshes have been supplied, override the prime mesh
    if (present(mesh_3d)) then
      diag_mesh_3d => mesh_3d
    else
      ! default to prime_mesh
      diag_mesh_3d => mesh_collection%get_mesh(prime_mesh_name)
    endif
    if (present(mesh_2d)) then
      diag_mesh_2d => mesh_2d
    else
      diag_mesh_2d => mesh_collection%get_mesh_variant(diag_mesh_3d, &
                                                       extrusion_id=TWOD)
    endif
#endif

    ! metadata lookup
    grid_ref = get_field_grid_ref(xios_id)
    call split_composite_grid_ref(grid_ref, axis_ref)
    domain_ref = get_field_domain_ref(xios_id)
    if (axis_ref == '') then
      ! no axis encoded in grid name - try to get it from XIOS
      axis_ref = get_field_axis_ref(xios_id)
    end if

    ! Currently just returns 0 unless input is non-zero
    order_h = get_field_order(xios_id, force_order_h)
    order_v = get_field_order(xios_id, force_order_v)

    ! derive function space and flavour from metadata
    fsenum = get_field_fsenum(xios_id, grid_ref, domain_ref)
    flavour = get_field_flavour(xios_id, grid_ref, domain_ref, axis_ref)

    ! select mesh
    if (present(force_mesh)) then
      this_mesh => force_mesh
    else if (                                                                 &
      flavour == radiation .or.                                               &
      flavour == tile .or.                                                    &
      flavour == planar) then
      this_mesh => diag_mesh_2d
    else if (flavour == vanilla .or. flavour == vanilla_multi) then
      this_mesh => diag_mesh_3d
    else
      call log_event( "unexpected flavour: " // flavour, log_level_error )
    end if

    ! derive ndata (multidata dimension)
    select case (flavour)
    case (radiation)
     ! radiation fields are treated as special multidata fields
      if (present(force_rad_levels)) then
        ndata = force_rad_levels
      else
        ndata = diag_mesh_3d%get_nlayers() + 1
      end if
    case (tile)
      ! genuine multidata field - derive dimension from axis metadata
      if (axis_ref == "") then
        ndata = 1
      else
        ndata = get_ndata(get_field_tile_id(xios_id, axis_ref))
      end if
    case (planar)
      ! scalar field on 2d mesh
      ndata = 1
    case (vanilla )
      ! scalar field on 3d mesh
      ndata = 1
    case (vanilla_multi)
      ! multidata field on 3d mesh
      ndata = get_ndata(get_field_tile_id(xios_id, axis_ref))
    case default
      call log_event("unexpected flavour: " // flavour, log_level_error)
    end select

    if (present(force_ndata)) then
      ! we are pre-initialising a field with the time window size,
      ! which is currently not supported for multidata fields;
      ! should the need arise, one could multiply ndata with force_ndata
      if (ndata /= 1) then
        call log_event(                                                      &
          "cannot override ndata for multidata field with a time window",    &
        log_level_error)
      end if
      ndata = force_ndata
    end if

    ! set up function space - needed by psyclone even for inactive fields
    vector_space => function_space_collection%get_fs(                         &
      this_mesh,                                                              &
      order_h,                                                                &
      order_v,                                                                &
      fsenum,                                                                 &
      ndata)
#ifndef UNIT_TEST
    write(log_scratch_space,                                                  &
       '("field: ", A, ", ", A, '                                             &
       // '", mesh: ", A, ", order_h: ", I2, ", order_v: ", I2, '             &
       // '", fs: ", A, ", flavour: ", A, ", ndata: ", I5)')                  &
      trim(xios_id), trim(status), trim(this_mesh%get_mesh_name()), order_h,  &
      order_v, trim(name_from_functionspace(fsenum)), trim(flavour), ndata
    call log_event(log_scratch_space, log_level_trace)
#endif
    ! paranoia
    nullify(this_mesh)
end function space_from_metadata

end module space_from_metadata_mod
