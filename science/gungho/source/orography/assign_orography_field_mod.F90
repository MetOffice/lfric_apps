!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Module to assign the values of the surface height to model
!> coordinates using either an analytic orography function or from
!> a surface_altitude field.
!> Note that unlike other algorithms, this is breaks encapsulation in order to
!> write to the chi field.  This is an exception and only allowed in the set up
!> phase of the model.  Generally, the chi field is read only (and this is
!> enforced through PSyClone).
!-------------------------------------------------------------------------------
module assign_orography_field_mod

  use constants_mod,                  only : r_def, i_def, l_def
  use orography_config_mod,           only : orog_init_option, &
                                             orog_init_option_analytic, &
                                             orog_init_option_ancil, &
                                             orog_init_option_none, &
                                             orog_init_option_start_dump
  use base_mesh_config_mod,           only : geometry, &
                                             geometry_spherical
  use finite_element_config_mod,      only : coord_system,            &
                                             coord_order,             &
                                             coord_system_xyz,        &
                                             coord_system_alphabetaz, &
                                             coord_system_lonlatz
  use mesh_collection_mod,            only : mesh_collection
  use coord_transform_mod,            only : xyz2llr, llr2xyz, alphabetar2llr
  use orography_helper_functions_mod, only : z2eta_linear, &
                                             eta2z_linear, &
                                             eta2z_smooth
  use analytic_orography_mod,         only : orography_profile
  use extrusion_config_mod,           only : stretching_height, &
                                             stretching_method, &
                                             stretching_method_linear
  use log_mod,                        only : log_event,      &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use fs_continuity_mod,              only : W0, Wchi
  use function_space_mod,             only : BASIS
  use surface_altitude_alg_mod,       only : surface_altitude_alg

  implicit none

  private

  public :: assign_orography_field
  ! These are made public only for unit-testing
  public :: analytic_orography_spherical_xyz
  public :: analytic_orography_spherical_alphabetaz
  public :: analytic_orography_spherical_lonlatz
  public :: analytic_orography_cartesian
  public :: ancil_orography_spherical_xyz
  public :: ancil_orography_spherical_sph
  public :: ancil_orography_cartesian

  interface

    subroutine analytic_orography_interface(nlayers,                     &
                                            ndf_chi, undf_chi, map_chi,  &
                                            ndf_pid, undf_pid, map_pid,  &
                                            domain_surface, domain_top,  &
                                            chi_1, chi_2, chi_3, panel_id)
      import :: r_def, i_def
      implicit none
      integer(kind=i_def), intent(in) :: nlayers, undf_chi, undf_pid
      integer(kind=i_def), intent(in) :: ndf_chi, ndf_pid
      integer(kind=i_def), intent(in) :: map_chi(ndf_chi), map_pid(ndf_pid)
      real(kind=r_def),    intent(in) :: domain_surface, domain_top
      real(kind=r_def), intent(inout) :: chi_1(undf_chi), chi_2(undf_chi), chi_3(undf_chi)
      real(kind=r_def),    intent(in) :: panel_id(undf_pid)
    end subroutine analytic_orography_interface

  end interface

  interface

    subroutine ancil_orography_interface(nlayers,                    &
                                         chi_1, chi_2, chi_3,        &
                                         panel_id,                   &
                                         surface_altitude,           &
                                         domain_surface, domain_top, &
                                         ndf_chi, undf_chi,          &
                                         map_chi,                    &
                                         ndf_pid, undf_pid,          &
                                         map_pid,                    &
                                         ndf, undf,                  &
                                         map, basis                  &
                                         )
      import :: r_def, i_def
      implicit none
      integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_chi, ndf_pid
      integer(kind=i_def), intent(in) :: undf, undf_chi, undf_pid
      integer(kind=i_def), dimension(ndf),     intent(in)   :: map
      integer(kind=i_def), dimension(ndf_chi), intent(in)   :: map_chi
      integer(kind=i_def), dimension(ndf_pid), intent(in)   :: map_pid
      real(kind=r_def), intent(in), dimension(ndf, ndf_chi) :: basis
      real(kind=r_def), dimension(undf_chi), intent(inout)  :: chi_1, chi_2, chi_3
      real(kind=r_def), dimension(undf_pid), intent(in)     :: panel_id
      real(kind=r_def), dimension(undf),     intent(in)     :: surface_altitude
      real(kind=r_def), intent(in)                          :: domain_surface, domain_top
    end subroutine ancil_orography_interface

  end interface

contains

  !=============================================================================
  !> @brief Updates model vertical coordinate using selected analytic orography.
  !>
  !> @details Model coordinate array of size 3 for the type field is passed in
  !> to be updated. The field proxy is used to break encapsulation and access
  !> the function space and the data attributes of the field so that its values
  !> can be updated. Model coordinates are updated by calling single column
  !> subroutines, one for spherical and the other for Cartesian domain. These
  !> routines calculate analytic orography from horizontal coordinates or else
  !> use the surface_altitude field and then update the vertical coordinate.
  !>
  !> @param[in,out] chi_inventory      Contains all of the model's coordinate
  !!                                   fields, itemised by mesh
  !> @param[in]     panel_id_inventory Contains all of the model's panel ID
  !!                                   fields, itemised by mesh
  !> @param[in]     mesh               Mesh to apply orography to
  !> @param[in]     surface_altitude   Field containing the surface altitude
  !=============================================================================
  subroutine assign_orography_field(chi_inventory, panel_id_inventory, &
                                    mesh, surface_altitude)

    use inventory_by_mesh_mod,          only : inventory_by_mesh_type
    use field_mod,                      only : field_type, field_proxy_type
    use mesh_mod,                       only : mesh_type
    use domain_mod,                     only : domain_type
    use orography_helper_functions_mod, only : set_horizontal_domain_size
    use function_space_collection_mod,  only : function_space_collection
    use orography_config_mod,           only : orog_init_option, &
                                               orog_init_option_ancil, &
                                               orog_init_option_start_dump
    use field_minmax_alg_mod,           only : log_field_minmax

    implicit none

    ! Arguments
    type( inventory_by_mesh_type ),  intent( inout )       :: chi_inventory
    type( inventory_by_mesh_type ),  intent( in )          :: panel_id_inventory
    type( mesh_type ),               intent( in ), pointer :: mesh

    ! We keep the surface_altitude as an optional argument since it is
    ! not needed for miniapps that only want analytic orography
    type( field_type ),  intent( in ), optional :: surface_altitude

    ! Local variables
    type( field_type ),  pointer :: chi(:) => null()
    type( field_type ),  pointer :: panel_id => null()
    type( field_proxy_type )     :: chi_proxy(3)
    type( field_proxy_type )     :: panel_id_proxy
    type( domain_type )          :: domain
    real(kind=r_def)             :: domain_top, domain_surface
    integer(kind=i_def)          :: cell
    integer(kind=i_def)          :: undf_chi, ndf_chi, nlayers
    integer(kind=i_def)          :: undf_pid, ndf_pid
    integer(kind=i_def)          :: undf_sf, ndf_sf
    integer(kind=i_def), pointer :: map_chi(:,:) => null()
    integer(kind=i_def), pointer :: map_pid(:,:) => null()
    integer(kind=i_def), pointer :: map_sf(:,:) => null()

    integer(kind=i_def)          :: surface_order
    type( mesh_type ), pointer   :: sf_mesh
    type( field_type )           :: surface_altitude_w0
    type( field_proxy_type )     :: sfc_alt_proxy

    real(kind=r_def), pointer :: nodes(:,:) => null()
    integer(kind=i_def) :: dim_sf, df, df_sf

    ! Procedure pointer
    procedure(analytic_orography_interface), pointer :: analytic_orography => null()
    procedure(ancil_orography_interface), pointer    :: ancil_orography => null()

    real(kind=r_def), allocatable :: basis_sf_on_chi(:,:,:)

    if (coord_order == 0 .and. orog_init_option/=orog_init_option_none) then
      call log_event( "assign_orography_field: "// &
         "Orography assignment is currently only available with coord_order > 0.", &
         LOG_LEVEL_ERROR )
    end if

    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)

    ! Get domain size
    domain = mesh%get_domain()
    call set_horizontal_domain_size( domain )

    ! Get physical height of flat domain surface from the domain_size object
    domain_surface = domain%get_base_height()

    ! Get domain top from the mesh object and domain_surface
    domain_top = mesh%get_domain_top() + domain_surface

    if (orog_init_option==orog_init_option_none) then

      call log_event( "assign_orography_field: "// &
         "Flat surface requested.", LOG_LEVEL_INFO )

    else if (orog_init_option==orog_init_option_analytic) then

      call log_event( "assign_orography_field: "// &
         "Assigning analytic orography.", LOG_LEVEL_INFO )

      ! Point to appropriate procedure to assign orography
      if ( geometry == geometry_spherical ) then
        if ( coord_system == coord_system_xyz ) then
          analytic_orography => analytic_orography_spherical_xyz
        else if ( coord_system == coord_system_alphabetaz ) then
          analytic_orography => analytic_orography_spherical_alphabetaz
        else if ( coord_system == coord_system_lonlatz ) then
          analytic_orography => analytic_orography_spherical_lonlatz
        else
          call log_event("Error: this coordinate system is not " // &
                         "implemented with analytic orography", &
                         LOG_LEVEL_ERROR)
        end if
      else
        analytic_orography => analytic_orography_cartesian
      end if

      ! Break encapsulation and get the proxy
      chi_proxy(1)   = chi(1)%get_proxy()
      chi_proxy(2)   = chi(2)%get_proxy()
      chi_proxy(3)   = chi(3)%get_proxy()
      undf_chi       = chi_proxy(1)%vspace%get_undf()
      ndf_chi        = chi_proxy(1)%vspace%get_ndf()
      panel_id_proxy = panel_id%get_proxy()
      undf_pid       = panel_id_proxy%vspace%get_undf()
      ndf_pid        = panel_id_proxy%vspace%get_ndf()
      nlayers        = chi_proxy(1)%vspace%get_nlayers()

      map_chi => chi_proxy(1)%vspace%get_whole_dofmap()
      map_pid => panel_id_proxy%vspace%get_whole_dofmap()

      ! Call column procedure
      do cell = 1,chi_proxy(1)%vspace%get_ncell()

        call analytic_orography(nlayers,           &
                                ndf_chi,           &
                                undf_chi,          &
                                map_chi(:,cell),   &
                                ndf_pid,           &
                                undf_pid,          &
                                map_pid(:,cell),   &
                                domain_surface,    &
                                domain_top,        &
                                chi_proxy(1)%data, &
                                chi_proxy(2)%data, &
                                chi_proxy(3)%data, &
                                panel_id_proxy%data )
      end do


    else if  (orog_init_option==orog_init_option_ancil .or. &
              orog_init_option==orog_init_option_start_dump) then

      call log_event( "assign_orography_field: "// &
         "Assigning orography from surface_altitude field.", LOG_LEVEL_INFO )

      ! Point to appropriate procedure to assign orography
      if ( geometry == geometry_spherical ) then
        if ( coord_system == coord_system_xyz ) then
          ancil_orography => ancil_orography_spherical_xyz
        else if ( (coord_system == coord_system_alphabetaz) .or. &
                  (coord_system == coord_system_lonlatz) ) then
          ancil_orography => ancil_orography_spherical_sph
        else
          call log_event("Error: this coordinate system is not " //  &
                         "implemented with ancillary orography", &
                         LOG_LEVEL_ERROR)
        end if
      else
        ancil_orography => ancil_orography_cartesian
      end if

      if ( present(surface_altitude) ) then

        ! Set up the surface altitude field on W0 points
        sf_mesh =>  surface_altitude%get_mesh()
        surface_order = surface_altitude%get_element_order()
        call surface_altitude_w0%initialise( vector_space =  &
           function_space_collection%get_fs(sf_mesh, surface_order, W0), &
                                             halo_depth = sf_mesh%get_halo_depth() )

        if (surface_altitude%which_function_space()==W0) then
          call surface_altitude%copy_field_serial(surface_altitude_w0)
        else
          call surface_altitude_alg( surface_altitude_w0, surface_altitude )
        end if

        call log_field_minmax( LOG_LEVEL_INFO, 'srf_alt', surface_altitude )
        call log_field_minmax( LOG_LEVEL_INFO, 'srf_alt_w0', &
                                                      surface_altitude_w0 )

        nullify ( sf_mesh )
      end if

      ! Break encapsulation and get the proxy
      chi_proxy(1) = chi(1)%get_proxy()
      chi_proxy(2) = chi(2)%get_proxy()
      chi_proxy(3) = chi(3)%get_proxy()
      panel_id_proxy = panel_id%get_proxy()
      sfc_alt_proxy = surface_altitude_w0%get_proxy()

      undf_chi = chi_proxy(1)%vspace%get_undf()
      ndf_chi  = chi_proxy(1)%vspace%get_ndf()
      undf_pid = panel_id_proxy%vspace%get_undf()
      ndf_pid  = panel_id_proxy%vspace%get_ndf()
      nlayers  = chi_proxy(1)%vspace%get_nlayers()
      undf_sf  = sfc_alt_proxy%vspace%get_undf()
      ndf_sf   = sfc_alt_proxy%vspace%get_ndf()

      map_chi => chi_proxy(1)%vspace%get_whole_dofmap()
      map_sf => sfc_alt_proxy%vspace%get_whole_dofmap()
      map_pid => panel_id_proxy%vspace%get_whole_dofmap()

      dim_sf = sfc_alt_proxy%vspace%get_dim_space()
      nodes => chi_proxy(1)%vspace%get_nodes()
      allocate(basis_sf_on_chi(dim_sf, ndf_sf, ndf_chi))
      !
      ! Compute basis/diff-basis arrays
      !
      do df=1,ndf_chi
        do df_sf=1,ndf_sf
          basis_sf_on_chi(:,df_sf,df) = &
             sfc_alt_proxy%vspace%call_function(BASIS,df_sf,nodes(:,df))
        end do
      end do

      ! Ensure halo is clean
      call sfc_alt_proxy%halo_exchange( &
                             depth=sfc_alt_proxy%get_field_proxy_halo_depth())

      ! Call column procedure
      do cell = 1,chi_proxy(1)%vspace%get_ncell()

        call ancil_orography(nlayers,                    &
                             chi_proxy(1)%data,          &
                             chi_proxy(2)%data,          &
                             chi_proxy(3)%data,          &
                             panel_id_proxy%data,        &
                             sfc_alt_proxy%data,         &
                             domain_surface, domain_top, &
                             ndf_chi, undf_chi,          &
                             map_chi(:,cell),            &
                             ndf_pid, undf_pid,          &
                             map_pid(:,cell),            &
                             ndf_sf, undf_sf,            &
                             map_sf(:,cell),             &
                             basis_sf_on_chi)
      end do

      deallocate(basis_sf_on_chi)

    end if

    nullify(chi, panel_id)

  end subroutine assign_orography_field

  !=============================================================================
  !> @brief Updates spherical vertical coordinate for a single column using
  !>        selected analytic orography.
  !>
  !> @details Calculates analytic orography from chi_1 and chi_2 horizontal
  !>          coordinates and then updates chi_3. As model coordinates for
  !>          spherical domain are currently (x,y,z) form they first need to be
  !>          converted to (long,lat,r) to assign orography to the model surface.
  !>          After evaluation of the new surface height chi_3 is updated using
  !>          its nondimensional eta coordinate and then transformed back to
  !>          (x,y,z) form.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in]     ndf_chi        Array size and loop bound for map_chi
  !> @param[in]     undf_chi       Column coordinates' array size and loop bound
  !> @param[in]     map_chi        Indirection map for coordinate field
  !> @param[in]     ndf_pid        Array size and loop bound for map_pid
  !> @param[in]     undf_pid       Panel ID array size and loop bound
  !> @param[in]     map_pid        Indirection map for panel_id
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in,out] chi_1          1st coordinate field in Wchi
  !> @param[in,out] chi_2          2nd coordinate field in Wchi
  !> @param[in,out] chi_3          3rd coordinate field in Wchi
  !> @param[in]     panel_id       Field giving the ID for mesh panels
  !=============================================================================
  subroutine analytic_orography_spherical_xyz(nlayers, ndf_chi, undf_chi, map_chi, &
                                              ndf_pid, undf_pid, map_pid,          &
                                              domain_surface, domain_top,          &
                                              chi_1, chi_2, chi_3, panel_id)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers, undf_chi, undf_pid
    integer(kind=i_def), intent(in)    :: ndf_chi, ndf_pid
    integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
    integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
    real(kind=r_def),    intent(in)    :: domain_surface, domain_top
    real(kind=r_def),    intent(inout) :: chi_1(undf_chi), chi_2(undf_chi), chi_3(undf_chi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    ! Internal variables
    integer(kind=i_def) :: k, df, dfk
    real(kind=r_def)    :: chi_3_r
    real(kind=r_def)    :: surface_height, domain_depth
    real(kind=r_def)    :: eta
    real(kind=r_def)    :: longitude, latitude, r

    domain_depth = domain_top - domain_surface

    ! Calculate orography and update chi_3
    do df = 1, ndf_chi
      do k = 0, nlayers-1
        dfk = map_chi(df)+k

        ! Model coordinates for spherical domain are in (x,y,z) form so they need
        ! to be converted to (long,lat,r) first
        call xyz2llr(chi_1(dfk), chi_2(dfk), chi_3(dfk), longitude, latitude, r)

        ! Calculate surface height for each DoF using selected analytic orography
        surface_height = orography_profile%analytic_orography(longitude, latitude)

        ! Calculate nondimensional coordinate from current height coordinate
        ! (chi_3) with flat domain_surface
        eta = z2eta_linear(r, domain_surface, domain_top)

        select case(stretching_method)
        case(stretching_method_linear)
          chi_3_r = eta2z_linear(eta, domain_surface + surface_height, domain_top)
        case default
          chi_3_r = domain_surface + &
             eta2z_smooth(eta, surface_height, domain_depth, stretching_height)
        end select

        ! Convert spherical coordinates back to model (x,y,z) form
        call llr2xyz(longitude, latitude, chi_3_r, &
                     chi_1(dfk), chi_2(dfk), chi_3(dfk))

      end do
    end do

  end subroutine analytic_orography_spherical_xyz

  !=============================================================================
  !> @brief Updates (alpha, beta, height) vertical coordinates for a single
  !>        column using selected analytic orography.
  !>
  !> @details Calculates analytic orography from chi_1 and chi_2 horizontal
  !>          coordinates and then updates chi_3. This works directly on the
  !>          cubed sphere (alpha,beta,r) coordinates.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in]     ndf_chi        Array size and loop bound for map_chi
  !> @param[in]     undf_chi       Column coordinates' array size and loop bound
  !> @param[in]     map_chi        Indirection map for coordinate field
  !> @param[in]     ndf_pid        Array size and loop bound for map_pid
  !> @param[in]     undf_pid       Panel ID array size and loop bound
  !> @param[in]     map_pid         Indirection map for panel_id
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in,out] chi_1          1st coordinate field in Wchi
  !> @param[in,out] chi_2          2nd coordinate field in Wchi
  !> @param[in,out] chi_3          3rd coordinate field in Wchi
  !> @param[in]     panel_id       Field giving the ID for mesh panels
  !=============================================================================
  subroutine analytic_orography_spherical_alphabetaz(nlayers,        &
                                                     ndf_chi,        &
                                                     undf_chi,       &
                                                     map_chi,        &
                                                     ndf_pid,        &
                                                     undf_pid,       &
                                                     map_pid,        &
                                                     domain_surface, &
                                                     domain_top,     &
                                                     chi_1,          &
                                                     chi_2,          &
                                                     chi_3,          &
                                                     panel_id)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers, undf_chi, undf_pid
    integer(kind=i_def), intent(in)    :: ndf_chi, ndf_pid
    integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
    integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
    real(kind=r_def),    intent(in)    :: domain_surface, domain_top
    real(kind=r_def),    intent(inout) :: chi_1(undf_chi), chi_2(undf_chi), chi_3(undf_chi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    ! Internal variables
    integer(kind=i_def) :: k, df, dfk, ipanel
    real(kind=r_def)    :: surface_height, domain_depth
    real(kind=r_def)    :: eta
    real(kind=r_def)    :: longitude, latitude, radius

    domain_depth = domain_top - domain_surface

    ipanel = int(panel_id(map_pid(1)), i_def)

    ! Calculate orography and update chi_3
    do df = 1, ndf_chi
      do k = 0, nlayers-1
        dfk = map_chi(df)+k

        ! Model coordinates for spherical domain are in (alpha,beta,r) form
        ! They need to be converted to (long,lat,r) for reading analytic orog
        radius = chi_3(dfk) + domain_surface
        call alphabetar2llr(chi_1(dfk), chi_2(dfk), radius, &
                            ipanel, longitude, latitude)

        ! Calculate surface height for each DoF using selected analytic orography
        surface_height = orography_profile%analytic_orography(longitude, latitude)

        ! Calculate nondimensional coordinate from current height coordinate
        ! (chi_3) with flat domain_surface
        eta = z2eta_linear(chi_3(dfk), 0.0_r_def, domain_depth)

        select case(stretching_method)
        case(stretching_method_linear)
          chi_3(dfk) = eta2z_linear(eta, surface_height, domain_depth)
        case default
          chi_3(dfk) = eta2z_smooth(eta, surface_height, domain_depth, stretching_height)
        end select

      end do
    end do

  end subroutine analytic_orography_spherical_alphabetaz

  !=============================================================================
  !> @brief Updates (longitude, latitude, h) vertical coordinates for a single
  !>        column using selected analytic orography.
  !>
  !> @details Calculates analytic orography from chi_1 and chi_2 horizontal
  !>          coordinates and then updates chi_3. This works directly on the
  !>          (longitude, latitude) coordinates.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in]     ndf_chi        Array size and loop bound for map_chi
  !> @param[in]     undf_chi       Column coordinates' array size and loop bound
  !> @param[in]     map_chi        Indirection map for coordinate field
  !> @param[in]     ndf_pid        Array size and loop bound for map_pid
  !> @param[in]     undf_pid       Panel ID array size and loop bound
  !> @param[in]     map_pid        Indirection map for panel_id
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in,out] chi_1          1st coordinate field in Wchi
  !> @param[in,out] chi_2          2nd coordinate field in Wchi
  !> @param[in,out] chi_3          3rd coordinate field in Wchi
  !> @param[in]     panel_id       Field giving the ID for mesh panels
  !=============================================================================
  subroutine analytic_orography_spherical_lonlatz(nlayers,        &
                                                  ndf_chi,        &
                                                  undf_chi,       &
                                                  map_chi,        &
                                                  ndf_pid,        &
                                                  undf_pid,       &
                                                  map_pid,        &
                                                  domain_surface, &
                                                  domain_top,     &
                                                  chi_1,          &
                                                  chi_2,          &
                                                  chi_3,          &
                                                  panel_id)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers, undf_chi, undf_pid
    integer(kind=i_def), intent(in)    :: ndf_chi, ndf_pid
    integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
    integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
    real(kind=r_def),    intent(in)    :: domain_surface, domain_top
    real(kind=r_def),    intent(inout) :: chi_1(undf_chi), chi_2(undf_chi), chi_3(undf_chi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    ! Internal variables
    integer(kind=i_def) :: k, df, dfk
    real(kind=r_def)    :: surface_height, domain_depth
    real(kind=r_def)    :: eta

    domain_depth = domain_top - domain_surface

    ! Calculate orography and update chi_3
    do df = 1, ndf_chi
      do k = 0, nlayers-1
        dfk = map_chi(df)+k

        ! Model coordinates for spherical domain are in (lon,lat,r) form
        ! Must be in (lon,lat,r) for reading analytic orography so don't convert
        ! Calculate surface height for each DoF using selected analytic orography
        surface_height = orography_profile%analytic_orography(chi_1(dfk), chi_2(dfk))

        ! Calculate nondimensional coordinate from current height coordinate
        ! (chi_3) with flat domain_surface
        eta = z2eta_linear(chi_3(dfk), 0.0_r_def, domain_depth)

        select case(stretching_method)
        case(stretching_method_linear)
          chi_3(dfk) = eta2z_linear(eta, surface_height, domain_depth)
        case default
          chi_3(dfk) = eta2z_smooth(eta, surface_height, domain_depth, stretching_height)
        end select

      end do
    end do

  end subroutine analytic_orography_spherical_lonlatz

  !=============================================================================
  !> @brief Updates Cartesian vertical coordinate for a single column using
  !>        selected analytic orography.
  !>
  !> @details Calculates analytic orography from chi_1 and chi_2 horizontal
  !>          coordinates and then updates chi_3. After evaluation of the new
  !>          surface height chi_3 is updated using its nondimensional eta
  !>          coordinate.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in]     ndf_chi        Array size and loop bound for map_chi
  !> @param[in]     undf_chi       Column coordinates' array size and loop bound
  !> @param[in]     map_chi        Indirection map for coordinate field
  !> @param[in]     ndf_pid        Array size and loop bound for map_pid
  !> @param[in]     undf_pid       Panel ID array size and loop bound
  !> @param[in]     map_pid        Indirection map for panel_id
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in,out] chi_1          1st coordinate field in Wchi
  !> @param[in,out] chi_2          2nd coordinate field in Wchi
  !> @param[in,out] chi_3          3rd coordinate field in Wchi
  !> @param[in]     panel_id       Field giving the ID for mesh panels
  !=============================================================================
  subroutine analytic_orography_cartesian(nlayers, ndf_chi, undf_chi, map_chi, &
                                          ndf_pid, undf_pid, map_pid,          &
                                          domain_surface, domain_top,          &
                                          chi_1, chi_2, chi_3, panel_id)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers, undf_chi, undf_pid
    integer(kind=i_def), intent(in)    :: ndf_chi, ndf_pid
    integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
    integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
    real(kind=r_def),    intent(in)    :: domain_surface, domain_top
    real(kind=r_def),    intent(inout) :: chi_1(undf_chi), chi_2(undf_chi), chi_3(undf_chi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: k, df, dfk
    real(kind=r_def)    :: surface_height, domain_depth
    real(kind=r_def)    :: eta

    domain_depth = domain_top - domain_surface

    ! Calculate orography and update chi_3
    do df = 1, ndf_chi
      do k = 0, nlayers-1
        dfk = map_chi(df)+k

        ! Calculate surface height for each DoF using selected analytic orography
        surface_height = orography_profile%analytic_orography(chi_1(dfk), chi_2(dfk))

        ! Calculate nondimensional coordinate from current height coordinate
        ! (chi_3) with flat domain_surface
        eta = z2eta_linear(chi_3(dfk), domain_surface, domain_top)

        select case(stretching_method)
        case(stretching_method_linear)
          chi_3(dfk) = eta2z_linear(eta, domain_surface + surface_height, domain_top)
        case default
          chi_3(dfk) = domain_surface + &
             eta2z_smooth(eta, surface_height, domain_depth, stretching_height)
        end select
      end do
    end do

    return
  end subroutine analytic_orography_cartesian

  !=============================================================================
  !> @brief Modify vertical coordinate based on the input surface_altitude field.
  !>        For spherical geometries with a Cartesian coordinate system.
  !> Note that this routine assumes the chi coordinates in a column are
  !> associated with a flat domain on input and then modified on output.
  !> Therefore it will not operate correctly with a horizontally continuous chi
  !> field.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in,out] chi_1          1st coordinate field in Wchi
  !> @param[in,out] chi_2          2nd coordinate field in Wchi
  !> @param[in,out] chi_3          3rd coordinate field in Wchi
  !> @param[in]     panel_id       Field giving the ID for mesh panels
  !> @param[in]     surface_altitude Surface altitude field data
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in]     ndf_chi        Array size and loop bound for map_chi
  !> @param[in]     undf_chi       Column coordinates' array size and loop bound
  !> @param[in]     map_chi        Indirection map for coordinate field
  !> @param[in]     ndf_pid        Array size and loop bound for map_pid
  !> @param[in]     undf_pid       Panel ID array size and loop bound
  !> @param[in]     map_pid        Indirection map for pid
  !> @param[in]     ndf            Array size and loop bound for surface altitude field
  !> @param[in]     undf           Total number of dofs for surface altitude field
  !> @param[in]     map            Indirection map for surface altitude field
  !> @param[in]     basis          Basis functions for surface altitude field
  !=============================================================================
  subroutine ancil_orography_spherical_xyz(nlayers,                    &
                                           chi_1, chi_2, chi_3,        &
                                           panel_id,                   &
                                           surface_altitude,           &
                                           domain_surface, domain_top, &
                                           ndf_chi, undf_chi,          &
                                           map_chi,                    &
                                           ndf_pid, undf_pid,          &
                                           map_pid,                    &
                                           ndf, undf,                  &
                                           map, basis                  &
                                          )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf, undf_chi, undf_pid

  integer(kind=i_def), dimension(ndf),      intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi),  intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),  intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(ndf, ndf_chi) :: basis

  real(kind=r_def), dimension(undf_chi), intent(inout) :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id
  real(kind=r_def), dimension(undf),     intent(in)    :: surface_altitude
  real(kind=r_def), intent(in)                         :: domain_surface, domain_top

  ! Internal variables
  integer(kind=i_def) :: k, df, dfchi, dfk
  real(kind=r_def)    :: chi_3_r, eta
  real(kind=r_def)    :: surface_height(ndf_chi), domain_depth
  real(kind=r_def)    :: longitude, latitude, r

  domain_depth = domain_top - domain_surface

  ! Calculate new surface_height at each chi dof
  surface_height(:) = 0.0_r_def
  do dfchi = 1, ndf_chi
    do df = 1, ndf
      surface_height(dfchi) = surface_height(dfchi) + surface_altitude(map(df))*basis(df,dfchi)
    end do
  end do

  ! Update chi
  do df = 1, ndf_chi
    do k = 0, nlayers-1
      dfk = map_chi(df)+k

      ! Model coordinates for spherical domain are in (x,y,z) form so they need
      ! to be converted to (long,lat,r) first
      call xyz2llr(chi_1(dfk), chi_2(dfk), chi_3(dfk), longitude, latitude, r)

      ! Calculate nondimensional coordinate from current flat height coordinate
      ! (chi_3) with flat domain_surface
      eta = z2eta_linear(r, domain_surface, domain_top)

      ! Calculate new height spherical coordinate (chi_3_r) from its
      ! nondimensional coordinate eta and surface_height
      select case(stretching_method)
      case(stretching_method_linear)
        chi_3_r = eta2z_linear(eta, domain_surface+surface_height(df), domain_top)
      case default
        chi_3_r = domain_surface + &
           eta2z_smooth(eta, surface_height(df), domain_depth, stretching_height)
      end select

      ! Convert spherical coordinates back to model (x,y,z) form
      call llr2xyz(longitude, latitude, chi_3_r, &
                   chi_1(dfk), chi_2(dfk), chi_3(dfk))
    end do
  end do

  end subroutine ancil_orography_spherical_xyz

  !=============================================================================
  !> @brief Modify vertical coordinate based on the input surface_altitude field.
  !>        For spherical geometries with (alpha,beta,r) or (lon,lat,r)
  !>        coordinate systems.
  !> Note that this routine assumes the chi coordinates in a column are
  !> associated with a flat domain on input and then modified on output.
  !> Therefore it will not operate correctly with a horizontally continuous chi
  !> field.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in,out] chi_1          1st coordinate field in Wchi
  !> @param[in,out] chi_2          2nd coordinate field in Wchi
  !> @param[in,out] chi_3          3rd coordinate field in Wchi
  !> @param[in]     panel_id       Field giving the ID for mesh panels
  !> @param[in]     surface_altitude Surface altitude field data
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in]     ndf_chi        Array size and loop bound for map_chi
  !> @param[in]     undf_chi       Column coordinates' array size and loop bound
  !> @param[in]     map_chi        Indirection map for coordinate field
  !> @param[in]     ndf_pid        Array size and loop bound for map_pid
  !> @param[in]     undf_pid       Panel ID array size and loop bound
  !> @param[in]     map_pid        Indirection map for panel ID
  !> @param[in]     ndf            Array size and loop bound for surface altitude field
  !> @param[in]     undf           Total number of dofs for surface altitude field
  !> @param[in]     map            Indirection map for surface altitude field
  !> @param[in]     basis          Basis functions for surface altitude field
  !=============================================================================
  subroutine ancil_orography_spherical_sph(nlayers,                    &
                                           chi_1, chi_2, chi_3,        &
                                           panel_id,                   &
                                           surface_altitude,           &
                                           domain_surface, domain_top, &
                                           ndf_chi, undf_chi,          &
                                           map_chi,                    &
                                           ndf_pid, undf_pid,          &
                                           map_pid,                    &
                                           ndf, undf,                  &
                                           map, basis                  &
                                          )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf, undf_chi, undf_pid

  integer(kind=i_def), dimension(ndf),      intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi),  intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),  intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(ndf, ndf_chi) :: basis

  real(kind=r_def), dimension(undf_chi), intent(inout) :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid),  intent(in)   :: panel_id
  real(kind=r_def), dimension(undf),     intent(in)    :: surface_altitude
  real(kind=r_def), intent(in)                         :: domain_surface, domain_top

  ! Internal variables
  integer(kind=i_def) :: k, df, dfchi, dfk
  real(kind=r_def)    :: eta
  real(kind=r_def)    :: surface_height(ndf_chi), domain_depth

  domain_depth = domain_top - domain_surface

  ! Calculate new surface_height at each chi dof
  surface_height(:) = 0.0_r_def
  do dfchi = 1, ndf_chi
    do df = 1, ndf
      surface_height(dfchi) = surface_height(dfchi) + surface_altitude(map(df))*basis(df,dfchi)
    end do
  end do

  ! Update chi
  do df = 1, ndf_chi
    do k = 0, nlayers-1
      dfk = map_chi(df)+k

      ! Calculate nondimensional coordinate from current flat height coordinate
      ! (chi_3) with flat domain_surface
      eta = z2eta_linear(chi_3(dfk), 0.0_r_def, domain_depth)

      ! Calculate new height coordinate from its nondimensional coordinate
      ! eta and surface_height
      select case(stretching_method)
      case(stretching_method_linear)
        chi_3(dfk) = eta2z_linear(eta, surface_height(df), domain_depth)
      case default
        chi_3(dfk) = eta2z_smooth(eta, surface_height(df), domain_depth, stretching_height)
      end select
    end do
  end do

end subroutine ancil_orography_spherical_sph

  !=============================================================================
  !> @brief Modify vertical coordinate based on the input surface_altitude field.
  !>        For Cartesian geometries.
  !> Note that this routine assumes the chi coordinates in a column are
  !> associated with a flat domain on input and then modified on output.
  !> Therefore it will not operate correctly with a horizontally continuous chi
  !> field.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in,out] chi_1          1st coordinate field in Wchi
  !> @param[in,out] chi_2          2nd coordinate field in Wchi
  !> @param[in,out] chi_3          3rd coordinate field in Wchi
  !> @param[in]     panel_id       Field giving the ID for mesh panels
  !> @param[in]     surface_altitude Surface altitude field data
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in]     ndf_chi        Array size and loop bound for map_chi
  !> @param[in]     undf_chi       Column coordinates' array size and loop bound
  !> @param[in]     map_chi        Indirection map for coordinate field
  !> @param[in]     ndf_pid        Array size and loop bound for map_pid
  !> @param[in]     undf_pid       Panel ID array size and loop bound
  !> @param[in]     map_pid        Indirection map for panel_id
  !> @param[in]     ndf            Array size and loop bound for surface altitude field
  !> @param[in]     undf           Total number of dofs for surface altitude field
  !> @param[in]     map            Indirection map for surface altitude field
  !> @param[in]     basis          Basis functions for surface altitude field
  !=============================================================================
  subroutine ancil_orography_cartesian(nlayers,                    &
                                       chi_1, chi_2, chi_3,        &
                                       panel_id,                   &
                                       surface_altitude,           &
                                       domain_surface, domain_top, &
                                       ndf_chi, undf_chi,          &
                                       map_chi,                    &
                                       ndf_pid, undf_pid,          &
                                       map_pid,                    &
                                       ndf, undf,                  &
                                       map, basis                  &
                                       )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf, undf_chi, undf_pid

  integer(kind=i_def), dimension(ndf),      intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi),  intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),  intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(ndf, ndf_chi) :: basis

  real(kind=r_def), dimension(undf_chi), intent(inout) :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid),  intent(in)   :: panel_id
  real(kind=r_def), dimension(undf),     intent(in)    :: surface_altitude
  real(kind=r_def), intent(in)                         :: domain_surface, domain_top

  ! Internal variables
  integer(kind=i_def) :: k, df, dfchi, dfk
  real(kind=r_def)    :: eta
  real(kind=r_def)    :: surface_height(ndf_chi), domain_depth

  domain_depth = domain_top - domain_surface

  ! Calculate new surface_height at each chi dof
  surface_height(:) = 0.0_r_def
  do dfchi = 1, ndf_chi
    do df = 1, ndf
      surface_height(dfchi) = surface_height(dfchi) + surface_altitude(map(df))*basis(df,dfchi)
    end do
  end do

  ! Update chi_3
  do df = 1, ndf_chi
    do k = 0, nlayers-1
      dfk = map_chi(df)+k

      ! Calculate nondimensional coordinate from current height coordinate
      ! (chi_3) with flat domain_surface
      eta = z2eta_linear(chi_3(dfk), domain_surface, domain_top)

      ! Calculate new height coordinate from its nondimensional coordinate
      ! eta and surface_height
      select case(stretching_method)
      case(stretching_method_linear)
        chi_3(dfk) = eta2z_linear(eta, domain_surface+surface_height(df), domain_top)
      case default
        chi_3(dfk) = domain_surface + &
           eta2z_smooth(eta, surface_height(df), domain_depth, stretching_height)
      end select
    end do
  end do

end subroutine ancil_orography_cartesian

end module assign_orography_field_mod
