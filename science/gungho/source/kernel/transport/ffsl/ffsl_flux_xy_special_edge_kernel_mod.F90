!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Calculates the horizontal mass flux for FFSL, with special treatment
!!       at edges of cubed-sphere panels.
!> @details This kernel computes the flux in the x and y directions. A choice of
!!          constant, Nirvana or PPM is used to compute the edge reconstruction.
!!          This is multiplied by the fractional velocity to give
!!          the flux. For CFL > 1 the field values multiplied by the mass
!!          are summed between the flux point and the departure cell.
!!          Special treatment is used at the edges of cubed-sphere panels, to
!!          shift the reconstruction to improve accuracy.
!!          This kernel is used for the FFSL horizontal transport scheme.
!!
!> @note This kernel only works when field is a W3 field at lowest order since
!!       it is assumed that ndf_w3 = 1 with stencil_map(1,:) containing the
!!       relevant dofmaps.

module ffsl_flux_xy_special_edge_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    CELL_COLUMN, GH_WRITE,     &
                                    GH_READ, GH_SCALAR,        &
                                    STENCIL, GH_INTEGER,       &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    ANY_DISCONTINUOUS_SPACE_4, &
                                    ANY_DISCONTINUOUS_SPACE_5, &
                                    CROSS2D, GH_LOGICAL
  use constants_mod,         only : i_def, r_tran, r_def, l_def
  use fs_continuity_mod,     only : W3
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : W, E, N, S

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: ffsl_flux_xy_special_edge_kernel_type
    private
    type(arg_type) :: meta_args(21) = (/                                       &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! flux
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! field_for_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! dry_mass_for_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! field_for_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! dry_mass_for_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! frac_dry_flux
        arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! dep_lowest_k
        arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! dep_highest_k
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! panel_id
        arg_type(GH_FIELD*2, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! panel_edge_dist_x
        arg_type(GH_FIELD*2, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! panel_edge_dist_y
        arg_type(GH_FIELD,  GH_INTEGER,  GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ew
        arg_type(GH_FIELD,  GH_INTEGER,  GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ns
        arg_type(GH_SCALAR, GH_LOGICAL,  GH_READ),                             & ! high_order_edges
        arg_type(GH_SCALAR, GH_INTEGER,  GH_READ),                             & ! order
        arg_type(GH_SCALAR, GH_INTEGER,  GH_READ),                             & ! monotone
        arg_type(GH_SCALAR, GH_INTEGER,  GH_READ),                             & ! edge_monotone
        arg_type(GH_SCALAR, GH_REAL,     GH_READ),                             & ! min_val
        arg_type(GH_SCALAR, GH_INTEGER,  GH_READ),                             & ! ndep
        arg_type(GH_SCALAR, GH_REAL,     GH_READ)                              & ! dt
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_flux_xy_special_edge_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public  :: ffsl_flux_xy_special_edge_code
  private :: ffsl_flux_xy_special_edge_1d
contains

  !> @brief Compute the horizontal fluxes for FFSL.
  !> @param[in]     nlayers             Number of layers
  !> @param[in,out] flux                The output flux
  !> @param[in]     field_for_x         Field to use in evaluating x-flux
  !> @param[in]     stencil_sizes_x     Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_x       Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_x       Map of DoFs in the stencil for x-field
  !> @param[in]     dry_mass_for_x      Volume or dry mass field at W3 points
  !!                                    for use in evaluating x-flux
  !> @param[in]     stencil_sizes_mx    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_mx      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_mx      Map of DoFs in the stencil for x-mass
  !> @param[in]     field_for_y         Field to use in evaluating x-flux
  !> @param[in]     stencil_sizes_y     Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_y       Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_y       Map of DoFs in the stencil for y-field
  !> @param[in]     dry_mass_for_y      Volume or dry mass field at W3 points
  !!                                    for use in evaluating x-flux
  !> @param[in]     stencil_sizes_my    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_my      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_my      Map of DoFs in the stencil for y-mass
  !> @param[in]     dep_dist            Horizontal departure distances
  !> @param[in]     frac_dry_flux       Fractional part of the dry flux or wind
  !> @param[in]     dep_lowest_k        2D integer multidata W2H field, storing
  !!                                    the lowest model level to use in each
  !!                                    integer flux sum, for each column
  !> @param[in]     dep_highest_k       2D integer multidata W2H field, storing
  !!                                    the highest model level to use in each
  !!                                    integer flux sum, for each column
  !> @param[in]     panel_id            Field containing IDs of mesh panels
  !> @param[in]     panel_edge_dist_W   2D field containing the distance of each
  !!                                    column from the panel edge to the West
  !> @param[in]     panel_edge_dist_E   2D field containing the distance of each
  !!                                    column from the panel edge to the East
  !> @param[in]     stencil_sizes_cx    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_cx      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_cx      Map of DoFs in the stencil for x-case
  !> @param[in]     panel_edge_dist_S   2D field containing the distance of each
  !!                                    column from the panel edge to the South
  !> @param[in]     panel_edge_dist_N   2D field containing the distance of each
  !!                                    column from the panel edge to the North
  !> @param[in]     stencil_sizes_cy    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_cy      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_cy      Map of DoFs in the stencil for y-case
  !> @param[in]     face_selector_ew    2D field indicating which W/E faces to
  !!                                    loop over for this column
  !> @param[in]     face_selector_ns    2D field indicating which N/S faces to
  !!                                    loop over for this column
  !> @param[in]     high_order_edges    Flag to use high order reconstruction
  !!                                    at panel edges
  !> @param[in]     order               Order of reconstruction
  !> @param[in]     monotone            Horizontal monotone option for FFSL
  !> @param[in]     edge_monotone       Monotone option for FFSL at panel edges
  !> @param[in]     min_val             Minimum value to enforce when using
  !!                                    quasi-monotone limiter
  !> @param[in]     ndep                Number of multidata points for departure
  !!                                    index fields
  !> @param[in]     dt                  Time step
  !> @param[in]     ndf_w2h             Num of DoFs for W2h per cell
  !> @param[in]     undf_w2h            Num of DoFs for W2h in this partition
  !> @param[in]     map_w2h             Map for W2h
  !> @param[in]     ndf_w3              Num of DoFs for W3 per cell
  !> @param[in]     undf_w3             Num of DoFs for W3 in this partition
  !> @param[in]     map_w3              Map for W3
  !> @param[in]     ndf_depk            Num of DoFs for dep idx fields per cell
  !> @param[in]     undf_depk           Num of DoFs for this partition for dep
  !!                                    idx fields
  !> @param[in]     map_depk            Map for departure index fields
  !> @param[in]     ndf_pid             Num of DoFs for panel ID per cell
  !> @param[in]     undf_pid            Num DoFs for this partition for panel ID
  !> @param[in]     map_pid             Map for panel ID
  !> @param[in]     ndf_w3_2d           Num of DoFs for 2D W3 per cell
  !> @param[in]     undf_w3_2d          Num of DoFs for this partition for 2D W3
  !> @param[in]     map_w3_2d           Map for 2D W3
  subroutine ffsl_flux_xy_special_edge_code( nlayers,             &
                                             flux,                &
                                             field_for_x,         &
                                             stencil_sizes_x,     &
                                             stencil_max_x,       &
                                             stencil_map_x,       &
                                             dry_mass_for_x,      &
                                             stencil_sizes_mx,    &
                                             stencil_max_mx,      &
                                             stencil_map_mx,      &
                                             field_for_y,         &
                                             stencil_sizes_y,     &
                                             stencil_max_y,       &
                                             stencil_map_y,       &
                                             dry_mass_for_y,      &
                                             stencil_sizes_my,    &
                                             stencil_max_my,      &
                                             stencil_map_my,      &
                                             dep_dist,            &
                                             frac_dry_flux,       &
                                             dep_lowest_k,        &
                                             dep_highest_k,       &
                                             panel_id,            &
                                             panel_edge_dist_W,   &
                                             panel_edge_dist_E,   &
                                             stencil_sizes_cx,    &
                                             stencil_max_cx,      &
                                             stencil_map_cx,      &
                                             panel_edge_dist_S,   &
                                             panel_edge_dist_N,   &
                                             stencil_sizes_cy,    &
                                             stencil_max_cy,      &
                                             stencil_map_cy,      &
                                             face_selector_ew,    &
                                             face_selector_ns,    &
                                             high_order_edges,    &
                                             order,               &
                                             monotone,            &
                                             edge_monotone,       &
                                             min_val,             &
                                             ndep,                &
                                             dt,                  &
                                             ndf_w2h,             &
                                             undf_w2h,            &
                                             map_w2h,             &
                                             ndf_w3,              &
                                             undf_w3,             &
                                             map_w3,              &
                                             ndf_depk,            &
                                             undf_depk,           &
                                             map_depk,            &
                                             ndf_pid,             &
                                             undf_pid,            &
                                             map_pid,             &
                                             ndf_w3_2d,           &
                                             undf_w3_2d,          &
                                             map_w3_2d )

    use ffsl_flux_xy_kernel_mod, only: ffsl_flux_xy_1d
    use panel_edge_support_mod,  only: crosses_panel_edge

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_pid
    integer(kind=i_def), intent(in) :: undf_pid
    integer(kind=i_def), intent(in) :: ndf_depk
    integer(kind=i_def), intent(in) :: undf_depk
    integer(kind=i_def), intent(in) :: ndep
    integer(kind=i_def), intent(in) :: stencil_max_x
    integer(kind=i_def), intent(in) :: stencil_max_mx
    integer(kind=i_def), intent(in) :: stencil_max_y
    integer(kind=i_def), intent(in) :: stencil_max_my
    integer(kind=i_def), intent(in) :: stencil_max_cx
    integer(kind=i_def), intent(in) :: stencil_max_cy
    integer(kind=i_def), intent(in) :: stencil_sizes_x(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_mx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_y(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_my(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_cx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_cy(4)

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
    integer(kind=i_def), intent(in) :: map_depk(ndf_depk)
    integer(kind=i_def), intent(in) :: stencil_map_x(ndf_w3,stencil_max_x,4)
    integer(kind=i_def), intent(in) :: stencil_map_mx(ndf_w3,stencil_max_mx,4)
    integer(kind=i_def), intent(in) :: stencil_map_y(ndf_w3,stencil_max_y,4)
    integer(kind=i_def), intent(in) :: stencil_map_my(ndf_w3,stencil_max_my,4)
    integer(kind=i_def), intent(in) :: stencil_map_cx(ndf_pid,stencil_max_cx,4)
    integer(kind=i_def), intent(in) :: stencil_map_cy(ndf_pid,stencil_max_cy,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: field_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: field_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)
    real(kind=r_tran),   intent(in)    :: frac_dry_flux(undf_w2h)
    logical(kind=l_def), intent(in)    :: high_order_edges
    integer(kind=i_def), intent(in)    :: order
    integer(kind=i_def), intent(in)    :: monotone
    integer(kind=i_def), intent(in)    :: edge_monotone
    real(kind=r_tran),   intent(in)    :: min_val
    real(kind=r_tran),   intent(in)    :: dt
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
    integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: dep_lowest_k(undf_depk)
    integer(kind=i_def), intent(in)    :: dep_highest_k(undf_depk)

    ! Internal arguments
    integer(kind=i_def) :: i, ipanel, recon_size
    integer(kind=i_def) :: edge_dist_E, edge_dist_W, edge_dist_S, edge_dist_N
    logical(kind=l_def) :: use_edge_treatment
    integer(kind=i_def) :: stencil_extent_xl, stencil_extent_xr
    integer(kind=i_def) :: stencil_extent_yl, stencil_extent_yr
    integer(kind=i_def) :: stencil_map_x_1d(-stencil_max_x:stencil_max_x)
    integer(kind=i_def) :: stencil_map_y_1d(-stencil_max_y:stencil_max_y)
    integer(kind=i_def) :: stencil_map_x_pid(-stencil_max_x:stencil_max_x)
    integer(kind=i_def) :: stencil_map_y_pid(-stencil_max_y:stencil_max_y)

    ! Set size of array for reconstruction
    if (high_order_edges) then
      ! Increase stencil size by one for special edges
      recon_size  = 3 + 2*order
    else
      ! Use same stencil size as without special edges
      recon_size  = 1 + 2*order
    end if

    ! Form X and Y 1D stencils
    stencil_extent_xl = stencil_sizes_x(1) - 1
    stencil_extent_xr = stencil_sizes_x(3) - 1
    stencil_extent_yl = stencil_sizes_y(2) - 1
    stencil_extent_yr = stencil_sizes_y(4) - 1

    do i = -stencil_extent_xl, 0
      stencil_map_x_1d(i) = stencil_map_x(1, 1-i, 1)
      stencil_map_x_pid(i) = stencil_map_cx(1, 1-i, 1)
    end do
    do i = 1, stencil_extent_xr
      stencil_map_x_1d(i) = stencil_map_x(1, i+1, 3)
      stencil_map_x_pid(i) = stencil_map_cx(1, i+1, 3)
    end do

    do i = -stencil_extent_yl, 0
      stencil_map_y_1d(i) = stencil_map_y(1, 1-i, 2)
      stencil_map_y_pid(i) = stencil_map_cy(1, 1-i, 2)
    end do
    do i = 1, stencil_extent_yr
      stencil_map_y_1d(i) = stencil_map_y(1, i+1, 4)
      stencil_map_y_pid(i) = stencil_map_cy(1, i+1, 4)
    end do

    ipanel = INT(panel_id(map_pid(1)), i_def)
    edge_dist_E = panel_edge_dist_E(map_pid(1))
    edge_dist_W = panel_edge_dist_W(map_pid(1))
    edge_dist_S = panel_edge_dist_S(map_pid(1))
    edge_dist_N = panel_edge_dist_N(map_pid(1))

    ! X direction ==============================================================
    use_edge_treatment = crosses_panel_edge(                                   &
        edge_dist_W, edge_dist_E, MAX(stencil_extent_xl, stencil_extent_xr),   &
        order, face_selector_ew(map_w3_2d(1)), .false., ipanel, 1, W, E,       &
        dep_highest_k, ndep, ndf_depk, undf_depk, map_depk                     &
    )

    if (use_edge_treatment) then
      call ffsl_flux_xy_special_edge_1d( nlayers,             &
                                         .true.,              &
                                         flux,                &
                                         field_for_x,         &
                                         field_for_y,         &
                                         stencil_extent_xl,   &
                                         stencil_extent_xr,   &
                                         stencil_max_x,       &
                                         stencil_map_x_1d,    &
                                         dry_mass_for_x,      &
                                         dry_mass_for_y,      &
                                         dep_dist,            &
                                         frac_dry_flux,       &
                                         dep_lowest_k,        &
                                         dep_highest_k,       &
                                         ipanel,              &
                                         panel_edge_dist_W,   &
                                         panel_edge_dist_E,   &
                                         panel_edge_dist_S,   &
                                         panel_edge_dist_N,   &
                                         stencil_map_x_pid,   &
                                         face_selector_ew,    &
                                         high_order_edges,    &
                                         recon_size,          &
                                         order,               &
                                         monotone,            &
                                         edge_monotone,       &
                                         min_val,             &
                                         ndep,                &
                                         dt,                  &
                                         ndf_w2h,             &
                                         undf_w2h,            &
                                         map_w2h,             &
                                         ndf_w3,              &
                                         undf_w3,             &
                                         map_w3,              &
                                         ndf_depk,            &
                                         undf_depk,           &
                                         map_depk,            &
                                         ndf_pid,             &
                                         undf_pid,            &
                                         map_pid,             &
                                         ndf_w3_2d,           &
                                         undf_w3_2d,          &
                                         map_w3_2d )
    else
      ! Not near a panel edge so default to standard FFSL
      call ffsl_flux_xy_1d( nlayers,             &
                            .true.,              &
                            flux,                &
                            field_for_x,         &
                            stencil_extent_xl,   &
                            stencil_extent_xr,   &
                            stencil_max_x,       &
                            stencil_map_x_1d,    &
                            dry_mass_for_x,      &
                            dep_dist,            &
                            frac_dry_flux,       &
                            dep_lowest_k,        &
                            dep_highest_k,       &
                            face_selector_ew,    &
                            order,               &
                            monotone,            &
                            min_val,             &
                            ndep,                &
                            dt,                  &
                            ndf_w2h,             &
                            undf_w2h,            &
                            map_w2h,             &
                            ndf_w3,              &
                            undf_w3,             &
                            map_w3,              &
                            ndf_depk,            &
                            undf_depk,           &
                            map_depk,            &
                            ndf_w3_2d,           &
                            undf_w3_2d,          &
                            map_w3_2d )
    end if

    ! Y direction ==============================================================
    use_edge_treatment = crosses_panel_edge(                                   &
        edge_dist_S, edge_dist_N, MAX(stencil_extent_yl, stencil_extent_yr),   &
        order, face_selector_ns(map_w3_2d(1)), .false., ipanel, 2, S, N,       &
        dep_highest_k, ndep, ndf_depk, undf_depk, map_depk                     &
    )

    if (use_edge_treatment) then
      call ffsl_flux_xy_special_edge_1d( nlayers,             &
                                         .false.,             &
                                         flux,                &
                                         field_for_y,         &
                                         field_for_x,         &
                                         stencil_extent_yl,   &
                                         stencil_extent_yr,   &
                                         stencil_max_y,       &
                                         stencil_map_y_1d,    &
                                         dry_mass_for_y,      &
                                         dry_mass_for_x,      &
                                         dep_dist,            &
                                         frac_dry_flux,       &
                                         dep_lowest_k,        &
                                         dep_highest_k,       &
                                         ipanel,              &
                                         panel_edge_dist_S,   &
                                         panel_edge_dist_N,   &
                                         ! NB: swap E and W here to simplify
                                         !     the logic in the 1D routine
                                         panel_edge_dist_E,   &
                                         panel_edge_dist_W,   &
                                         stencil_map_y_pid,   &
                                         face_selector_ns,    &
                                         high_order_edges,    &
                                         recon_size,          &
                                         order,               &
                                         monotone,            &
                                         edge_monotone,       &
                                         min_val,             &
                                         ndep,                &
                                         dt,                  &
                                         ndf_w2h,             &
                                         undf_w2h,            &
                                         map_w2h,             &
                                         ndf_w3,              &
                                         undf_w3,             &
                                         map_w3,              &
                                         ndf_depk,            &
                                         undf_depk,           &
                                         map_depk,            &
                                         ndf_pid,             &
                                         undf_pid,            &
                                         map_pid,             &
                                         ndf_w3_2d,           &
                                         undf_w3_2d,          &
                                         map_w3_2d )
    else
      ! Not near a panel edge so default to standard FFSL
      call ffsl_flux_xy_1d( nlayers,             &
                            .false.,             &
                            flux,                &
                            field_for_y,         &
                            stencil_extent_yl,   &
                            stencil_extent_yr,   &
                            stencil_max_y,       &
                            stencil_map_y_1d,    &
                            dry_mass_for_y,      &
                            dep_dist,            &
                            frac_dry_flux,       &
                            dep_lowest_k,        &
                            dep_highest_k,       &
                            face_selector_ns,    &
                            order,               &
                            monotone,            &
                            min_val,             &
                            ndep,                &
                            dt,                  &
                            ndf_w2h,             &
                            undf_w2h,            &
                            map_w2h,             &
                            ndf_w3,              &
                            undf_w3,             &
                            map_w3,              &
                            ndf_depk,            &
                            undf_depk,           &
                            map_depk,            &
                            ndf_w3_2d,           &
                            undf_w3_2d,          &
                            map_w3_2d )
    end if

  end subroutine ffsl_flux_xy_special_edge_code

! ============================================================================ !
! SINGLE UNDERLYING 1D ROUTINE
! ============================================================================ !

  !> @brief Compute the flux with shifted edge stencils for either the
  !!        x- or y-direction
  subroutine ffsl_flux_xy_special_edge_1d( nlayers,                   &
                                           x_direction,               &
                                           flux,                      &
                                           field,                     &
                                           field_swapped,             &
                                           stencil_extent_l,          &
                                           stencil_extent_r,          &
                                           stencil_max,               &
                                           stencil_map,               &
                                           dry_mass,                  &
                                           dry_mass_swapped,          &
                                           dep_dist,                  &
                                           frac_dry_flux,             &
                                           dep_lowest_k,              &
                                           dep_highest_k,             &
                                           ipanel,                    &
                                           panel_edge_dist_l,         &
                                           panel_edge_dist_r,         &
                                           panel_edge_dist_swapped_l, &
                                           panel_edge_dist_swapped_r, &
                                           stencil_map_pid,           &
                                           face_selector,             &
                                           high_order_edges,          &
                                           recon_size,                &
                                           order,                     &
                                           monotone,                  &
                                           edge_monotone,             &
                                           min_val,                   &
                                           ndep,                      &
                                           dt,                        &
                                           ndf_w2h,                   &
                                           undf_w2h,                  &
                                           map_w2h,                   &
                                           ndf_w3,                    &
                                           undf_w3,                   &
                                           map_w3,                    &
                                           ndf_depk,                  &
                                           undf_depk,                 &
                                           map_depk,                  &
                                           ndf_pid,                   &
                                           undf_pid,                  &
                                           map_pid,                   &
                                           ndf_w3_2d,                 &
                                           undf_w3_2d,                &
                                           map_w3_2d )

    use panel_edge_support_mod,         only: rotated_panel_neighbour, FAR_AWAY
    use subgrid_horizontal_support_mod, only: nirvana_special_edge,            &
                                              linear_special_edge,             &
                                              fourth_order_special_edge,       &
                                              fourth_nirvana_special_edge
    use subgrid_common_support_mod,     only: monotonic_edge,                  &
                                              subgrid_quadratic_recon

    implicit none

    ! Arguments: Function space metadata
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_pid
    integer(kind=i_def), intent(in) :: undf_pid
    integer(kind=i_def), intent(in) :: ndf_depk
    integer(kind=i_def), intent(in) :: undf_depk
    integer(kind=i_def), intent(in) :: ndep
    integer(kind=i_def), intent(in) :: stencil_extent_l
    integer(kind=i_def), intent(in) :: stencil_extent_r
    integer(kind=i_def), intent(in) :: stencil_max
    integer(kind=i_def), intent(in) :: ipanel
    integer(kind=i_def), intent(in) :: recon_size

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
    integer(kind=i_def), intent(in) :: map_depk(ndf_depk)
    integer(kind=i_def), intent(in) :: stencil_map(-stencil_max:stencil_max)
    integer(kind=i_def), intent(in) :: stencil_map_pid(-stencil_max:stencil_max)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)
    real(kind=r_tran),   intent(in)    :: frac_dry_flux(undf_w2h)
    integer(kind=i_def), intent(in)    :: face_selector(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: dep_lowest_k(undf_depk)
    integer(kind=i_def), intent(in)    :: dep_highest_k(undf_depk)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_l(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_r(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_swapped_l(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_swapped_r(undf_pid)

    ! Arguments: Targets
    real(kind=r_tran), target, intent(in) :: field(undf_w3)
    real(kind=r_tran), target, intent(in) :: field_swapped(undf_w3)
    real(kind=r_tran), target, intent(in) :: dry_mass(undf_w3)
    real(kind=r_tran), target, intent(in) :: dry_mass_swapped(undf_w3)

    ! Arguments: Scalars
    logical(kind=l_def), intent(in) :: high_order_edges
    integer(kind=i_def), intent(in) :: order
    integer(kind=i_def), intent(in) :: monotone
    integer(kind=i_def), intent(in) :: edge_monotone
    real(kind=r_tran),   intent(in) :: min_val
    real(kind=r_tran),   intent(in) :: dt
    logical(kind=l_def), intent(in) :: x_direction

    ! Local arrays
    integer(kind=i_def) :: int_cell_idx(nlayers)
    integer(kind=i_def) :: int_cell_switch(nlayers)
    integer(kind=i_def) :: sign_disp(nlayers)
    integer(kind=i_def) :: rel_dep_cell_idx(nlayers)
    integer(kind=i_def) :: rel_dep_cell_idx_bounded(nlayers)
    integer(kind=i_def) :: special_edge_case(nlayers)
    real(kind=r_tran)   :: displacement(nlayers)
    real(kind=r_tran)   :: frac_dist(nlayers)
    real(kind=r_tran)   :: recon_field(nlayers)
    real(kind=r_tran)   :: field_edge_left(nlayers)
    real(kind=r_tran)   :: field_edge_right(nlayers)
    real(kind=r_tran)   :: swap_switch(nlayers)
    real(kind=r_tran)   :: field_local(nlayers,recon_size)
    real(kind=r_tran)   :: field_local_tmp(1,recon_size)
    integer(kind=i_def) :: edge_dist_local(-stencil_max:stencil_max)
    integer(kind=i_def) :: mono_option(nlayers)
    integer(kind=i_def) :: rel_idx(nlayers)
    integer(kind=i_def) :: loc_idx(nlayers)
    integer(kind=i_def) :: edge_dist(nlayers)
    integer(kind=i_def) :: local_dofs(2)

    ! Local scalars
    integer(kind=i_def) :: dof_iterator
    integer(kind=i_def) :: col_idx, pid_idx, rel_idx_k
    integer(kind=i_def) :: dof_offset
    integer(kind=i_def) :: k, j, w2h_idx, idx_dep
    integer(kind=i_def) :: k_low, k_high, ndep_half
    integer(kind=i_def) :: rel_idx_min, rel_idx_max
    integer(kind=i_def) :: order_offset
    integer(kind=i_def) :: rot_edge_dist_l, rot_edge_dist_r
    integer(kind=i_def) :: edge_dist_l, edge_dist_r
    integer(kind=i_def) :: rotated_panel_l, rotated_panel_r
    logical(kind=l_def) :: monotone_same
    real(kind=r_tran)   :: direction
    real(kind=r_tran)   :: inv_dt

    real(kind=r_tran),   pointer :: field_ptr(:)
    real(kind=r_tran),   pointer :: mass_ptr(:)

    inv_dt = 1.0_r_tran/dt

    if (x_direction) then
      local_dofs = (/ W, E /)
      direction = 1.0_r_tran
    else
      ! y-direction
      local_dofs = (/ S, N /)
      direction = -1.0_r_tran
    end if

    ! Set stencil info ---------------------------------------------------------
    if (high_order_edges) then
      ! The order offset determines which 1D scheme to use
      ! and = 1 indicates higher order edges
      order_offset = 1
    else
      ! The order offset determines which 1D scheme to use
      order_offset = 0
    end if

    ndep_half = (ndep - 1_i_def) / 2_i_def
    rel_idx_min = MAX(-stencil_extent_l, -ndep_half)
    rel_idx_max = MIN(stencil_extent_r, ndep_half)

    ! Pre-determine aspects of crossing panel boundary -------------------------
    pid_idx = map_pid(1)
    rotated_panel_l = rotated_panel_neighbour(ipanel, local_dofs(1))
    rotated_panel_r = rotated_panel_neighbour(ipanel, local_dofs(2))

    ! Set distance on each side to a panel edge.
    ! If there is no rotation in panel, set the distance to be very large
    if (ABS(rotated_panel_l) > 0) then
      rot_edge_dist_l = -ABS(panel_edge_dist_l(pid_idx))
    else
      rot_edge_dist_l = -FAR_AWAY
    end if
    if (ABS(rotated_panel_r) > 0) then
      rot_edge_dist_r = ABS(panel_edge_dist_r(pid_idx))
    else
      rot_edge_dist_r = FAR_AWAY
    end if

    ! Fill a local array which contains the distance to the nearest panel edge
    ! for each column in the stencil, in the direction of the stencil. The
    ! problem is that there may be different nearest panel edges for each side
    ! of the stencil, so we need to create a unified array for each side

    ! This table gives the edge_dist field to use when crossing a panel edge
    ! Direction | Stencil Side | No Rotation | Clockwise | Anti-clockwise
    ! -------------------------------------------------------------------
    !     X     |      L       |      W      |     S     |       N
    !     X     |      R       |      E      |     N     |       S
    !     Y     |      L       |      S      |     E     |       W
    !     Y     |      R       |      N      |     W     |       E

    ! Loop over the left side of the stencil -----------------------------------
    do rel_idx_k = -stencil_extent_l, -1
      pid_idx = stencil_map_pid(rel_idx_k)
      if (rotated_panel_l == 1 .and. rel_idx_k < rot_edge_dist_l + 1) then
        ! Clockwise rotation to the left
        edge_dist_l = -ABS(panel_edge_dist_swapped_l(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_swapped_r(pid_idx))
      else if (rotated_panel_l == -1 .and. rel_idx_k < rot_edge_dist_l + 1) then
        ! Anti-clockwise rotation to the left, so swap the directions
        edge_dist_l = -ABS(panel_edge_dist_swapped_r(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_swapped_l(pid_idx))
      else
        ! No rotation, so use the standard distances for this direction
        edge_dist_l = -ABS(panel_edge_dist_l(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_r(pid_idx))
      end if

      ! Set the nearest distance for the left-most column
      if (ABS(edge_dist_l) < ABS(edge_dist_r)) then
        edge_dist_local(rel_idx_k) = edge_dist_l
      else
        edge_dist_local(rel_idx_k) = edge_dist_r
      end if
    end do

    ! In the central cell, there is no rotation --------------------------------
    pid_idx = map_pid(1)
    edge_dist_l = -ABS(panel_edge_dist_l(pid_idx))
    edge_dist_r = ABS(panel_edge_dist_r(pid_idx))

    ! Set the nearest distance for the left-most column
    if (ABS(edge_dist_l) < ABS(edge_dist_r)) then
      edge_dist_local(rel_idx_k) = edge_dist_l
    else
      edge_dist_local(rel_idx_k) = edge_dist_r
    end if

    ! Look at the right-most column in the stencil -----------------------------
    do rel_idx_k = 1, stencil_extent_r
      pid_idx = stencil_map_pid(rel_idx_k)
      if (rotated_panel_r == 1 .and. rel_idx_k > rot_edge_dist_r - 1) then
        ! Clockwise rotation to the right
        edge_dist_l = -ABS(panel_edge_dist_swapped_l(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_swapped_r(pid_idx))
      else if (rotated_panel_r == -1 .and. rel_idx_k > rot_edge_dist_r - 1) then
        ! Anti-clockwise rotation to the right, so swap the directions
        edge_dist_l = -ABS(panel_edge_dist_swapped_r(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_swapped_l(pid_idx))
      else
        ! No rotation, so use the standard distances for this direction
        edge_dist_l = -ABS(panel_edge_dist_l(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_r(pid_idx))
      end if

      ! Set the nearest distance for the right-most column
      if (ABS(edge_dist_l) < ABS(edge_dist_r)) then
        edge_dist_local(rel_idx_k) = edge_dist_l
      else
        edge_dist_local(rel_idx_k) = edge_dist_r
      end if
    end do

    ! If the panel edge is not within the set of cells immediately around
    ! the departure cell, then do not need to do anything special.
    ! Therefore set the distance to be zero
    do rel_idx_k = -stencil_extent_l, stencil_extent_r
      if (ABS(edge_dist_local(rel_idx_k)) > order) then
        edge_dist_local(rel_idx_k) = 0
      end if
    end do

    ! Loop over the direction dofs to compute flux at each dof -----------------
    do dof_iterator = 1, face_selector(map_w3_2d(1))

      ! Pull out index to avoid multiple indirections
      w2h_idx = map_w2h(local_dofs(dof_iterator))
      idx_dep = map_depk(local_dofs(dof_iterator)) + ndep_half

      ! Set a local offset, dependent on the face we are looping over
      select case (local_dofs(dof_iterator))
      case ( W, S )
        dof_offset = 0
      case ( E, N )
        dof_offset = 1
      end select

      ! ====================================================================== !
      ! Extract departure info
      ! ====================================================================== !
      ! Pull out departure point, and separate into integer / frac parts
      ! NB: minus sign in 'direction' because the Y1D stencil runs from S to N
      displacement(:) = direction*dep_dist(w2h_idx:w2h_idx + nlayers-1)
      frac_dist(:) = ABS(displacement(:) - REAL(INT(displacement(:), i_def), r_tran))
      sign_disp(:) = INT(SIGN(1.0_r_tran, displacement(:)))

      ! The relative index of the departure cell
      rel_dep_cell_idx(:) = (                                                  &
        dof_offset - INT(displacement(:), i_def)                               &
        + (1 - sign_disp(:)) / 2 - 1                                           &
      )

      ! ====================================================================== !
      ! INTEGER FLUX
      ! ====================================================================== !
      ! Set flux to be zero initially
      flux(w2h_idx : w2h_idx + nlayers-1) = 0.0_r_tran

      ! Loop over columns in stencil
      do rel_idx_k = rel_idx_min, rel_idx_max
        ! Extract indices for looping over levels. We loop from the lowest
        ! relevant level to the highest
        k_low = dep_lowest_k(idx_dep + rel_idx_k)
        k_high = dep_highest_k(idx_dep + rel_idx_k)

        ! Determine if we have crossed into a rotated panel
        if (rel_idx_k < rot_edge_dist_l + 1 .or. rel_idx_k > rot_edge_dist_r - 1) then
          field_ptr => field_swapped
          mass_ptr => dry_mass_swapped
        else
          field_ptr => field
          mass_ptr => dry_mass
        end if

        ! If k_high is -1, there are no levels to loop over for this column
        if (k_high > -1_i_def) then
          col_idx = stencil_map(rel_idx_k)

          ! Adjust relative index based on sign for each level
          int_cell_idx(k_low+1 : k_high+1) = (                                 &
              1 + sign_disp(k_low+1 : k_high+1)                                &
              * (dof_offset - rel_idx_k - 1                                    &
                + (1 - sign_disp(k_low+1 : k_high+1)) / 2)                     &
          )
          ! Factor for whether to include mass for each cell in calculation
          ! This is one for cells between the departure point and flux point,
          ! but zero for cells outside
          int_cell_switch(k_low+1 : k_high+1) = (                              &
              1 + SIGN(1, ABS(INT(displacement(k_low+1 : k_high+1), i_def))    &
                - int_cell_idx(k_low+1 : k_high+1))                            &
              * SIGN(1, int_cell_idx(k_low+1 : k_high+1) - 1)                  &
          ) / 2

          ! Add integer mass to flux
          flux(w2h_idx+k_low : w2h_idx+k_high) =                               &
              flux(w2h_idx+k_low : w2h_idx+k_high)                             &
              + int_cell_switch(k_low+1 : k_high+1)                            &
              * field_ptr(col_idx+k_low : col_idx+k_high)                      &
              * mass_ptr(col_idx+k_low : col_idx+k_high)
        end if
      end do

      ! ====================================================================== !
      ! Determine panel edge case for each layer
      ! ====================================================================== !
      monotone_same = .true.

      ! For each layer, determine distance from departure cell to panel edge
      rel_dep_cell_idx_bounded(:) = (                                          &
        MIN(stencil_extent_r, MAX(-stencil_extent_l, rel_dep_cell_idx(:)))     &
      )
      do k = 1, nlayers
        edge_dist(k) = edge_dist_local(rel_dep_cell_idx_bounded(k))
      end do

      ! Convert all edge distances into a number describing how to shift the
      ! stencil of cells used in the reconstruction.
      ! This is a number in the range 1:order, and is adjusted based on the
      ! sign of the wind to match the arrangement of data in field_local
      special_edge_case(:) = (                                                 &
          (1 + sign_disp(:))/2 * (edge_dist(:) + order + 1)                    &
          + (1 - sign_disp(:))/2 * (order + 1 - edge_dist(:))                  &
      )

      ! Determine which monotonicity to use
      ! If special_edge_case == order + 1 then no special treatment
      do k = 1, nlayers
        if (special_edge_case(k) /= order + 1) then
          mono_option(k) = edge_monotone
        else
          mono_option(k) = monotone
        end if

        if (k > 1) then
          if (mono_option(k) /= mono_option(k-1)) monotone_same = .false.
        end if
      end do

      ! ====================================================================== !
      ! Populate local arrays for fractional flux calculations
      ! ====================================================================== !
      do j = 1, recon_size
        ! If this column has idx 0, find relative index for the column of the
        ! departure cell, between -stencil_extent_l and stencil_extent_r
        rel_idx(:) = MIN(stencil_extent_r, MAX(-stencil_extent_l,              &
            rel_dep_cell_idx(:) + j - order - order_offset - 1                 &
        ))

        ! Determine whether this column has crossed a panel edge
        ! This factor is 1 if the column has crossed a panel edge, and 0 if it
        ! has not. It is made from (1 - switch_l*switch_r) where:
        ! switch_l is 0 when on a rotated left panel, 1 otherwise
        ! switch_r is 0 when on a rotated right panel, 1 otherwise
        swap_switch(:) = REAL(1                                                &
            - (1 + SIGN(1, rel_idx - rot_edge_dist_l - 1))                     &
            ! 0 if crossed to rotated panel, 1 if not
            * (1 + SIGN(1, rot_edge_dist_r - rel_idx - 1)) / 4, r_tran         &
        )

        ! Swap the order of the local array depending on sign of wind
        ! j if wind > 0, recon_size - j + 1 if wind < 0
        loc_idx = (                                                            &
            (1 + sign_disp)/2*j + (1 - sign_disp)/2*(recon_size - j + 1)       &
        )

        ! Extract reconstruction data
        do k = 1, nlayers
          col_idx = stencil_map(rel_idx(k))
          field_local(k, loc_idx(k)) = (                                       &
            (1.0_r_tran - swap_switch(k))*field(col_idx+k-1)                   &
            + swap_switch(k)*field_swapped(col_idx+k-1)                        &
          )
        end do
      end do

      ! ====================================================================== !
      ! EDGE RECONSTRUCTION
      ! ====================================================================== !
      select case ( order )
      case ( 1 )
        ! Nirvana reconstruction
        if ( high_order_edges ) then
          ! Compute edge values using shifted stencil Nirvana at panel edges
          call nirvana_special_edge(                                           &
                  field_local, special_edge_case,                              &
                  field_edge_left, field_edge_right, nlayers                   &
          )
        else
          ! Compute edge values reverting to linear interpolation at panel edges
          call linear_special_edge(                                            &
                  field_local, special_edge_case,                              &
                  field_edge_left, field_edge_right, nlayers                   &
          )
        end if

      case ( 2 )
        ! PPM reconstruction
        if ( high_order_edges ) then
          ! Compute edge values using shifted stencil PPM at panel edges
          call fourth_order_special_edge(                                      &
                  field_local, special_edge_case,                              &
                  field_edge_left, field_edge_right, nlayers                   &
          )
        else
          ! Compute edge values reverting to Nirvana at panel edges
          call fourth_nirvana_special_edge(                                    &
                  field_local, special_edge_case,                              &
                  field_edge_left, field_edge_right, nlayers                   &
          )
        end if
      end select

      ! ====================================================================== !
      ! FRACTIONAL FLUX RECONSTRUCTION
      ! ====================================================================== !
      if (order == 0) then
        ! Constant reconstruction
        recon_field(:) = field_local(:,1+order_offset)

      else if (monotone_same) then
        ! Can apply monotonicity and reconstruction routines to whole array
        ! at once, as the monotonicity is the same on all levels
        ! ----------------------------------------------------------------------
        ! Apply monotonicity to edges if required
        call monotonic_edge(                                                   &
                field_local, mono_option(1), min_val,                          &
                field_edge_left, field_edge_right,                             &
                order+order_offset, 1, nlayers                                 &
        )
        ! Compute reconstruction using field edge values
        ! and quadratic subgrid reconstruction
        call subgrid_quadratic_recon(                                          &
                recon_field, frac_dist, field_local,                           &
                field_edge_left, field_edge_right, mono_option(1),             &
                order+order_offset, nlayers                                    &
        )
      else
        ! Need to perform monotonicity and reconstruction on each level, as
        ! the monotonicity options are different on different levels
        ! ----------------------------------------------------------------------
        do k = 1, nlayers
          ! Apply monotonicity to edges if required
          field_local_tmp(1,:) = field_local(k,:)
          call monotonic_edge(                                                 &
                  field_local_tmp, mono_option(k), min_val,                    &
                  field_edge_left(k), field_edge_right(k),                     &
                  order+order_offset, 1, 1                                     &
          )
          ! Compute reconstruction using field edge values
          ! and quadratic subgrid reconstruction
          call subgrid_quadratic_recon(                                        &
                  recon_field(k), frac_dist(k), field_local_tmp,               &
                  field_edge_left(k), field_edge_right(k), mono_option(k),     &
                  order+order_offset, 1                                        &
          )
        end do
      end if

      ! ====================================================================== !
      ! Assign flux
      ! ====================================================================== !
      ! NB: minus sign from 'direction' before integer component of flux
      ! returns us to the usual y-direction for the rest of the model
      ! This is not needed for the fractional part, as frac_dry_flux has
      ! the correct sign already
      flux(w2h_idx : w2h_idx+nlayers-1) = inv_dt * (                           &
        recon_field(:) * frac_dry_flux(w2h_idx : w2h_idx+nlayers-1)            &
        + direction * sign_disp(:) * flux(w2h_idx : w2h_idx+nlayers-1)         &
      )
    end do ! dof_iterator

  end subroutine ffsl_flux_xy_special_edge_1d

end module ffsl_flux_xy_special_edge_kernel_mod
