!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the horizontal mass flux for FFSL, with special treatment
!!       at edges of cubed-sphere panels..
!> @details This kernel computes the flux in the x and y directions. A choice of
!!          constant, Nirvana, or PPM is used to compute the edge reconstruction.
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

module ffsl_flux_xy_spt_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    CELL_COLUMN, GH_WRITE,     &
                                    GH_READ, GH_SCALAR,        &
                                    STENCIL, X1D, GH_INTEGER,  &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    Y1D
  use constants_mod,         only : i_def, r_tran, r_def
  use fs_continuity_mod,     only : W3, W2h
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : W, E, N, S

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: ffsl_flux_xy_spt_kernel_type
    private
    type(arg_type) :: meta_args(15) = (/                                       &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),          & ! field_for_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),          & ! dry_mass_for_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(Y1D)),          & ! field_for_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(Y1D)),          & ! dry_mass_for_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! frac_dry_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                                STENCIL(X1D)), & ! panel_id_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                                STENCIL(Y1D)), & ! panel_id_y
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ew
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ns
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                        & ! order
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                        & ! monotone
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                        & ! extent_size
         arg_type(GH_SCALAR, GH_REAL,    GH_READ     )                         & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_flux_xy_spt_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: ffsl_flux_xy_spt_code

contains

  !> @brief Compute the horizontal fluxes for FFSL.
  !> @param[in]     nlayers             Number of layers
  !> @param[in,out] flux                The output flux
  !> @param[in]     field_for_x         Field to use in evaluating x-flux
  !> @param[in]     stencil_size_x      Local length of field_for_x W3 stencil
  !> @param[in]     stencil_map_x       Dofmap for the field_for_x stencil
  !> @param[in]     dry_mass_for_x      Volume or dry mass field at W3 points
  !!                                    for use in evaluating x-flux
  !> @param[in]     stencil_size_mass_x Local length of mass x-stencil
  !> @param[in]     stencil_map_mass_x  Dofmap for the mass x-stencil
  !> @param[in]     field_for_y         Field to use in evaluating x-flux
  !> @param[in]     stencil_size_y      Local length of field_for_y W3 stencil
  !> @param[in]     stencil_map_y       Dofmap for the field_for_y stencil
  !> @param[in]     dry_mass_for_y      Volume or dry mass field at W3 points
  !!                                    for use in evaluating x-flux
  !> @param[in]     stencil_size_mass_y Local length of mass y-stencil
  !> @param[in]     stencil_map_mass_y  Dofmap for the mass y-stencil
  !> @param[in]     dep_dist            Horizontal departure distances
  !> @param[in]     frac_dry_flux       Fractional part of the dry flux or wind
  !> @param[in]     panel_id_x          Panel ID used for the x-flux
  !> @param[in]     stencil_size_px     Local length of panel_id_x stencil
  !> @param[in]     stencil_map_px      Dofmap for the panel_id_x stencil
  !> @param[in]     panel_id_y          Panel ID used for the y-flux
  !> @param[in]     stencil_size_py     Local length of panel_id_y stencil
  !> @param[in]     stencil_map_py      Dofmap for the panel_id_y stencil
  !> @param[in]     face_selector_ew    2D field indicating which W/E faces to
  !!                                    loop over for this column
  !> @param[in]     face_selector_ns    2D field indicating which N/S faces to
  !!                                    loop over for this column
  !> @param[in]     order               Order of reconstruction
  !> @param[in]     monotone            Horizontal monotone option for FFSL
  !> @param[in]     extent_size         Stencil extent needed for the LAM edge
  !> @param[in]     dt                  Time step
  !> @param[in]     ndf_w2h             Num of DoFs for W2h per cell
  !> @param[in]     undf_w2h            Num of DoFs for W2h in this partition
  !> @param[in]     map_w2h             Map for W2h
  !> @param[in]     ndf_w3              Num of DoFs for W3 per cell
  !> @param[in]     undf_w3             Num of DoFs for W3 in this partition
  !> @param[in]     map_w3              Map for W3
  !> @param[in]     ndf_w3_2d           Num of DoFs for 2D W3 per cell
  !> @param[in]     undf_w3_2d          Num of DoFs for this partition for 2D W3
  !> @param[in]     map_w3_2d           Map for 2D W3
  subroutine ffsl_flux_xy_spt_code( nlayers,             &
                                    flux,                &
                                    field_for_x,         &
                                    stencil_size_x,      &
                                    stencil_map_x,       &
                                    dry_mass_for_x,      &
                                    stencil_size_mass_x, &
                                    stencil_map_mass_x,  &
                                    field_for_y,         &
                                    stencil_size_y,      &
                                    stencil_map_y,       &
                                    dry_mass_for_y,      &
                                    stencil_size_mass_y, &
                                    stencil_map_mass_y,  &
                                    dep_dist,            &
                                    frac_dry_flux,       &
                                    panel_id_x,          &
                                    stencil_size_px,     &
                                    stencil_map_px,      &
                                    panel_id_y,          &
                                    stencil_size_py,     &
                                    stencil_map_py,      &
                                    face_selector_ew,    &
                                    face_selector_ns,    &
                                    order,               &
                                    monotone,            &
                                    extent_size,         &
                                    dt,                  &
                                    ndf_w2h,             &
                                    undf_w2h,            &
                                    map_w2h,             &
                                    ndf_w3,              &
                                    undf_w3,             &
                                    map_w3,              &
                                    ndf_pid,             &
                                    undf_pid,            &
                                    map_pid,             &
                                    ndf_w3_2d,           &
                                    undf_w3_2d,          &
                                    map_w3_2d )

    use subgrid_horizontal_support_mod, only: horizontal_nirvana_recon_spt_edges, &
                                              horizontal_ppm_recon_spt_edges

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
    integer(kind=i_def), intent(in) :: stencil_size_x
    integer(kind=i_def), intent(in) :: stencil_size_mass_x
    integer(kind=i_def), intent(in) :: stencil_size_y
    integer(kind=i_def), intent(in) :: stencil_size_mass_y
    integer(kind=i_def), intent(in) :: stencil_size_px
    integer(kind=i_def), intent(in) :: stencil_size_py

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
    integer(kind=i_def), intent(in) :: stencil_map_x(ndf_w3,stencil_size_x)
    integer(kind=i_def), intent(in) :: stencil_map_mass_x(ndf_w3,stencil_size_mass_x)
    integer(kind=i_def), intent(in) :: stencil_map_y(ndf_w3,stencil_size_y)
    integer(kind=i_def), intent(in) :: stencil_map_mass_y(ndf_w3,stencil_size_mass_y)
    integer(kind=i_def), intent(in) :: stencil_map_px(ndf_pid,stencil_size_px)
    integer(kind=i_def), intent(in) :: stencil_map_py(ndf_pid,stencil_size_py)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: field_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: field_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)
    real(kind=r_tran),   intent(in)    :: frac_dry_flux(undf_w2h)
    integer(kind=i_def), intent(in)    :: order
    integer(kind=i_def), intent(in)    :: monotone
    integer(kind=i_def), intent(in)    :: extent_size
    real(kind=r_tran),   intent(in)    :: dt
    real(kind=r_def),    intent(in)    :: panel_id_x(undf_pid)
    real(kind=r_def),    intent(in)    :: panel_id_y(undf_pid)
    integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)

    ! Integers
    integer(kind=i_def) :: local_dofs_x(2), local_dofs_y(2)
    integer(kind=i_def) :: dof_iterator
    integer(kind=i_def) :: stencil_idx, col_idx, pid_idx
    integer(kind=i_def) :: sign_displacement, int_displacement
    integer(kind=i_def) :: rel_idx, dep_cell_idx, rel_dep_cell_idx
    integer(kind=i_def) :: dof_offset, sign_offset
    integer(kind=i_def) :: k, j, idx_3d
    integer(kind=i_def) :: stencil_size, stencil_half, lam_edge_size, recon_size

    ! Reals
    real(kind=r_tran)   :: displacement, frac_dist
    real(kind=r_tran)   :: recon_field, int_mass

    real(kind=r_tran),   allocatable :: field_local(:)
    integer(kind=i_def), allocatable :: ipanel_local(:)

    ! x-direction
    local_dofs_x = (/ W, E /)
    local_dofs_y = (/ S, N /)

    ! set size of array for reconstruction
    recon_size = 3 + 2*order
    allocate(field_local(recon_size))
    allocate(ipanel_local(recon_size))

    ! Get half stencil size
    stencil_size = stencil_size_x
    stencil_half = (stencil_size + 1_i_def) / 2_i_def

    ! Get size the stencil should be to check if we are at the edge of a LAM domain
    lam_edge_size = 2_i_def*extent_size+1_i_def

    if ( lam_edge_size > stencil_size) then

      ! At edge of LAM, so set output to zero ----------------------------------
      do k = 0, nlayers - 1
        do dof_iterator = 1, face_selector_ew(map_w3_2d(1))
          flux( map_w2h(local_dofs_x(dof_iterator)) + k ) = 0.0_r_tran
        end do
        do dof_iterator = 1, face_selector_ns(map_w3_2d(1))
          flux( map_w2h(local_dofs_y(dof_iterator)) + k ) = 0.0_r_tran
        end do
      end do

    else

      ! Not at edge of LAM so compute fluxes -----------------------------------

      ! Loop over the x direction dofs to compute flux at each dof
      do dof_iterator = 1, face_selector_ew(map_w3_2d(1))

        ! Pull out index to avoid multiple indirections
        idx_3d = map_w2h(local_dofs_x(dof_iterator))

        ! Set a local offset, dependent on the face we are looping over
        select case (local_dofs_x(dof_iterator))
        case ( W, S )
          dof_offset = 0
        case ( E, N )
          dof_offset = 1
        end select

        ! Loop over vertical levels
        do k = 0, nlayers - 1

          ! Flux calculation for a single face =================================

          ! Pull out departure point, and separate into integer / frac parts
          displacement = dep_dist(idx_3d + k)
          int_displacement = INT(displacement, i_def)
          frac_dist = displacement - REAL(int_displacement, r_tran)
          sign_displacement = INT(SIGN(1.0_r_tran, displacement))

          ! Set an offset for the stencil index, based on dep point sign
          sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

          ! The local index of the departure cell
          rel_dep_cell_idx = dof_offset - int_displacement + sign_offset - 1

          ! Integer part =======================================================
          int_mass = 0.0_r_tran
          do j = 1, ABS(int_displacement)
            ! If this column has idx 0, find relative index along column of
            ! the departure cell, between - stencil_half and stencil_half
            rel_idx = dof_offset - (j - 1) * sign_displacement + sign_offset - 1

            ! Determine the index in the stencil from rel_idx
            ! e.g. for extent 4:
            ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
            ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
            stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2

            col_idx = stencil_map_x(1,stencil_idx)

            int_mass = int_mass + field_for_x(col_idx + k) * dry_mass_for_x(col_idx + k)

          end do

          ! Fractional part ====================================================
          ! Extract reconstruction data
          do j = 1, recon_size
            ! If this column has idx 0, find relative index along column of
            ! the departure cell, between - stencil_half and stencil_half
            rel_idx = rel_dep_cell_idx + j - order - 2

            ! Determine the index in the stencil from rel_idx
            ! e.g. for extent 4:
            ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
            ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
            stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2

            col_idx = stencil_map_x(1,stencil_idx)
            pid_idx = stencil_map_px(1,stencil_idx)

            ! Populate small local array of field values
            field_local(j) = field_for_x(col_idx + k)
            ipanel_local(j) = int(panel_id_x(pid_idx), i_def)

            ! Store index of departure cell for later
            if (j == order + 2) dep_cell_idx = col_idx

          end do

          select case ( order )
          case ( 0 )
            ! Constant reconstruction
            recon_field = field_local(2)
          case ( 1 )
            ! Nirvana reconstruction
            call horizontal_nirvana_recon_spt_edges(recon_field, frac_dist,    &
                                                    field_local, ipanel_local, &
                                                    monotone)
          case ( 2 )
            ! PPM
            call horizontal_ppm_recon_spt_edges(recon_field, frac_dist,        &
                                                field_local, ipanel_local,     &
                                                monotone)
          end select

          ! Assign flux ========================================================
          flux(idx_3d + k) = (recon_field * frac_dry_flux(idx_3d + k)          &
                              + int_mass * SIGN(1.0_r_tran, displacement)) / dt
        end do ! vertical levels k

      end do ! dof_iterator

      ! ----------------------------------------------------------------------- !

      ! Loop over the y direction dofs to compute flux at each dof
      do dof_iterator = 1, face_selector_ns(map_w3_2d(1))

        ! Pull out index to avoid multiple indirections
        idx_3d = map_w2h(local_dofs_y(dof_iterator))

        ! Set a local offset, dependent on the face we are looping over
        select case (local_dofs_y(dof_iterator))
        case ( W, S )
          dof_offset = 0
        case ( E, N )
          dof_offset = 1
        end select

        ! Loop over vertical levels
        do k = 0, nlayers - 1

          ! Flux calculation for a single face =================================

          ! Pull out departure point, and separate into integer / frac parts
          ! NB: minus signs are needed because the Y1D stencil runs from S to N
          ! but winds blowing S to N have negative sign
          displacement = -dep_dist(idx_3d + k)     ! NB: minus sign
          int_displacement = INT(displacement, i_def)
          frac_dist = displacement - REAL(int_displacement, r_tran)
          sign_displacement = INT(SIGN(1.0_r_tran, displacement))

          ! Set an offset for the stencil index, based on dep point sign
          sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

          ! The local index of the departure cell
          rel_dep_cell_idx = dof_offset - int_displacement + sign_offset - 1

          ! Integer part =======================================================
          int_mass = 0.0_r_tran
          do j = 1, ABS(int_displacement)
            ! If this column has idx 0, find relative index along column of
            ! the departure cell, between - stencil_half and stencil_half
            rel_idx = dof_offset - (j - 1) * sign_displacement + sign_offset - 1

            ! Determine the index in the stencil from rel_idx
            ! e.g. for extent 4:
            ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
            ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
            stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2

            col_idx = stencil_map_y(1,stencil_idx)

            int_mass = int_mass + field_for_y(col_idx + k) * dry_mass_for_y(col_idx + k)

          end do

          ! Fractional part ====================================================
          ! Extract reconstruction data
          do j = 1, recon_size
            ! If this column has idx 0, find relative index along column of
            ! the departure cell, between - stencil_half and stencil_half
            rel_idx = rel_dep_cell_idx + j - order - 2

            ! Determine the index in the stencil from rel_idx
            ! e.g. for extent 4:
            ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
            ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
            stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2

            col_idx = stencil_map_y(1,stencil_idx)
            pid_idx = stencil_map_py(1,stencil_idx)

            ! Populate small local array of field values
            field_local(j) = field_for_y(col_idx + k)
            ipanel_local(j) = int(panel_id_y(pid_idx), i_def)

            ! Store index of departure cell for later
            if (j == order + 2) dep_cell_idx = col_idx

          end do

          select case ( order )
          case ( 0 )
            ! Constant reconstruction
            recon_field = field_local(2)
          case ( 1 )
            ! Nirvana reconstruction
            call horizontal_nirvana_recon_spt_edges(recon_field, frac_dist,    &
                                                    field_local, ipanel_local, &
                                                    monotone)
          case ( 2 )
            ! PPM
            call horizontal_ppm_recon_spt_edges(recon_field, frac_dist,        &
                                                field_local, ipanel_local,     &
                                                monotone)
          end select

          ! Assign flux ========================================================
          ! NB: minus sign before integer component of flux returns us to the
          ! usual y-direction for the rest of the model
          ! This is not needed for the fractional part, as frac_dry_flux has
          ! the correct sign already
          flux(idx_3d + k) = (recon_field * frac_dry_flux(idx_3d + k)        &
                              - int_mass * SIGN(1.0_r_tran, displacement)) / dt
        end do ! vertical levels k

      end do ! dof_iterator

    end if

    deallocate(field_local)
    deallocate(ipanel_local)

  end subroutine ffsl_flux_xy_spt_code

end module ffsl_flux_xy_spt_kernel_mod
