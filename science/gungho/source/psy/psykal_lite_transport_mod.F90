!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief An implementation of the PSy layer for certain transport routines

!> @details Contains implementations of the PSy layer for routines used in
!!          transport methods, which for various reasons give optimisations or
!!          simplifications that are currently not available through PSyclone.
!!          This file is structured as follows:
!!          - Routines relating to the extended mesh. This changes the halo
!!            values of fields so that they correspond to extended mesh panels
!!            on the cubed sphere, to give improved interpolation accuracy by
!!            avoiding the kinks in the coordinate lines on this mesh.
!!          - Fused built-ins that are used in FFSL
!!          - An FFSL panel swap routine, which swaps the halo values of two
!!            fields, which greatly simplifies the horizontal FFSL code.
module psykal_lite_transport_mod

  use field_mod,           only : field_type, field_proxy_type
  use r_tran_field_mod,    only : r_tran_field_type, r_tran_field_proxy_type
  use integer_field_mod,   only : integer_field_type, integer_field_proxy_type
  use constants_mod,       only : r_def, i_def, r_tran, l_def
  use mesh_mod,            only : mesh_type

  implicit none
  public

contains

! ============================================================================ !
! EXTENDED MESH ROUTINES
! ============================================================================ !
! These need psykal_lite implementation because they:
! - loop over halo cells *and* use stencils. This is described by PSyclone issue
!   #2781
! It may be possible to implement some of these routines without psykal_lite
! code by doing redundant computation (ticket #4302)

!>@brief Remap a scalar field from the standard cubed sphere mesh onto an extended
!!       mesh
!!       This routine loops only over halo cells and uses stencils, which is not
!!       currently correctly supported by PSyclone (issue #2781 describes this)
!!       This routine can be removed when PSyclone PR #3089 are picked up
subroutine invoke_init_remap_on_extended_mesh_kernel_type(remap_weights, remap_indices, &
                                                          chi_ext, chi, chi_stencil_depth, &
                                                          panel_id, pid_stencil_depth, &
                                                          linear_remap, ndata)

  use init_remap_on_extended_mesh_kernel_mod, only: init_remap_on_extended_mesh_code
  use function_space_mod,                     only: BASIS, DIFF_BASIS
  use mesh_mod,                               only: mesh_type
  use stencil_2D_dofmap_mod,                  only: stencil_2D_dofmap_type, STENCIL_2D_CROSS

  implicit none

  type(r_tran_field_type), intent(in) :: remap_weights
  type(field_type), intent(in) :: chi_ext(3), chi(3), panel_id
  type(integer_field_type), intent(in) :: remap_indices
  logical(kind=l_def), intent(in) :: linear_remap
  integer(kind=i_def), intent(in) :: ndata
  integer(kind=i_def), intent(in) :: chi_stencil_depth, pid_stencil_depth
  integer(kind=i_def) :: cell
  integer(kind=i_def) :: df_nodal, df_wchi
  real(kind=r_def), allocatable :: basis_wchi(:,:,:)
  integer(kind=i_def) :: dim_wchi
  real(kind=r_def), pointer :: nodes_remap(:,:) => null()
  integer(kind=i_def) :: nlayers
  type(r_tran_field_proxy_type) :: remap_weights_proxy
  type(integer_field_proxy_type) :: remap_indices_proxy
  type(field_proxy_type) :: chi_ext_proxy(3), chi_proxy(3), panel_id_proxy
  integer(kind=i_def), pointer :: map_remap(:,:) => null(), map_panel_id(:,:) => null(), map_wchi(:,:) => null()
  integer(kind=i_def) :: ndf_remap, undf_remap, ndf_wchi, undf_wchi, ndf_panel_id, undf_panel_id
  type(mesh_type), pointer :: mesh => null()
  type(stencil_2d_dofmap_type), pointer :: stencil_map => null()
  integer(kind=i_def), pointer :: wchi_stencil_size(:,:) => null()
  integer(kind=i_def), pointer :: wchi_stencil_dofmap(:,:,:,:) => null()
  integer(kind=i_def)          :: wchi_stencil_max_branch_length
  integer(kind=i_def), pointer :: pid_stencil_size(:,:) => null()
  integer(kind=i_def), pointer :: pid_stencil_dofmap(:,:,:,:) => null()
  integer(kind=i_def)          :: pid_stencil_max_branch_length
  integer(kind=i_def)          :: cell_start, cell_end

  ! Initialise field and/or operator proxies
  remap_weights_proxy = remap_weights%get_proxy()
  remap_indices_proxy = remap_indices%get_proxy()
  chi_ext_proxy(1) = chi_ext(1)%get_proxy()
  chi_ext_proxy(2) = chi_ext(2)%get_proxy()
  chi_ext_proxy(3) = chi_ext(3)%get_proxy()
  chi_proxy(1) = chi(1)%get_proxy()
  chi_proxy(2) = chi(2)%get_proxy()
  chi_proxy(3) = chi(3)%get_proxy()
  panel_id_proxy = panel_id%get_proxy()

  ! Initialise number of layers
  nlayers = remap_weights_proxy%vspace%get_nlayers()

  ! Create a mesh object
  mesh => remap_weights_proxy%vspace%get_mesh()

  ! Initialise stencil dofmaps
  stencil_map => chi_ext_proxy(1)%vspace%get_stencil_2d_dofmap(STENCIL_2D_CROSS, chi_stencil_depth)
  wchi_stencil_max_branch_length = chi_stencil_depth + 1_i_def
  wchi_stencil_dofmap => stencil_map%get_whole_dofmap()
  wchi_stencil_size => stencil_map%get_stencil_sizes()

  stencil_map => panel_id_proxy%vspace%get_stencil_2d_dofmap(STENCIL_2D_CROSS, pid_stencil_depth)
  pid_stencil_max_branch_length = pid_stencil_depth + 1_i_def
  pid_stencil_dofmap => stencil_map%get_whole_dofmap()
  pid_stencil_size => stencil_map%get_stencil_sizes()

  ! Look-up dofmaps for each function space
  map_remap => remap_weights_proxy%vspace%get_whole_dofmap()
  map_wchi => chi_ext_proxy(1)%vspace%get_whole_dofmap()
  map_panel_id => panel_id_proxy%vspace%get_whole_dofmap()

  ! Initialise number of DoFs for remap
  ndf_remap = remap_weights_proxy%vspace%get_ndf()
  undf_remap = remap_weights_proxy%vspace%get_undf()

  ! Initialise number of DoFs for wchi
  ndf_wchi = chi_ext_proxy(1)%vspace%get_ndf()
  undf_wchi = chi_ext_proxy(1)%vspace%get_undf()

  ! Initialise number of DoFs for panel_id
  ndf_panel_id = panel_id_proxy%vspace%get_ndf()
  undf_panel_id = panel_id_proxy%vspace%get_undf()

  ! Initialise evaluator-related quantities for the target function spaces
  nodes_remap => remap_weights_proxy%vspace%get_nodes()

  ! Allocate basis/diff-basis arrays
  dim_wchi = chi_ext_proxy(1)%vspace%get_dim_space()
  allocate (basis_wchi(dim_wchi, ndf_wchi, ndf_remap))

  ! Compute basis/diff-basis arrays
  do df_nodal = 1,ndf_remap
    do df_wchi = 1,ndf_wchi
      basis_wchi(:,df_wchi,df_nodal) = chi_ext_proxy(1)%vspace%call_function(BASIS,df_wchi,nodes_remap(:,df_nodal))
    end do
  end do

  ! Call kernels and communication routines
  if (panel_id_proxy%is_dirty(depth=mesh%get_halo_depth())) THEN
    call panel_id_proxy%halo_exchange(depth=mesh%get_halo_depth())
  end if

  cell_start = mesh%get_last_edge_cell() + 1
  cell_end   = mesh%get_last_halo_cell(mesh%get_halo_depth())

  !$omp parallel default(shared), private(cell)
  !$omp do schedule(static)
  do cell = cell_start, cell_end
    call init_remap_on_extended_mesh_code(nlayers, &
                                      remap_weights_proxy%data, &
                                      remap_indices_proxy%data, &
                                      chi_ext_proxy(1)%data, &
                                      chi_ext_proxy(2)%data, &
                                      chi_ext_proxy(3)%data, &
                                      chi_proxy(1)%data, &
                                      chi_proxy(2)%data, &
                                      chi_proxy(3)%data, &
                                      wchi_stencil_size(:,cell), &
                                      wchi_stencil_dofmap(:,:,:,cell), &
                                      wchi_stencil_max_branch_length, &
                                      panel_id_proxy%data, &
                                      pid_stencil_size(:,cell), &
                                      pid_stencil_dofmap(:,:,:,cell), &
                                      pid_stencil_max_branch_length, &
                                      linear_remap, &
                                      ndata, &
                                      ndf_remap, &
                                      undf_remap, &
                                      map_remap(:,cell), &
                                      ndf_wchi, &
                                      undf_wchi,&
                                      map_wchi(:,cell), &
                                      basis_wchi, &
                                      ndf_panel_id, &
                                      undf_panel_id, map_panel_id(:,cell))
  end do
  !$omp end do

  ! Set halos dirty/clean for fields modified in the above loop
  !$omp master
  call remap_weights_proxy%set_clean(mesh%get_halo_depth())
  call remap_indices_proxy%set_clean(mesh%get_halo_depth())
  !$omp end master
  !
  !$omp end parallel

  ! Deallocate basis arrays
  deallocate (basis_wchi)

end subroutine invoke_init_remap_on_extended_mesh_kernel_type


!>@brief Remap a scalar field from the standard cubed sphere mesh onto an extended
!!       mesh
!!       This routine loops only over halo cells and uses stencils, which is not
!!       currently correctly supported by PSyclone (issue #2781 describes this)
!!       This routine can be removed when PSyclone PR #3089 are picked up
subroutine invoke_remap_on_extended_mesh_kernel_type(remap_field, field, stencil_depth, &
                                                     remap_weights, remap_indices, &
                                                     panel_id, &
                                                     ndata, &
                                                     monotone, enforce_minvalue, minvalue, &
                                                     halo_compute_depth )

  use remap_on_extended_mesh_kernel_mod, only: remap_on_extended_mesh_code
  use mesh_mod,                          only: mesh_type
  use stencil_2D_dofmap_mod,             only: stencil_2D_dofmap_type, STENCIL_2D_CROSS
  implicit none

  type(r_tran_field_type), intent(in) :: remap_field, field, remap_weights
  type(integer_field_type), intent(in) :: remap_indices
  type(field_type), intent(in) :: panel_id
  integer(kind=i_def), intent(in) :: ndata
  logical(kind=l_def), intent(in) :: monotone
  logical(kind=l_def), intent(in) :: enforce_minvalue
  real(kind=r_tran),   intent(in) :: minvalue
  integer(kind=i_def), intent(in) :: halo_compute_depth
  integer(kind=i_def) :: cell, stencil_depth
  integer(kind=i_def) :: nlayers
  type(r_tran_field_proxy_type) :: remap_field_proxy, field_proxy, remap_weights_proxy
  type(integer_field_proxy_type) :: remap_indices_proxy
  type(field_proxy_type) :: panel_id_proxy
  integer(kind=i_def), pointer :: map_remap_field(:,:) => null(), map_panel_id(:,:) => null(), map_remap(:,:) => null()
  integer(kind=i_def) :: ndf_remap_field, undf_remap_field, ndf_remap, undf_remap, ndf_panel_id, undf_panel_id
  type(mesh_type), pointer :: mesh => null()
  type(stencil_2d_dofmap_type), pointer :: stencil_map => null()
  integer(kind=i_def), pointer :: stencil_size(:,:) => null()
  integer(kind=i_def), pointer :: stencil_dofmap(:,:,:,:) => null()
  integer(kind=i_def)          :: stencil_max_branch_length
  integer(kind=i_def)          :: cell_start, cell_end

  ! Initialise field and/or operator proxies
  remap_field_proxy = remap_field%get_proxy()
  field_proxy = field%get_proxy()
  remap_weights_proxy = remap_weights%get_proxy()
  remap_indices_proxy = remap_indices%get_proxy()
  panel_id_proxy = panel_id%get_proxy()

  ! Initialise number of layers
  nlayers = remap_field_proxy%vspace%get_nlayers()

  ! Create a mesh object
  mesh => remap_field_proxy%vspace%get_mesh()

  ! Initialise stencil dofmaps
  stencil_map => field_proxy%vspace%get_stencil_2d_dofmap(STENCIL_2D_CROSS, stencil_depth)
  stencil_max_branch_length = stencil_depth + 1_i_def
  stencil_dofmap => stencil_map%get_whole_dofmap()
  stencil_size => stencil_map%get_stencil_sizes()

  ! Look-up dofmaps for each function space
  map_remap_field => remap_field_proxy%vspace%get_whole_dofmap()
  map_remap => remap_weights_proxy%vspace%get_whole_dofmap()
  map_panel_id => panel_id_proxy%vspace%get_whole_dofmap()

  ! Initialise number of DoFs for remap_field
  ndf_remap_field = remap_field_proxy%vspace%get_ndf()
  undf_remap_field = remap_field_proxy%vspace%get_undf()

  ! Initialise number of DoFs for interpolation fields
  ndf_remap = remap_weights_proxy%vspace%get_ndf()
  undf_remap = remap_weights_proxy%vspace%get_undf()

  ! Initialise number of DoFs for panel_id
  ndf_panel_id = panel_id_proxy%vspace%get_ndf()
  undf_panel_id = panel_id_proxy%vspace%get_undf()

  ! Call kernels and communication routines
  if (field_proxy%is_dirty(depth=mesh%get_halo_depth())) THEN
    call field_proxy%halo_exchange(depth=mesh%get_halo_depth())
  end if
  if (panel_id_proxy%is_dirty(depth=halo_compute_depth)) THEN
    call panel_id_proxy%halo_exchange(depth=halo_compute_depth)
  end if
  cell_start = mesh%get_last_edge_cell() + 1
  cell_end   = mesh%get_last_halo_cell(halo_compute_depth)

  !$omp parallel default(shared), private(cell)
  !$omp do schedule(static)
  do cell = cell_start, cell_end
    call remap_on_extended_mesh_code(nlayers, &
                                      remap_field_proxy%data, &
                                      field_proxy%data, &
                                      stencil_size(:,cell), &
                                      stencil_dofmap(:,:,:,cell), &
                                      stencil_max_branch_length, &
                                      remap_weights_proxy%data, &
                                      remap_indices_proxy%data, &
                                      panel_id_proxy%data, &
                                      ndata, &
                                      monotone, &
                                      enforce_minvalue, &
                                      minvalue, &
                                      ndf_remap_field, &
                                      undf_remap_field, &
                                      map_remap_field(:,cell), &
                                      ndf_remap, &
                                      undf_remap, &
                                      map_remap(:,cell), &
                                      ndf_panel_id, &
                                      undf_panel_id, map_panel_id(:,cell))
  end do
  !$omp end do

  ! Set halos dirty/clean for fields modified in the above loop
  !$omp master
  call remap_field_proxy%set_clean(halo_compute_depth)
  !$omp end master
  !
  !$omp end parallel
end subroutine invoke_remap_on_extended_mesh_kernel_type

!> @brief Computes X_times_Y into the halo cells. Requires a psykal_lite
!!        implementation because:
!!        (a) it needs to run into the halos
!!        (b) the field needs to be marked as clean afterwards
!! Ticket #4302 will investigate replacing this with redundant computation
SUBROUTINE invoke_deep_X_times_Y(field_1, field_2, field_3)
  TYPE(r_tran_field_type), intent(in) :: field_1, field_2, field_3
  INTEGER(KIND=i_def) depth, clean_depth, max_depth
  INTEGER df
  INTEGER(KIND=i_def) loop_start, loop_stop
  TYPE(r_tran_field_proxy_type) field_1_proxy, field_2_proxy, field_3_proxy
  !
  ! Initialise field and/or operator proxies
  !
  field_1_proxy = field_1%get_proxy()
  field_2_proxy = field_2%get_proxy()
  field_3_proxy = field_3%get_proxy()
  !
  ! Determine depth to which original field is clean
  !
  max_depth = MIN(field_1%get_field_halo_depth(), &
                  field_2%get_field_halo_depth(), &
                  field_3%get_field_halo_depth())
  do depth = 1, max_depth
    if (field_2_proxy%is_dirty(depth=depth) .or. &
        field_3_proxy%is_dirty(depth=depth)) then
      ! Dirty
      clean_depth = depth - 1
      exit
    else
      ! Clean
      clean_depth = depth
    end if
  end do
  !
  ! Set-up all of the loop bounds
  !
  loop_start = 1
  loop_stop = field_1_proxy%vspace%get_last_dof_halo(clean_depth)
  !
  DO df = loop_start, loop_stop
    field_1_proxy%data(df) = field_2_proxy%data(df) * field_3_proxy%data(df)
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  CALL field_1_proxy%set_dirty()
  CALL field_1_proxy%set_clean(clean_depth)

END SUBROUTINE invoke_deep_X_times_Y

!> @brief Computes X_divideby_Y into the halo cells. Requires a psykal_lite
!!        implementation because:
!!        (a) it needs to run into the halos
!!        (b) the field needs to be marked as clean afterwards
!! Ticket #4302 will investigate replacing this with redundant computation
SUBROUTINE invoke_deep_X_divideby_Y(field_1, field_2, field_3)
  TYPE(r_tran_field_type), intent(in) :: field_1, field_2, field_3
  INTEGER(KIND=i_def) depth, clean_depth, max_depth
  INTEGER df
  INTEGER(KIND=i_def) loop_start, loop_stop
  TYPE(r_tran_field_proxy_type) field_1_proxy, field_2_proxy, field_3_proxy
  !
  ! Initialise field and/or operator proxies
  !
  field_1_proxy = field_1%get_proxy()
  field_2_proxy = field_2%get_proxy()
  field_3_proxy = field_3%get_proxy()
  !
  ! Determine depth to which original field is clean
  !
  max_depth = MIN(field_1%get_field_halo_depth(), &
                  field_2%get_field_halo_depth(), &
                  field_3%get_field_halo_depth())
  do depth = 1, max_depth
    if (field_2_proxy%is_dirty(depth=depth) .or. &
        field_3_proxy%is_dirty(depth=depth)) then
      ! Dirty
      clean_depth = depth - 1
      exit
    else
      ! Clean
      clean_depth = depth
    end if
  end do
  !
  ! Set-up all of the loop bounds
  !
  loop_start = 1
  loop_stop = field_1_proxy%vspace%get_last_dof_halo(clean_depth)
  !
  ! Call kernels and communication routines
  !
  !
  DO df = loop_start, loop_stop
    field_1_proxy%data(df) = field_2_proxy%data(df) / field_3_proxy%data(df)
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  CALL field_1_proxy%set_dirty()
  CALL field_1_proxy%set_clean(clean_depth)

END SUBROUTINE invoke_deep_X_divideby_Y

!> @brief Computes the shifting of a mass field into the halo cells. Requires
!!        a psykal_lite implementation because:
!!        (a) it needs to run into the halos
!!        (b) the field needs to be marked as clean afterwards
!! Ticket #4302 will investigate replacing this with redundant computation
SUBROUTINE invoke_deep_shift_mass(mass_shifted, mass_prime)
  USE sci_shift_mass_w3_kernel_mod, ONLY: shift_mass_w3_code
  TYPE(r_tran_field_type), intent(in) :: mass_shifted, mass_prime
  INTEGER(KIND=i_def) depth, clean_depth, max_depth
  INTEGER(KIND=i_def) cell
  INTEGER(KIND=i_def) loop0_start, loop0_stop
  INTEGER(KIND=i_def) nlayers
  TYPE(r_tran_field_proxy_type) mass_shifted_proxy, mass_prime_proxy
  INTEGER(KIND=i_def), pointer :: map_adspc3_mass_shifted(:,:) => null(), map_w3(:,:) => null()
  INTEGER(KIND=i_def) ndf_adspc3_mass_shifted, undf_adspc3_mass_shifted, ndf_w3, undf_w3
  TYPE(mesh_type), pointer :: mesh => null()
  !
  ! Initialise field and/or operator proxies
  !
  mass_shifted_proxy = mass_shifted%get_proxy()
  mass_prime_proxy = mass_prime%get_proxy()
  !
  ! Initialise number of layers
  !
  nlayers = mass_shifted_proxy%vspace%get_nlayers()
  !
  ! Create a mesh object
  !
  mesh => mass_shifted_proxy%vspace%get_mesh()
  !
  ! Look-up dofmaps for each function space
  !
  map_adspc3_mass_shifted => mass_shifted_proxy%vspace%get_whole_dofmap()
  map_w3 => mass_prime_proxy%vspace%get_whole_dofmap()
  !
  ! Initialise number of DoFs for adspc3_mass_shifted
  !
  ndf_adspc3_mass_shifted = mass_shifted_proxy%vspace%get_ndf()
  undf_adspc3_mass_shifted = mass_shifted_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for w3
  !
  ndf_w3 = mass_prime_proxy%vspace%get_ndf()
  undf_w3 = mass_prime_proxy%vspace%get_undf()
  !
  ! Determine depth to which original field is clean
  !
  max_depth = MIN(mass_shifted%get_field_halo_depth(), &
                  mass_prime%get_field_halo_depth())
  do depth = 1, max_depth
    if (mass_prime_proxy%is_dirty(depth=depth)) then
      ! Dirty
      clean_depth = depth - 1
      exit
    else
      ! Clean
      clean_depth = depth
    end if
  end do
  !
  ! Set-up all of the loop bounds
  !
  loop0_start = 1
  loop0_stop = mesh%get_last_halo_cell(clean_depth)
  !
  ! Call kernels and communication routines
  !
  DO cell=loop0_start,loop0_stop
    !
    CALL shift_mass_w3_code(nlayers, mass_shifted_proxy%data, mass_prime_proxy%data, &
&ndf_adspc3_mass_shifted, undf_adspc3_mass_shifted, map_adspc3_mass_shifted(:,cell), &
&ndf_w3, undf_w3, map_w3(:,cell))
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  CALL mass_shifted_proxy%set_dirty()
  CALL mass_shifted_proxy%set_clean(clean_depth)
  !
END SUBROUTINE invoke_deep_shift_mass

! ============================================================================ !
! PANEL EDGE REMAPPING ROUTINES
! ============================================================================ !

!> @brief Perform remapping of field near panel edges
!!       This routine loops only over halo cells and uses stencils, which is not
!!       currently correctly supported by PSyclone (issue #2781 describes this)
!!       This routine can be removed when PSyclone PR #3089 are picked up
SUBROUTINE invoke_panel_edge_remap_kernel_type(                                &
              remapped_in_x, remapped_in_y,                                    &
              field_for_x_ptr, cross_depth_x,                                  &
              field_for_y_ptr,  cross_depth_y,                                 &
              panel_edge_weights_x, panel_edge_weights_y,                      &
              panel_edge_indices_x, panel_edge_indices_y, panel_edge_dist,     &
              halo_mask_x, halo_mask_y, compute_all_halo,                      &
              remap_depth, ndata, monotone, enforce_minvalue, minvalue,        &
              ffsl_depth)
  USE panel_edge_remap_kernel_mod, ONLY: panel_edge_remap_code
  USE mesh_mod, ONLY: mesh_type
  USE stencil_2D_dofmap_mod, ONLY: stencil_2D_dofmap_type, STENCIL_2D_CROSS
  REAL(KIND=r_tran), intent(in) :: minvalue
  INTEGER(KIND=i_def), intent(in) :: remap_depth, ndata
  LOGICAL(KIND=l_def), intent(in) :: monotone, enforce_minvalue, compute_all_halo
  TYPE(r_tran_field_type), intent(in) :: remapped_in_x, remapped_in_y, &
                                         field_for_x_ptr, field_for_y_ptr, &
                                         panel_edge_weights_x, panel_edge_weights_y
  TYPE(integer_field_type), intent(in) :: panel_edge_indices_x, panel_edge_indices_y, &
                                          panel_edge_dist(4), halo_mask_x, halo_mask_y
  INTEGER(KIND=i_def), intent(in) :: cross_depth_x, cross_depth_y
  INTEGER, intent(in) :: ffsl_depth
  INTEGER(KIND=i_def) cell
  INTEGER(KIND=i_def) loop0_start, loop0_stop
  INTEGER(KIND=i_def) nlayers_remapped_in_x
  INTEGER(KIND=i_def), pointer, dimension(:) :: panel_edge_dist_1_data => null(), &
                                                panel_edge_dist_2_data => null(), &
                                                panel_edge_dist_3_data => null(), &
                                                panel_edge_dist_4_data => null()
  INTEGER(KIND=i_def), pointer, dimension(:) :: panel_edge_indices_y_data => null()
  INTEGER(KIND=i_def), pointer, dimension(:) :: panel_edge_indices_x_data => null()
  INTEGER(KIND=i_def), pointer, dimension(:) :: halo_mask_x_data => null()
  INTEGER(KIND=i_def), pointer, dimension(:) :: halo_mask_y_data => null()
  TYPE(integer_field_proxy_type) panel_edge_indices_x_proxy, panel_edge_indices_y_proxy, &
                                 panel_edge_dist_proxy(4), halo_mask_x_proxy, halo_mask_y_proxy
  REAL(KIND=r_tran), pointer, dimension(:) :: panel_edge_weights_y_data => null()
  REAL(KIND=r_tran), pointer, dimension(:) :: panel_edge_weights_x_data => null()
  REAL(KIND=r_tran), pointer, dimension(:) :: field_for_y_ptr_data => null()
  REAL(KIND=r_tran), pointer, dimension(:) :: field_for_x_ptr_data => null()
  REAL(KIND=r_tran), pointer, dimension(:) :: remapped_in_y_data => null()
  REAL(KIND=r_tran), pointer, dimension(:) :: remapped_in_x_data => null()
  TYPE(r_tran_field_proxy_type) remapped_in_x_proxy, remapped_in_y_proxy, &
                                field_for_x_ptr_proxy, field_for_y_ptr_proxy, &
                                panel_edge_weights_x_proxy, panel_edge_weights_y_proxy
  INTEGER(KIND=i_def), pointer :: map_adspc1_remapped_in_x(:,:) => null(), &
                                  map_adspc3_panel_edge_dist(:,:) => null(), &
                                  map_adspc5_panel_edge_weights_x(:,:) => null(), &
                                  map_adspc7_field_for_x_ptr(:,:) => null()
  INTEGER(KIND=i_def) ndf_adspc1_remapped_in_x, undf_adspc1_remapped_in_x, &
                      ndf_adspc7_field_for_x_ptr, undf_adspc7_field_for_x_ptr, &
                      ndf_adspc5_panel_edge_weights_x, undf_adspc5_panel_edge_weights_x, &
                      ndf_adspc3_panel_edge_dist, undf_adspc3_panel_edge_dist
  INTEGER(KIND=i_def) max_halo_depth_mesh
  TYPE(mesh_type), pointer :: mesh => null()
  INTEGER(KIND=i_def) field_for_x_ptr_max_branch_length, field_for_y_ptr_max_branch_length
  INTEGER(KIND=i_def), pointer :: field_for_x_ptr_stencil_size(:,:) => null()
  INTEGER(KIND=i_def), pointer :: field_for_y_ptr_stencil_size(:,:) => null()
  INTEGER(KIND=i_def), pointer :: field_for_x_ptr_stencil_dofmap(:,:,:,:) => null()
  INTEGER(KIND=i_def), pointer :: field_for_y_ptr_stencil_dofmap(:,:,:,:) => null()
  TYPE(stencil_2D_dofmap_type), pointer :: field_for_x_ptr_stencil_map => null()
  TYPE(stencil_2D_dofmap_type), pointer :: field_for_y_ptr_stencil_map => null()

  !
  ! Initialise field and/or operator proxies
  !
  remapped_in_x_proxy = remapped_in_x%get_proxy()
  remapped_in_x_data => remapped_in_x_proxy%data
  remapped_in_y_proxy = remapped_in_y%get_proxy()
  remapped_in_y_data => remapped_in_y_proxy%data
  field_for_x_ptr_proxy = field_for_x_ptr%get_proxy()
  field_for_x_ptr_data => field_for_x_ptr_proxy%data
  field_for_y_ptr_proxy = field_for_y_ptr%get_proxy()
  field_for_y_ptr_data => field_for_y_ptr_proxy%data
  panel_edge_weights_x_proxy = panel_edge_weights_x%get_proxy()
  panel_edge_weights_x_data => panel_edge_weights_x_proxy%data
  panel_edge_weights_y_proxy = panel_edge_weights_y%get_proxy()
  panel_edge_weights_y_data => panel_edge_weights_y_proxy%data
  panel_edge_indices_x_proxy = panel_edge_indices_x%get_proxy()
  panel_edge_indices_x_data => panel_edge_indices_x_proxy%data
  panel_edge_indices_y_proxy = panel_edge_indices_y%get_proxy()
  panel_edge_indices_y_data => panel_edge_indices_y_proxy%data
  panel_edge_dist_proxy(1) = panel_edge_dist(1)%get_proxy()
  panel_edge_dist_1_data => panel_edge_dist_proxy(1)%data
  panel_edge_dist_proxy(2) = panel_edge_dist(2)%get_proxy()
  panel_edge_dist_2_data => panel_edge_dist_proxy(2)%data
  panel_edge_dist_proxy(3) = panel_edge_dist(3)%get_proxy()
  panel_edge_dist_3_data => panel_edge_dist_proxy(3)%data
  panel_edge_dist_proxy(4) = panel_edge_dist(4)%get_proxy()
  panel_edge_dist_4_data => panel_edge_dist_proxy(4)%data
  halo_mask_x_proxy = halo_mask_x%get_proxy()
  halo_mask_x_data => halo_mask_x_proxy%data
  halo_mask_y_proxy = halo_mask_y%get_proxy()
  halo_mask_y_data => halo_mask_y_proxy%data
  !
  ! Initialise number of layers
  !
  nlayers_remapped_in_x = remapped_in_x_proxy%vspace%get_nlayers()
  !
  ! Create a mesh object
  !
  mesh => remapped_in_x_proxy%vspace%get_mesh()
  max_halo_depth_mesh = mesh%get_halo_depth()
  !
  ! Initialise stencil dofmaps
  !
  field_for_x_ptr_stencil_map => field_for_x_ptr_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,cross_depth_x)
  field_for_x_ptr_max_branch_length = cross_depth_x + 1_i_def
  field_for_x_ptr_stencil_dofmap => field_for_x_ptr_stencil_map%get_whole_dofmap()
  field_for_x_ptr_stencil_size => field_for_x_ptr_stencil_map%get_stencil_sizes()
  field_for_y_ptr_stencil_map => field_for_y_ptr_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,cross_depth_y)
  field_for_y_ptr_max_branch_length = cross_depth_y + 1_i_def
  field_for_y_ptr_stencil_dofmap => field_for_y_ptr_stencil_map%get_whole_dofmap()
  field_for_y_ptr_stencil_size => field_for_y_ptr_stencil_map%get_stencil_sizes()
  !
  ! Look-up dofmaps for each function space
  !
  map_adspc1_remapped_in_x => remapped_in_x_proxy%vspace%get_whole_dofmap()
  map_adspc7_field_for_x_ptr => field_for_x_ptr_proxy%vspace%get_whole_dofmap()
  map_adspc5_panel_edge_weights_x => panel_edge_weights_x_proxy%vspace%get_whole_dofmap()
  map_adspc3_panel_edge_dist => panel_edge_dist_proxy(1)%vspace%get_whole_dofmap()
  !
  ! Initialise number of DoFs for adspc1_remapped_in_x
  !
  ndf_adspc1_remapped_in_x = remapped_in_x_proxy%vspace%get_ndf()
  undf_adspc1_remapped_in_x = remapped_in_x_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc7_field_for_x_ptr
  !
  ndf_adspc7_field_for_x_ptr = field_for_x_ptr_proxy%vspace%get_ndf()
  undf_adspc7_field_for_x_ptr = field_for_x_ptr_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc5_panel_edge_weights_x
  !
  ndf_adspc5_panel_edge_weights_x = panel_edge_weights_x_proxy%vspace%get_ndf()
  undf_adspc5_panel_edge_weights_x = panel_edge_weights_x_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc3_panel_edge_dist
  !
  ndf_adspc3_panel_edge_dist = panel_edge_dist_proxy(1)%vspace%get_ndf()
  undf_adspc3_panel_edge_dist = panel_edge_dist_proxy(1)%vspace%get_undf()
  !
  ! Set-up all of the loop bounds
  !
  loop0_start = 1
  loop0_stop = mesh%get_last_halo_cell(ffsl_depth)
  !
  ! Call kernels and communication routines
  !
  IF (field_for_x_ptr_proxy%is_dirty(depth=ffsl_depth)) THEN
    CALL field_for_x_ptr_proxy%halo_exchange(depth=ffsl_depth)
  END IF
  IF (field_for_y_ptr_proxy%is_dirty(depth=ffsl_depth)) THEN
    CALL field_for_y_ptr_proxy%halo_exchange(depth=ffsl_depth)
  END IF
  IF (panel_edge_weights_x_proxy%is_dirty(depth=ffsl_depth)) THEN
    CALL panel_edge_weights_x_proxy%halo_exchange(depth=ffsl_depth)
  END IF
  IF (panel_edge_weights_y_proxy%is_dirty(depth=ffsl_depth)) THEN
    CALL panel_edge_weights_y_proxy%halo_exchange(depth=ffsl_depth)
  END IF
  IF (panel_edge_indices_x_proxy%is_dirty(depth=ffsl_depth)) THEN
    CALL panel_edge_indices_x_proxy%halo_exchange(depth=ffsl_depth)
  END IF
  IF (panel_edge_indices_y_proxy%is_dirty(depth=ffsl_depth)) THEN
    CALL panel_edge_indices_y_proxy%halo_exchange(depth=ffsl_depth)
  END IF
  IF (halo_mask_x_proxy%is_dirty(depth=ffsl_depth)) THEN
    CALL halo_mask_x_proxy%halo_exchange(depth=ffsl_depth)
  END IF
  IF (halo_mask_y_proxy%is_dirty(depth=ffsl_depth)) THEN
    CALL halo_mask_y_proxy%halo_exchange(depth=ffsl_depth)
  END IF
  IF (panel_edge_dist_proxy(1)%is_dirty(depth=ffsl_depth)) THEN
    CALL panel_edge_dist_proxy(1)%halo_exchange(depth=ffsl_depth)
  END IF
  IF (panel_edge_dist_proxy(2)%is_dirty(depth=ffsl_depth)) THEN
    CALL panel_edge_dist_proxy(2)%halo_exchange(depth=ffsl_depth)
  END IF
  IF (panel_edge_dist_proxy(3)%is_dirty(depth=ffsl_depth)) THEN
    CALL panel_edge_dist_proxy(3)%halo_exchange(depth=ffsl_depth)
  END IF
  IF (panel_edge_dist_proxy(4)%is_dirty(depth=ffsl_depth)) THEN
    CALL panel_edge_dist_proxy(4)%halo_exchange(depth=ffsl_depth)
  END IF
  DO cell = loop0_start, loop0_stop, 1
    CALL panel_edge_remap_code(nlayers_remapped_in_x, remapped_in_x_data, remapped_in_y_data, field_for_x_ptr_data, &
&field_for_x_ptr_stencil_size(:,cell), field_for_x_ptr_max_branch_length, field_for_x_ptr_stencil_dofmap(:,:,:,cell), &
&field_for_y_ptr_data, field_for_y_ptr_stencil_size(:,cell), field_for_y_ptr_max_branch_length, &
&field_for_y_ptr_stencil_dofmap(:,:,:,cell), panel_edge_weights_x_data, panel_edge_weights_y_data, panel_edge_indices_x_data, &
&panel_edge_indices_y_data, panel_edge_dist_1_data, panel_edge_dist_2_data, panel_edge_dist_3_data, panel_edge_dist_4_data, &
&halo_mask_x_data, halo_mask_y_data, compute_all_halo, &
&remap_depth, ndata, monotone, enforce_minvalue, minvalue, ndf_adspc1_remapped_in_x, undf_adspc1_remapped_in_x, &
&map_adspc1_remapped_in_x(:,cell), ndf_adspc7_field_for_x_ptr, undf_adspc7_field_for_x_ptr, map_adspc7_field_for_x_ptr(:,cell), &
&ndf_adspc5_panel_edge_weights_x, undf_adspc5_panel_edge_weights_x, map_adspc5_panel_edge_weights_x(:,cell), &
&ndf_adspc3_panel_edge_dist, undf_adspc3_panel_edge_dist, map_adspc3_panel_edge_dist(:,cell))
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  CALL remapped_in_x_proxy%set_dirty()
  CALL remapped_in_x_proxy%set_clean(ffsl_depth)
  CALL remapped_in_y_proxy%set_dirty()
  CALL remapped_in_y_proxy%set_clean(ffsl_depth)
  !
  !

END SUBROUTINE invoke_panel_edge_remap_kernel_type

! ============================================================================ !
! ROUTINE TO SET HALO MASK
! ============================================================================ !

  ! This routine requires a psykal lite implementation because it loops over
  ! halo cells and also uses a stencil, which is not currently supported by
  ! PSyclone (issue #2781 describes this)
  ! This routine can be removed when PSyclone PR #3089 are picked up
  SUBROUTINE invoke_halo_mask_xy_kernel_type(halo_mask_x, halo_mask_y, halo_mask, cross_depth, max_halo_depth)
    USE halo_mask_xy_kernel_mod, ONLY: halo_mask_xy_code
    USE mesh_mod, ONLY: mesh_type
    USE stencil_2D_dofmap_mod, ONLY: stencil_2D_dofmap_type, STENCIL_2D_CROSS
    TYPE(integer_field_type), intent(in) :: halo_mask_x, halo_mask_y, halo_mask
    INTEGER(KIND=i_def), intent(in) :: cross_depth
    INTEGER, intent(in) :: max_halo_depth
    INTEGER(KIND=i_def) cell
    INTEGER(KIND=i_def) loop0_start, loop0_stop
    INTEGER(KIND=i_def) nlayers_halo_mask_x
    INTEGER(KIND=i_def), pointer, dimension(:) :: halo_mask_data => null()
    INTEGER(KIND=i_def), pointer, dimension(:) :: halo_mask_y_data => null()
    INTEGER(KIND=i_def), pointer, dimension(:) :: halo_mask_x_data => null()
    TYPE(integer_field_proxy_type) halo_mask_x_proxy, halo_mask_y_proxy, halo_mask_proxy
    INTEGER(KIND=i_def), pointer :: map_adspc3_halo_mask_x(:,:) => null()
    INTEGER(KIND=i_def) ndf_adspc3_halo_mask_x, undf_adspc3_halo_mask_x
    INTEGER(KIND=i_def) max_halo_depth_mesh
    TYPE(mesh_type), pointer :: mesh => null()
    INTEGER(KIND=i_def) halo_mask_max_branch_length
    INTEGER(KIND=i_def), pointer :: halo_mask_stencil_size(:,:) => null()
    INTEGER(KIND=i_def), pointer :: halo_mask_stencil_dofmap(:,:,:,:) => null()
    TYPE(stencil_2D_dofmap_type), pointer :: halo_mask_stencil_map => null()
    !
    ! Initialise field and/or operator proxies
    !
    halo_mask_x_proxy = halo_mask_x%get_proxy()
    halo_mask_x_data => halo_mask_x_proxy%data
    halo_mask_y_proxy = halo_mask_y%get_proxy()
    halo_mask_y_data => halo_mask_y_proxy%data
    halo_mask_proxy = halo_mask%get_proxy()
    halo_mask_data => halo_mask_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_halo_mask_x = halo_mask_x_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => halo_mask_x_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Initialise stencil dofmaps
    !
    halo_mask_stencil_map => halo_mask_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,cross_depth)
    halo_mask_max_branch_length = cross_depth + 1_i_def
    halo_mask_stencil_dofmap => halo_mask_stencil_map%get_whole_dofmap()
    halo_mask_stencil_size => halo_mask_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc3_halo_mask_x => halo_mask_x_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for adspc3_halo_mask_x
    !
    ndf_adspc3_halo_mask_x = halo_mask_x_proxy%vspace%get_ndf()
    undf_adspc3_halo_mask_x = halo_mask_x_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = mesh%get_last_edge_cell()+1
    loop0_stop = mesh%get_last_halo_cell(max_halo_depth)
    !
    ! Call kernels and communication routines
    !
    IF (halo_mask_proxy%is_dirty(depth=max_halo_depth)) THEN
      CALL halo_mask_proxy%halo_exchange(depth=max_halo_depth)
    END IF
    DO cell = loop0_start, loop0_stop, 1
      CALL halo_mask_xy_code(nlayers_halo_mask_x, halo_mask_x_data, halo_mask_y_data, halo_mask_data, &
&halo_mask_stencil_size(:,cell), halo_mask_max_branch_length, halo_mask_stencil_dofmap(:,:,:,cell), ndf_adspc3_halo_mask_x, &
&undf_adspc3_halo_mask_x, map_adspc3_halo_mask_x(:,cell))
    END DO
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    CALL halo_mask_x_proxy%set_dirty()
    CALL halo_mask_x_proxy%set_clean(max_halo_depth)
    CALL halo_mask_y_proxy%set_dirty()
    CALL halo_mask_y_proxy%set_clean(max_halo_depth)
    !
    !
  END SUBROUTINE invoke_halo_mask_xy_kernel_type

end module psykal_lite_transport_mod
