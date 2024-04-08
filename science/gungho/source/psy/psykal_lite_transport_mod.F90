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
! - loop over only halo cells
! - need to mark the fields as clean afterwards, so that the values in the halo
!   cells don't get overwritten
! It may be possible to implement some of these routines without psykal_lite
! code by doing redundant computation (ticket #4302)

!>@brief Remap a scalar field from the standard cubed sphere mesh onto an extended
!!       mesh
!!       This routine loops only over halo cells (and always to the full
!!       depth of the halo). It is not clear if this functionality will ever
!!       be supported by psyclone or if it should be treated as a special
!!       case (as with computation of coordinate fields). Issue 2300 has been
!!       opened to investigate this.
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
!!       This routine loops only over halo cells (and always to the full
!!       depth of the halo). It is not clear if this functionality will ever
!!       be supported by psyclone or if it should be treated as a special
!!       case (as with computation of coordinate fields). Issue 2300 has been
!!       opened to investigate this.
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

!------------------------------------------------------------------------------

!> @brief Computes the Det(J) at W3 field on the extended mesh, which requires a
!!        psykal_lite implementation because:
!!        (a) it needs to run up to the full halo depth
!!        (b) the field needs to be marked as clean afterwards
SUBROUTINE invoke_extended_detj_at_w3_kernel_type(detj_at_w3_r_tran, chi_list, panel_id_list)
  USE calc_detj_at_w3_kernel_mod, ONLY: calc_detj_at_w3_code
  USE function_space_mod, ONLY: BASIS, DIFF_BASIS
  TYPE(r_tran_field_type), intent(in) :: detj_at_w3_r_tran
  TYPE(field_type), intent(in) :: chi_list(3), panel_id_list
  INTEGER(KIND=i_def) cell
  INTEGER(KIND=i_def) loop0_start, loop0_stop
  INTEGER(KIND=i_def) df_aspc1_chi_list, df_nodal
  REAL(KIND=r_def), allocatable :: basis_aspc1_chi_list_on_w3(:,:,:), diff_basis_aspc1_chi_list_on_w3(:,:,:)
  INTEGER(KIND=i_def) dim_aspc1_chi_list, diff_dim_aspc1_chi_list
  REAL(KIND=r_def), pointer :: nodes_w3(:,:) => null()
  INTEGER(KIND=i_def) nlayers
  TYPE(field_proxy_type) chi_list_proxy(3), panel_id_list_proxy
  TYPE(r_tran_field_proxy_type) detj_at_w3_r_tran_proxy
  INTEGER(KIND=i_def), pointer :: map_adspc3_panel_id_list(:,:) => null(), map_aspc1_chi_list(:,:) => null(), &
&map_w3(:,:) => null()
  INTEGER(KIND=i_def) ndf_w3, undf_w3, ndf_aspc1_chi_list, undf_aspc1_chi_list, ndf_adspc3_panel_id_list, &
&undf_adspc3_panel_id_list
  INTEGER(KIND=i_def) max_halo_depth_mesh
  TYPE(mesh_type), pointer :: mesh => null()
  !
  ! Initialise field and/or operator proxies
  !
  detj_at_w3_r_tran_proxy = detj_at_w3_r_tran%get_proxy()
  chi_list_proxy(1) = chi_list(1)%get_proxy()
  chi_list_proxy(2) = chi_list(2)%get_proxy()
  chi_list_proxy(3) = chi_list(3)%get_proxy()
  panel_id_list_proxy = panel_id_list%get_proxy()
  !
  ! Initialise number of layers
  !
  nlayers = detj_at_w3_r_tran_proxy%vspace%get_nlayers()
  !
  ! Create a mesh object
  !
  mesh => detj_at_w3_r_tran_proxy%vspace%get_mesh()
  max_halo_depth_mesh = mesh%get_halo_depth()
  !
  ! Look-up dofmaps for each function space
  !
  map_w3 => detj_at_w3_r_tran_proxy%vspace%get_whole_dofmap()
  map_aspc1_chi_list => chi_list_proxy(1)%vspace%get_whole_dofmap()
  map_adspc3_panel_id_list => panel_id_list_proxy%vspace%get_whole_dofmap()
  !
  ! Initialise number of DoFs for w3
  !
  ndf_w3 = detj_at_w3_r_tran_proxy%vspace%get_ndf()
  undf_w3 = detj_at_w3_r_tran_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for aspc1_chi_list
  !
  ndf_aspc1_chi_list = chi_list_proxy(1)%vspace%get_ndf()
  undf_aspc1_chi_list = chi_list_proxy(1)%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc3_panel_id_list
  !
  ndf_adspc3_panel_id_list = panel_id_list_proxy%vspace%get_ndf()
  undf_adspc3_panel_id_list = panel_id_list_proxy%vspace%get_undf()
  !
  ! Initialise evaluator-related quantities for the target function spaces
  !
  nodes_w3 => detj_at_w3_r_tran_proxy%vspace%get_nodes()
  !
  ! Allocate basis/diff-basis arrays
  !
  dim_aspc1_chi_list = chi_list_proxy(1)%vspace%get_dim_space()
  diff_dim_aspc1_chi_list = chi_list_proxy(1)%vspace%get_dim_space_diff()
  ALLOCATE (basis_aspc1_chi_list_on_w3(dim_aspc1_chi_list, ndf_aspc1_chi_list, ndf_w3))
  ALLOCATE (diff_basis_aspc1_chi_list_on_w3(diff_dim_aspc1_chi_list, ndf_aspc1_chi_list, ndf_w3))
  !
  ! Compute basis/diff-basis arrays
  !
  DO df_nodal=1,ndf_w3
    DO df_aspc1_chi_list=1,ndf_aspc1_chi_list
      basis_aspc1_chi_list_on_w3(:,df_aspc1_chi_list,df_nodal) = &
&chi_list_proxy(1)%vspace%call_function(BASIS,df_aspc1_chi_list,nodes_w3(:,df_nodal))
    END DO
  END DO
  DO df_nodal=1,ndf_w3
    DO df_aspc1_chi_list=1,ndf_aspc1_chi_list
      diff_basis_aspc1_chi_list_on_w3(:,df_aspc1_chi_list,df_nodal) = &
&chi_list_proxy(1)%vspace%call_function(DIFF_BASIS,df_aspc1_chi_list,nodes_w3(:,df_nodal))
    END DO
  END DO
  !
  ! Set-up all of the loop bounds
  !
  loop0_start = 1
  loop0_stop = mesh%get_last_halo_cell(mesh%get_halo_depth())
  !
  ! Call kernels and communication routines
  !
  DO cell=loop0_start,loop0_stop
    !
    CALL calc_detj_at_w3_code(nlayers, detj_at_w3_r_tran_proxy%data, chi_list_proxy(1)%data, chi_list_proxy(2)%data, &
&chi_list_proxy(3)%data, panel_id_list_proxy%data, ndf_w3, undf_w3, map_w3(:,cell), ndf_aspc1_chi_list, undf_aspc1_chi_list, &
&map_aspc1_chi_list(:,cell), basis_aspc1_chi_list_on_w3, diff_basis_aspc1_chi_list_on_w3, ndf_adspc3_panel_id_list, &
&undf_adspc3_panel_id_list, map_adspc3_panel_id_list(:,cell))
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  call detj_at_w3_r_tran_proxy%set_clean(mesh%get_halo_depth())
  !
  !
  ! Deallocate basis arrays
  !
  DEALLOCATE (basis_aspc1_chi_list_on_w3, diff_basis_aspc1_chi_list_on_w3)
  !
END SUBROUTINE invoke_extended_detj_at_w3_kernel_type

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
  USE shift_mass_w3_kernel_mod, ONLY: shift_mass_w3_code
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

!============================================================================= !
! BUILT-IN FUSION
!============================================================================= !
! These are specific routines to optimise a chain of built-ins that are used
! in FFSL. In this case, we know that all fields are of the same length and
! so the operation can be contained in a single do-loop, rather than being
! performed in a long sequence of separate do-loops.
! It should be possible to generate these via loop fusion (ticket #4301)

!> @brief Compute the advective increment for a 1D FFSL step. Following
!!        Putman and Lin, JCP, 2007 this calculates:
!!        inc = 1/dt * (field - (field - dt*inc)/adv_one)
!! Ticket #4301 will investigate replacing this loop fusion
SUBROUTINE invoke_ffsl_advective_increment(increment, field, dt, adv_one)
  REAL(KIND=r_tran), intent(in) :: dt
  TYPE(r_tran_field_type), intent(in) :: increment, field, adv_one
  INTEGER df
  INTEGER(KIND=i_def) loop_start, loop_stop
  TYPE(r_tran_field_proxy_type) increment_proxy, field_proxy, adv_one_proxy
  !
  ! Initialise field and/or operator proxies
  !
  increment_proxy = increment%get_proxy()
  field_proxy = field%get_proxy()
  adv_one_proxy = adv_one%get_proxy()
  !
  ! Set-up all of the loop bounds
  !
  loop_start = 1
  loop_stop = increment_proxy%vspace%get_last_dof_annexed()
  !
  ! Call kernels and communication routines
  !
  DO df = loop_start, loop_stop
    increment_proxy%data(df) =                                                 &
      increment_proxy%data(df) / adv_one_proxy%data(df)                        &
      + field_proxy%data(df) / dt * (1.0_r_tran - 1.0_r_tran / adv_one_proxy%data(df) )
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  CALL increment_proxy%set_dirty()
  !
  !
END SUBROUTINE invoke_ffsl_advective_increment

!> @brief Compute the outer update of a tracer using SWIFT. Calculates:
!!        f_np1 = 1/rho_np1 * ( 0.5*(f_x*rho_x + f_y*rho_y)
!!                             - dt*(inc_x + inc_y) )
!!        but using mass this gives
!!        f_np1 = 1/dry_mass_np1 * ( 0.5*(f_x*dry_mass_x + f_y*dry_mass_y)
!!                                  - dt*(mass_inc_x + mass_inc_y) )
!! Ticket #4301 will investigate replacing this loop fusion
SUBROUTINE invoke_swift_outer_update_tracer(field_np1, field_x, field_y,       &
                                            dry_mass_np1, dry_mass_x,          &
                                            dry_mass_y, increment_x,           &
                                            increment_y, dt)
  REAL(KIND=r_tran), intent(in) :: dt
  TYPE(r_tran_field_type), intent(in) :: field_np1, field_x, field_y
  TYPE(r_tran_field_type), intent(in) :: dry_mass_np1, dry_mass_x, dry_mass_y
  TYPE(r_tran_field_type), intent(in) :: increment_x, increment_y
  INTEGER df
  INTEGER(KIND=i_def) loop_start, loop_stop
  TYPE(r_tran_field_proxy_type) field_np1_proxy, field_x_proxy, field_y_proxy
  TYPE(r_tran_field_proxy_type) dry_mass_np1_proxy, dry_mass_x_proxy, dry_mass_y_proxy
  TYPE(r_tran_field_proxy_type) increment_x_proxy, increment_y_proxy
  !
  ! Initialise field and/or operator proxies
  !
  increment_x_proxy = increment_x%get_proxy()
  increment_y_proxy = increment_y%get_proxy()
  field_np1_proxy = field_np1%get_proxy()
  field_x_proxy = field_x%get_proxy()
  field_y_proxy = field_y%get_proxy()
  dry_mass_np1_proxy = dry_mass_np1%get_proxy()
  dry_mass_x_proxy = dry_mass_x%get_proxy()
  dry_mass_y_proxy = dry_mass_y%get_proxy()
  !
  ! Set-up all of the loop bounds
  !
  loop_start = 1
  loop_stop = field_np1_proxy%vspace%get_last_dof_annexed()
  !
  ! Call kernels and communication routines
  !
  DO df = loop_start, loop_stop
    field_np1_proxy%data(df) =                                                 &
      ( field_x_proxy%data(df)*dry_mass_x_proxy%data(df)                       &
        + field_y_proxy%data(df)*dry_mass_y_proxy%data(df)                     &
        - dt * (increment_x_proxy%data(df) + increment_y_proxy%data(df)) )     &
      / dry_mass_np1_proxy%data(df) * 0.5_r_tran
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  CALL field_np1_proxy%set_dirty()
  !
  !
END SUBROUTINE invoke_swift_outer_update_tracer

!> @brief Compute the inner update of a tracer using SWIFT. Calculates:
!!        f_np1 = 1/rho_np1 * ( f_n*rho_n - dt*increment )
!!        but using mass, which gives
!!        f_np1 = 1/dry_mass_np1 * ( f_n*dry_mass_n - dt*mass_increment )
!! Ticket #4301 will investigate replacing this loop fusion
SUBROUTINE invoke_swift_inner_update_tracer(field_np1, field_n, dry_mass_np1,  &
                                            dry_mass_n, increment, dt)
  REAL(KIND=r_tran), intent(in) :: dt
  TYPE(r_tran_field_type), intent(in) :: field_np1, field_n
  TYPE(r_tran_field_type), intent(in) :: dry_mass_np1, dry_mass_n
  TYPE(r_tran_field_type), intent(in) :: increment
  INTEGER df
  INTEGER(KIND=i_def) loop_start, loop_stop
  TYPE(r_tran_field_proxy_type) field_np1_proxy, field_n_proxy
  TYPE(r_tran_field_proxy_type) dry_mass_np1_proxy, dry_mass_n_proxy
  TYPE(r_tran_field_proxy_type) increment_proxy
  !
  ! Initialise field and/or operator proxies
  !
  increment_proxy = increment%get_proxy()
  field_np1_proxy = field_np1%get_proxy()
  field_n_proxy = field_n%get_proxy()
  dry_mass_np1_proxy = dry_mass_np1%get_proxy()
  dry_mass_n_proxy = dry_mass_n%get_proxy()
  !
  ! Set-up all of the loop bounds
  !
  loop_start = 1
  loop_stop = increment_proxy%vspace%get_last_dof_annexed()
  !
  ! Call kernels and communication routines
  !
  DO df = loop_start, loop_stop
    field_np1_proxy%data(df) =                                                 &
      ( field_n_proxy%data(df)*dry_mass_n_proxy%data(df)                       &
      - dt * increment_proxy%data(df) ) / dry_mass_np1_proxy%data(df)
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  CALL field_np1_proxy%set_dirty()
  !
  !
END SUBROUTINE invoke_swift_inner_update_tracer

!============================================================================= !
! FFSL PANEL SWAP
!============================================================================= !
! A special routine to FFSL, which involves swapping the values of two fields in
! their halos. This greatly simplifies the code as outer FFSL fluxes can be
! computed with the same kernels as inner steps.
! Requires a psykal_lite implementation because:
! - it loops only through halo cells
! - fields need to be marked as clean afterwards

!> @brief Swap the halo values of two fields
SUBROUTINE invoke_ffsl_panel_swap_kernel_type(field_x, field_y, panel_id, &
                                              stencil_extent)
  USE ffsl_panel_swap_kernel_mod, ONLY: ffsl_panel_swap_code
  TYPE(r_tran_field_type), intent(in) :: field_x, field_y
  TYPE(field_type), intent(in) :: panel_id
  INTEGER(KIND=i_def), intent(in) :: stencil_extent
  INTEGER(KIND=i_def) cell, depth
  INTEGER(KIND=i_def) loop0_start, loop0_stop
  INTEGER(KIND=i_def) nlayers
  TYPE(field_proxy_type) panel_id_proxy
  TYPE(r_tran_field_proxy_type) field_x_proxy, field_y_proxy
  INTEGER(KIND=i_def), pointer :: map_adspc3_panel_id(:,:) => null(), map_w3(:,:) => null()
  INTEGER(KIND=i_def) ndf_w3, undf_w3, ndf_adspc3_panel_id, undf_adspc3_panel_id
  INTEGER(KIND=i_def) max_halo_depth_mesh
  TYPE(mesh_type), pointer :: mesh => null()
  !
  ! Initialise field and/or operator proxies
  !
  field_x_proxy = field_x%get_proxy()
  field_y_proxy = field_y%get_proxy()
  panel_id_proxy = panel_id%get_proxy()
  !
  ! Initialise number of layers
  !
  nlayers = field_x_proxy%vspace%get_nlayers()
  !
  ! Create a mesh object
  !
  mesh => field_x_proxy%vspace%get_mesh()
  max_halo_depth_mesh = mesh%get_halo_depth()
  depth = MIN(max_halo_depth_mesh, stencil_extent)
  !
  ! Look-up dofmaps for each function space
  !
  map_w3 => field_x_proxy%vspace%get_whole_dofmap()
  map_adspc3_panel_id => panel_id_proxy%vspace%get_whole_dofmap()
  !
  ! Initialise number of DoFs for w3
  !
  ndf_w3 = field_x_proxy%vspace%get_ndf()
  undf_w3 = field_x_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc3_panel_id
  !
  ndf_adspc3_panel_id = panel_id_proxy%vspace%get_ndf()
  undf_adspc3_panel_id = panel_id_proxy%vspace%get_undf()
  !
  ! Set-up all of the loop bounds
  !
  loop0_start = mesh%get_last_edge_cell() + 1
  loop0_stop = mesh%get_last_halo_cell(depth)
  !
  ! Call kernels and communication routines
  !
  IF (field_x_proxy%is_dirty(depth=depth)) THEN
    CALL field_x_proxy%halo_exchange(depth=depth)
  END IF
  IF (field_y_proxy%is_dirty(depth=depth)) THEN
    CALL field_y_proxy%halo_exchange(depth=depth)
  END IF
  !
  !
  DO cell=loop0_start,loop0_stop
    !
    CALL ffsl_panel_swap_code(nlayers, field_x_proxy%data, field_y_proxy%data, panel_id_proxy%data, ndf_w3, undf_w3, &
&map_w3(:,cell), ndf_adspc3_panel_id, undf_adspc3_panel_id, map_adspc3_panel_id(:,cell))
  END DO
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  CALL field_x_proxy%set_dirty()
  CALL field_y_proxy%set_dirty()
  CALL field_x_proxy%set_clean(depth)
  CALL field_y_proxy%set_clean(depth)
  !
  !
END SUBROUTINE invoke_ffsl_panel_swap_kernel_type

end module psykal_lite_transport_mod
