!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Psy-lite code for `adj_convert_hdiv_field_kernel`.
!> Required due to a mismatch in what the kernel is expecting for the BASIS argument,
!> versus what PSyclone provides, requiring manual intervention.
!> Please see PSyclone issue #2798 for further information.
module invoke_adj_cvt_hdiv_field_kernel_mod
  use constants_mod, only: r_def, i_def
  use field_mod,     only: field_type, field_proxy_type

  implicit none

  contains

  subroutine invoke_adj_convert_hdiv_field_kernel(physical_field3, computational_field, chi, panel_id)
      use adj_sci_convert_hdiv_field_kernel_mod, only: adj_convert_hdiv_field_code
      use function_space_mod,                    only: BASIS, DIFF_BASIS
      use mesh_mod,                              only: mesh_type

      type(field_type), intent(in) :: physical_field3(3)
      type(field_type), intent(in) :: computational_field
      type(field_type), intent(in) :: chi(3)
      type(field_type), intent(in) :: panel_id

      integer(kind=i_def)           :: cell
      integer(kind=i_def)           :: colour
      integer(kind=i_def)           :: loop1_start
      integer(kind=i_def)           :: loop0_start, loop0_stop
      integer(kind=i_def)           :: df_aspc2_computational_field, df_aspc9_chi, df_nodal
      real(kind=r_def), allocatable :: basis_aspc2_computational_field_on_aspc1_physical_field3(:,:,:), &
                                       basis_aspc9_chi_on_aspc1_physical_field3(:,:,:), &
                                       diff_basis_aspc9_chi_on_aspc1_physical_field3(:,:,:)

      integer(kind=i_def)           :: dim_aspc2_computational_field, &
                                       dim_aspc9_chi, &
                                       diff_dim_aspc9_chi
      real(kind=r_def), pointer     :: nodes_aspc2_computational_field(:,:) => null()
      real(kind=r_def), pointer     :: nodes_aspc1_physical_field3(:,:) => null()
      integer(kind=i_def)           :: nlayers
      real(kind=r_def), pointer, dimension(:) :: panel_id_data => null()
      real(kind=r_def), pointer, dimension(:) :: chi_1_data => null(), &
                                                 chi_2_data => null(), &
                                                 chi_3_data => null()
      real(kind=r_def), pointer, dimension(:) :: computational_field_data => null()
      real(kind=r_def), pointer, dimension(:) :: physical_field3_1_data => null(), &
                                                 physical_field3_2_data => null(), &
                                                 physical_field3_3_data => null()

      type(field_proxy_type)                  :: physical_field3_proxy(3), &
                                                 computational_field_proxy, &
                                                 chi_proxy(3), &
                                                 panel_id_proxy

      integer(kind=i_def), pointer            :: map_adspc3_panel_id(:,:) => null(), &
                                                 map_aspc1_physical_field3(:,:) => null(), &
                                                 map_aspc2_computational_field(:,:) => null(), &
                                                 map_aspc9_chi(:,:) => null()
      integer(kind=i_def)                     :: ndf_aspc1_physical_field3, &
                                                 undf_aspc1_physical_field3, &
                                                 ndf_aspc2_computational_field, &
                                                 undf_aspc2_computational_field, &
                                                 ndf_aspc9_chi, &
                                                 undf_aspc9_chi, &
                                                 ndf_adspc3_panel_id, &
                                                 undf_adspc3_panel_id
      integer(kind=i_def), allocatable        :: last_halo_cell_all_colours(:,:)
      integer(kind=i_def)                     :: ncolour
      integer(kind=i_def), pointer            :: cmap(:,:)
      integer(kind=i_def)                     :: max_halo_depth_mesh
      type(mesh_type), pointer                :: mesh => null()

      !
      ! Initialise field and/or operator proxies
      !
      physical_field3_proxy(1) = physical_field3(1)%get_proxy()
      physical_field3_1_data => physical_field3_proxy(1)%data
      physical_field3_proxy(2) = physical_field3(2)%get_proxy()
      physical_field3_2_data => physical_field3_proxy(2)%data
      physical_field3_proxy(3) = physical_field3(3)%get_proxy()
      physical_field3_3_data => physical_field3_proxy(3)%data
      computational_field_proxy = computational_field%get_proxy()
      computational_field_data => computational_field_proxy%data
      chi_proxy(1) = chi(1)%get_proxy()
      chi_1_data => chi_proxy(1)%data
      chi_proxy(2) = chi(2)%get_proxy()
      chi_2_data => chi_proxy(2)%data
      chi_proxy(3) = chi(3)%get_proxy()
      chi_3_data => chi_proxy(3)%data
      panel_id_proxy = panel_id%get_proxy()
      panel_id_data => panel_id_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = physical_field3_proxy(1)%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => physical_field3_proxy(1)%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Get the colourmap
      !
      ncolour = mesh%get_ncolours()
      cmap => mesh%get_colour_map()
      !
      ! Look-up dofmaps for each function space
      !
      map_aspc1_physical_field3 => physical_field3_proxy(1)%vspace%get_whole_dofmap()
      map_aspc2_computational_field => computational_field_proxy%vspace%get_whole_dofmap()
      map_aspc9_chi => chi_proxy(1)%vspace%get_whole_dofmap()
      map_adspc3_panel_id => panel_id_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for aspc1_physical_field3
      !
      ndf_aspc1_physical_field3 = physical_field3_proxy(1)%vspace%get_ndf()
      undf_aspc1_physical_field3 = physical_field3_proxy(1)%vspace%get_undf()
      !
      ! Initialise number of DoFs for aspc2_computational_field
      !
      ndf_aspc2_computational_field = computational_field_proxy%vspace%get_ndf()
      undf_aspc2_computational_field = computational_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for aspc9_chi
      !
      ndf_aspc9_chi = chi_proxy(1)%vspace%get_ndf()
      undf_aspc9_chi = chi_proxy(1)%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc3_panel_id
      !
      ndf_adspc3_panel_id = panel_id_proxy%vspace%get_ndf()
      undf_adspc3_panel_id = panel_id_proxy%vspace%get_undf()
      !
      ! Initialise evaluator-related quantities for the target function spaces
      !
      nodes_aspc2_computational_field => computational_field_proxy%vspace%get_nodes()
      nodes_aspc1_physical_field3 => physical_field3_proxy(1)%vspace%get_nodes()
      !
      ! Allocate basis/diff-basis arrays
      !
      dim_aspc2_computational_field = computational_field_proxy%vspace%get_dim_space()
      dim_aspc9_chi = chi_proxy(1)%vspace%get_dim_space()
      diff_dim_aspc9_chi = chi_proxy(1)%vspace%get_dim_space_diff()
      allocate(basis_aspc2_computational_field_on_aspc1_physical_field3(dim_aspc2_computational_field, &
                                                                        ndf_aspc2_computational_field, &
                                                                        ndf_aspc1_physical_field3))
      allocate(basis_aspc9_chi_on_aspc1_physical_field3(dim_aspc9_chi, ndf_aspc9_chi, ndf_aspc1_physical_field3))
      allocate(diff_basis_aspc9_chi_on_aspc1_physical_field3(diff_dim_aspc9_chi, ndf_aspc9_chi, ndf_aspc1_physical_field3))
      !
      ! Compute basis/diff-basis arrays
      !
      do df_nodal=1,ndf_aspc1_physical_field3
        do df_aspc2_computational_field=1,ndf_aspc2_computational_field
          basis_aspc2_computational_field_on_aspc1_physical_field3(:,df_aspc2_computational_field,df_nodal) = &
            computational_field_proxy%vspace%call_function(BASIS, &
                                                           df_aspc2_computational_field, &
                                                           nodes_aspc1_physical_field3(:,df_nodal))
        end do
      end do
      do df_nodal=1,ndf_aspc1_physical_field3
        do df_aspc9_chi=1,ndf_aspc9_chi
          basis_aspc9_chi_on_aspc1_physical_field3(:,df_aspc9_chi,df_nodal) = &
            chi_proxy(1)%vspace%call_function(BASIS, df_aspc9_chi, nodes_aspc1_physical_field3(:,df_nodal))
        end do
      end do
      do df_nodal=1,ndf_aspc1_physical_field3
        do df_aspc9_chi=1,ndf_aspc9_chi
          diff_basis_aspc9_chi_on_aspc1_physical_field3(:,df_aspc9_chi,df_nodal) = &
            chi_proxy(1)%vspace%call_function(DIFF_BASIS, df_aspc9_chi, nodes_aspc1_physical_field3(:,df_nodal))
        end do
      end do
      !
      ! Initialise mesh properties
      !
      last_halo_cell_all_colours = mesh%get_last_halo_cell_all_colours()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = ncolour
      loop1_start = 1
      !
      ! Call kernels and communication routines
      !
      if (physical_field3_proxy(1)%is_dirty(depth=1)) then
        call physical_field3_proxy(1)%halo_exchange(depth=1)
      end if
      if (physical_field3_proxy(2)%is_dirty(depth=1)) then
        call physical_field3_proxy(2)%halo_exchange(depth=1)
      end if
      if (physical_field3_proxy(3)%is_dirty(depth=1)) then
        call physical_field3_proxy(3)%halo_exchange(depth=1)
      end if
      if (chi_proxy(1)%is_dirty(depth=1)) then
        call chi_proxy(1)%halo_exchange(depth=1)
      end if
      if (chi_proxy(2)%is_dirty(depth=1)) then
        call chi_proxy(2)%halo_exchange(depth=1)
      end if
      if (chi_proxy(3)%is_dirty(depth=1)) then
        call chi_proxy(3)%halo_exchange(depth=1)
      end if
      if (panel_id_proxy%is_dirty(depth=1)) then
        call panel_id_proxy%halo_exchange(depth=1)
      end if
      do colour = loop0_start, loop0_stop, 1
        !$omp parallel default(shared), private(cell)
        !$omp do schedule(static)
        do cell = loop1_start, last_halo_cell_all_colours(colour,1), 1
        call adj_convert_hdiv_field_code(nlayers,                                                  &
                                         physical_field3_1_data,                                   &
                                         physical_field3_2_data,                                   &
                                         physical_field3_3_data,                                   &
                                         computational_field_data,                                 &
                                         chi_1_data,                                               &
                                         chi_2_data,                                               &
                                         chi_3_data,                                               &
                                         panel_id_data,                                            &
                                         ndf_aspc1_physical_field3,                                &
                                         undf_aspc1_physical_field3,                               &
                                         map_aspc1_physical_field3(:,cmap(colour,cell)),           &
                                         ndf_aspc2_computational_field,                            &
                                         undf_aspc2_computational_field,                           &
                                         map_aspc2_computational_field(:,cmap(colour,cell)),       &
                                         basis_aspc2_computational_field_on_aspc1_physical_field3, &
                                         ndf_aspc9_chi,                                            &
                                         undf_aspc9_chi,                                           &
                                         map_aspc9_chi(:,cmap(colour,cell)),                       &
                                         basis_aspc9_chi_on_aspc1_physical_field3,                 &
                                         diff_basis_aspc9_chi_on_aspc1_physical_field3,            &
                                         ndf_adspc3_panel_id,                                      &
                                         undf_adspc3_panel_id,                                     &
                                         map_adspc3_panel_id(:,cmap(colour,cell)))
        end do
        !$omp end do
        !$omp end parallel
      end do
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      call computational_field_proxy%set_dirty()
      !
      ! Deallocate basis arrays
      !
      deallocate(basis_aspc2_computational_field_on_aspc1_physical_field3, &
                 basis_aspc9_chi_on_aspc1_physical_field3, &
                 diff_basis_aspc9_chi_on_aspc1_physical_field3)
    end subroutine invoke_adj_convert_hdiv_field_kernel

end module invoke_adj_cvt_hdiv_field_kernel_mod
