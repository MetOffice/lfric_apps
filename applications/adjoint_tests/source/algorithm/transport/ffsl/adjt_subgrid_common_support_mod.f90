!------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief Module containing adjoint test for adj_subgrid_common_support_mod routines
module adjt_subgrid_common_support_mod

  use constants_mod,                  only : i_def, r_tran, EPS
  use transport_enumerated_types_mod, only : monotone_none
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_DEBUG

  implicit none

  public :: adjt_subgrid_quadratic_recon

  contains
  !=============================================================================
  !> @brief   Adjoint test for adj_subgrid_quadratic_recon.
  !> @details Passes if adjoint is transpose of tangent linear.
  !>          Determined by testing the equality of inner products <Mx, Mx> and <AMx, x>,
  !>          where M is the tangent linear and A is the adjoint.
  !> @param[in] nlayers     Number of layers

  subroutine adjt_subgrid_quadratic_recon(nlayers)

    use subgrid_common_support_mod,     only : subgrid_quadratic_recon
    use adj_subgrid_common_support_mod, only : adj_subgrid_quadratic_recon

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers

    ! Arguments for tl and adj calls
    real(kind=r_tran), dimension(:),   allocatable :: reconstruction
    real(kind=r_tran), dimension(:),   allocatable :: dep
    real(kind=r_tran), dimension(:,:), allocatable :: field
    real(kind=r_tran), dimension(:),   allocatable :: edge_left
    real(kind=r_tran), dimension(:),   allocatable :: edge_right
    integer(kind=i_def)                            :: monotone
    integer(kind=i_def)                            :: order

    ! Copies of input fields used in inner products
    real(kind=r_tran), dimension(:),   allocatable :: recon_inp
    real(kind=r_tran), dimension(:,:), allocatable :: fld_inp
    real(kind=r_tran), dimension(:),   allocatable :: edge_l_inp
    real(kind=r_tran), dimension(:),   allocatable :: edge_r_inp

    ! Inner products
    integer(kind=i_def) :: i_order
    real(kind=r_tran)   :: recon_inner_prod
    real(kind=r_tran)   :: field_inner_prod
    real(kind=r_tran)   :: edge_l_inner_prod
    real(kind=r_tran)   :: edge_r_inner_prod
    real(kind=r_tran)   :: recon_sf
    real(kind=r_tran)   :: field_sf
    real(kind=r_tran)   :: edge_l_sf
    real(kind=r_tran)   :: edge_r_sf
    real(kind=r_tran)   :: inner1
    real(kind=r_tran)   :: recon_recon_inp_inner_prod
    real(kind=r_tran)   :: fld_fld_inp_inner_prod
    real(kind=r_tran)   :: edge_l_edge_l_inp_inner_prod
    real(kind=r_tran)   :: edge_r_edge_r_inp_inner_prod
    real(kind=r_tran)   :: inner2

    ! Test parameters and variables
    real(kind=r_tran), parameter :: overall_tolerance = 1500.0_r_tran
    real(kind=r_tran)            :: machine_tol
    real(kind=r_tran)            :: relative_diff

    order = 0

    allocate(reconstruction(1:nlayers),           &
             dep(1:nlayers),                      &
             field(1:nlayers, 1:(2*order + 1)),   &
             edge_left(1:nlayers),                &
             edge_right(1:nlayers),               &
             recon_inp(1:nlayers),                &
             fld_inp(1:nlayers, 1:(2*order + 1)), &
             edge_l_inp(1:nlayers),               &
             edge_r_inp(1:nlayers))

    ! Initialise arguments and call the forward routine.
    call RANDOM_NUMBER(reconstruction)
    recon_inp = reconstruction
    call RANDOM_NUMBER(dep)
    call RANDOM_NUMBER(field)
    fld_inp = field
    call RANDOM_NUMBER(edge_left)
    edge_l_inp = edge_left
    call RANDOM_NUMBER(edge_right)
    edge_r_inp = edge_right
    monotone = monotone_none

    call subgrid_quadratic_recon(reconstruction, &
                                 dep,            &
                                 field,          &
                                 edge_left,      &
                                 edge_right,     &
                                 monotone,       &
                                 order,          &
                                 nlayers)

    recon_inner_prod = dot_product(reconstruction, reconstruction)
    field_inner_prod = 0.0_r_tran
    do i_order = 1, 2*order + 1
      field_inner_prod = field_inner_prod &
                       + dot_product(field(:, i_order), field(:, i_order))
    end do
    edge_l_inner_prod = dot_product(edge_left, edge_left)
    edge_r_inner_prod = dot_product(edge_right, edge_right)

    write( log_scratch_space, * ) "adjt_subgrid_quadratic_recon inner products:"
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "reconstruction inner product = ", recon_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "field inner product = ", field_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "edge_left inner product = ", edge_l_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "edge_right inner product = ", edge_r_inner_prod
    call log_event( log_scratch_space, log_level_debug )

    recon_sf = 1.0_r_tran / real(recon_inner_prod + EPS, kind = r_tran)
    field_sf = 1.0_r_tran / real(field_inner_prod + EPS, kind = r_tran)
    edge_l_sf = 1.0_r_tran / real(edge_l_inner_prod + EPS, kind = r_tran)
    edge_r_sf = 1.0_r_tran / real(edge_r_inner_prod + EPS, kind = r_tran)

    inner1 = 0.0_r_tran
    inner1 = inner1 + recon_inner_prod * recon_sf
    inner1 = inner1 + field_inner_prod * field_sf
    inner1 = inner1 + edge_l_inner_prod * edge_l_sf
    inner1 = inner1 + edge_r_inner_prod * edge_r_sf

    reconstruction(:) = reconstruction(:) * recon_sf
    field(:,:) = field(:,:) * field_sf
    edge_left(:) = edge_left(:) * edge_l_sf
    edge_right(:) = edge_right(:) * edge_r_sf

    call adj_subgrid_quadratic_recon(reconstruction, &
                                     dep,            &
                                     field,          &
                                     edge_left,      &
                                     edge_right,     &
                                     monotone,       &
                                     order,          &
                                     nlayers)

    recon_recon_inp_inner_prod = dot_product(reconstruction, recon_inp)
    fld_fld_inp_inner_prod = 0.0_r_tran
    do i_order = 1, 2*order + 1
      fld_fld_inp_inner_prod = fld_fld_inp_inner_prod &
                             + dot_product(field(:, i_order), fld_inp(:, i_order))
    end do
    edge_l_edge_l_inp_inner_prod = dot_product(edge_left, edge_l_inp)
    edge_r_edge_r_inp_inner_prod = dot_product(edge_right, edge_r_inp)

    inner2 = 0.0_r_tran
    inner2 = inner2 + recon_recon_inp_inner_prod
    inner2 = inner2 + fld_fld_inp_inner_prod
    inner2 = inner2 + edge_l_edge_l_inp_inner_prod
    inner2 = inner2 + edge_r_edge_r_inp_inner_prod

    deallocate(reconstruction, &
               dep,            &
               field,          &
               edge_left,      &
               edge_right,     &
               recon_inp,      &
               fld_inp,        &
               edge_l_inp,     &
               edge_r_inp)

    ! Test the inner-product values for equality, allowing for the precision of the active variables
    machine_tol = spacing(max(abs(inner1), abs(inner2)))
    relative_diff = abs(inner1 - inner2) / machine_tol
    if (relative_diff < overall_tolerance) then
      write(log_scratch_space, *) "PASSED subgrid_quadratic_recon:", inner1, inner2, relative_diff
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
    else
      write(log_scratch_space, *) "FAILED subgrid_quadratic_recon:", inner1, inner2, relative_diff
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if


  end subroutine adjt_subgrid_quadratic_recon

end module adjt_subgrid_common_support_mod
