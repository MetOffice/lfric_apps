!------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief Module containing adjoint test for adj_subgrid_vertical_support_mod routines
module adjt_subgrid_vertical_support_mod

  use constants_mod,                  only : i_def, r_tran,     &
                                             l_def, EPS
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_DEBUG

  implicit none

  public :: adjt_third_order_vertical_edge

  contains
  !=============================================================================
  !> @brief   Adjoint test for adj_third_order_vertical_edge.
  !> @details Passes if adjoint is transpose of tangent linear.
  !>          Determined by testing the equality of inner products <Mx, Mx> and <AMx, x>,
  !>          where M is the tangent linear and A is the adjoint.
  !> @param[in] nlayers     Number of layers

  subroutine adjt_third_order_vertical_edge(nlayers)

    use subgrid_vertical_support_mod,     only : third_order_vertical_edge
    use adj_subgrid_vertical_support_mod, only : adj_third_order_vertical_edge

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers

    ! Arguments for tl and adj calls
    real(kind=r_tran), dimension(:), allocatable :: field
    real(kind=r_tran), dimension(:), allocatable :: dz
    real(kind=r_tran), dimension(:), allocatable :: dla_dz_1
    real(kind=r_tran), dimension(:), allocatable :: dla_dz_2
    real(kind=r_tran), dimension(:), allocatable :: dla_dz_3
    real(kind=r_tran), dimension(:), allocatable :: dlb_dz_1
    real(kind=r_tran), dimension(:), allocatable :: dlb_dz_2
    real(kind=r_tran), dimension(:), allocatable :: dlb_dz_3
    real(kind=r_tran), dimension(:), allocatable :: edge_above
    real(kind=r_tran), dimension(:), allocatable :: edge_below
    logical(kind=l_def)                          :: log_space

    ! Copies of input fields used in inner products
    real(kind=r_tran), dimension(:), allocatable :: fld_inp
    real(kind=r_tran), dimension(:), allocatable :: edge_a_inp
    real(kind=r_tran), dimension(:), allocatable :: edge_b_inp

    ! Inner products
    real(kind=r_tran)   :: field_inner_prod
    real(kind=r_tran)   :: edge_a_inner_prod
    real(kind=r_tran)   :: edge_b_inner_prod
    real(kind=r_tran)   :: field_sf
    real(kind=r_tran)   :: edge_a_sf
    real(kind=r_tran)   :: edge_b_sf
    real(kind=r_tran)   :: inner1
    real(kind=r_tran)   :: fld_fld_inp_inner_prod
    real(kind=r_tran)   :: edge_a_edge_a_inp_inner_prod
    real(kind=r_tran)   :: edge_b_edge_b_inp_inner_prod
    real(kind=r_tran)   :: inner2

    ! Test parameters and variables
    real(kind=r_tran), parameter :: overall_tolerance = 1500.0_r_tran
    real(kind=r_tran)            :: machine_tol
    real(kind=r_tran)            :: relative_diff

    allocate(field(1:nlayers),      &
             dz(1:nlayers),         &
             dla_dz_1(1:nlayers),   &
             dla_dz_2(1:nlayers),   &
             dla_dz_3(1:nlayers),   &
             dlb_dz_1(1:nlayers),   &
             dlb_dz_2(1:nlayers),   &
             dlb_dz_3(1:nlayers),   &
             edge_above(1:nlayers), &
             edge_below(1:nlayers), &
             fld_inp(1:nlayers),    &
             edge_a_inp(1:nlayers), &
             edge_b_inp(1:nlayers))

    ! Initialise arguments and call the forward routine.
    call RANDOM_NUMBER(field)
    fld_inp = field
    call RANDOM_NUMBER(dz)
    call RANDOM_NUMBER(dla_dz_1)
    call RANDOM_NUMBER(dla_dz_2)
    call RANDOM_NUMBER(dla_dz_3)
    call RANDOM_NUMBER(dlb_dz_1)
    call RANDOM_NUMBER(dlb_dz_2)
    call RANDOM_NUMBER(dlb_dz_3)
    call RANDOM_NUMBER(edge_above)
    edge_a_inp = edge_above
    call RANDOM_NUMBER(edge_below)
    edge_b_inp = edge_below
    log_space = .false._l_def

    call third_order_vertical_edge(field,                        &
                                   dla_dz_1, dla_dz_2, dla_dz_3, &
                                   dlb_dz_1, dlb_dz_2, dlb_dz_3, &
                                   dz,                           &
                                   edge_above, edge_below,       &
                                   log_space, nlayers)


    field_inner_prod =  dot_product(field, field)
    edge_a_inner_prod = dot_product(edge_above, edge_above)
    edge_b_inner_prod = dot_product(edge_below, edge_below)

    write( log_scratch_space, * ) "adjt_third_order_vertical_edge inner products:"
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "field inner product = ", field_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "edge_above inner product = ", edge_a_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "edge_below inner product = ", edge_b_inner_prod
    call log_event( log_scratch_space, log_level_debug )

    field_sf = 1.0_r_tran / real(field_inner_prod + EPS, kind = r_tran)
    edge_a_sf = 1.0_r_tran / real(edge_a_inner_prod + EPS, kind = r_tran)
    edge_b_sf = 1.0_r_tran / real(edge_b_inner_prod + EPS, kind = r_tran)

    inner1 = 0.0_r_tran
    inner1 = inner1 + field_inner_prod * field_sf
    inner1 = inner1 + edge_a_inner_prod * edge_a_sf
    inner1 = inner1 + edge_b_inner_prod * edge_b_sf

    field(:) = field(:) * field_sf
    edge_above(:) = edge_above(:) * edge_a_sf
    edge_below(:) = edge_below(:) * edge_b_sf

    call adj_third_order_vertical_edge(field,                        &
                                       dla_dz_1, dla_dz_2, dla_dz_3, &
                                       dlb_dz_1, dlb_dz_2, dlb_dz_3, &
                                       dz,                           &
                                       edge_above, edge_below,       &
                                       log_space, nlayers)


    fld_fld_inp_inner_prod = dot_product(field, fld_inp)
    edge_a_edge_a_inp_inner_prod = dot_product(edge_above, edge_a_inp)
    edge_b_edge_b_inp_inner_prod = dot_product(edge_below, edge_b_inp)

    inner2 = 0.0_r_tran
    inner2 = inner2 + fld_fld_inp_inner_prod
    inner2 = inner2 + edge_a_edge_a_inp_inner_prod
    inner2 = inner2 + edge_b_edge_b_inp_inner_prod

    deallocate(field,      &
               dz,         &
               dla_dz_1,   &
               dla_dz_2,   &
               dla_dz_3,   &
               dlb_dz_1,   &
               dlb_dz_2,   &
               dlb_dz_3,   &
               edge_above, &
               edge_below, &
               fld_inp,    &
               edge_a_inp, &
               edge_b_inp)

    ! Test the inner-product values for equality, allowing for the precision of the active variables
    machine_tol = spacing(max(abs(inner1), abs(inner2)))
    relative_diff = abs(inner1 - inner2) / machine_tol
    if (relative_diff < overall_tolerance) then
      write(log_scratch_space, *) "PASSED third_order_vertical_edge:", inner1, inner2, relative_diff
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
    else
      write(log_scratch_space, *) "FAILED third_order_vertical_edge:", inner1, inner2, relative_diff
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if


  end subroutine adjt_third_order_vertical_edge

end module adjt_subgrid_vertical_support_mod
