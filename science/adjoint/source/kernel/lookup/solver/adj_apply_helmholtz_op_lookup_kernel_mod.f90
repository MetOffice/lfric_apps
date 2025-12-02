!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Apply adjoint Helmholtz operator stored as stencil of coefficients.

module adj_apply_helmholtz_op_lookup_kernel_mod

  use argument_mod,          only: arg_type,              &
                                   GH_FIELD, GH_REAL,     &
                                   GH_INTEGER, &
                                   GH_SCALAR, GH_LOGICAL, &
                                   GH_READ, GH_WRITE, GH_READWRITE,    &
                                   STENCIL, CROSS2D,      &
                                   CELL_COLUMN, &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,         only: r_solver, i_def, l_def
  use fs_continuity_mod,     only: W3
  use reference_element_mod, only: W, S, E, N
  use kernel_mod,            only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: adj_apply_helmholtz_op_lookup_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                    &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ, W3),                         &  ! y
         arg_type(GH_FIELD,   GH_REAL,    GH_READWRITE,  W3),                   &  ! x
         arg_type(GH_FIELD*9, GH_REAL,    GH_READ,  W3),                        &  ! Helmholtz
         arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &  ! lookup
         arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &  ! set_counts
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &  ! nindices
         arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)                              &  ! limited_area
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_apply_helmholtz_op_lookup_code
  end type adj_apply_helmholtz_op_lookup_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  public :: adj_apply_helmholtz_op_lookup_code

contains

!> @brief Apply the Helmholtz operator using an adjoint lookup table.
!> @param[in]     nlayers      Number of vertical levels to solve over
!> @param[in,out] y            Application of the operator to the pressure field
!> @param[in]     x            Pressure field to apply operator to
!> @param[in]     Helm_C       Diagonal entry to Helmholtz matrix
!> @param[in]     Helm_N       North (j+1) entry to Helmholtz matrix
!> @param[in]     Helm_E       East (i+1) entry to Helmholtz matrix
!> @param[in]     Helm_S       South (j-1) entry to Helmholtz matrix
!> @param[in]     Helm_W       West (j-1) entry to Helmholtz matrix
!> @param[in]     Helm_U       Upper (k+1) entry to Helmholtz matrix
!> @param[in]     Helm_UU      2nd Upper (k+2) entry to Helmholtz matrix
!> @param[in]     Helm_D       Lower (k-1) entry to Helmholtz matrix
!> @param[in]     Helm_DD      2nd Lower (k-2) entry to Helmholtz matrix
!> @param[in]     lookup       Adjoint lookup table.
!> @param[in]     set_counts   Adjoint lookup table index set counters.
!> @param[in]     nindices     Number of indices per index set in lookup table.
!> @param[in]     limited_area Switch to use code that can handle stencils at edges
!> @param[in]     ndf          Number of dofs per cell for all fields, should
!!                             be = 1
!> @param[in]     undf         Size of all field arrays
!> @param[in]     map          Array containing the address of the first dof in the
!> @param[in]     ndf_lu       Number of dofs per cell (lookup table)
!> @param[in]     undf_lu      Size of field arrays (lookup table)
!> @param[in]     map_lu       Lookup table dofmap.
!> @param[in]     ndf_sc       Number of dofs per cell (set counts)
!> @param[in]     undf_sc      Size of field arrays (set counts)
!> @param[in]     map_sc       Set counts dofmap.
subroutine adj_apply_helmholtz_op_lookup_code(nlayers, y, x, &
                                              Helm_C,         &
                                              Helm_N, Helm_E, Helm_S, Helm_W,   &
                                              Helm_U, Helm_UU, Helm_D, Helm_DD, &
                                              lookup, &
                                              set_counts, &
                                              nindices, &
                                              limited_area, &
                                              ndf, undf, map, &
                                              ndf_lu, undf_lu, map_lu, &
                                              ndf_sc, undf_sc, map_sc)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers, ndf, undf
  integer(kind=i_def), intent(in)    :: ndf_lu, undf_lu
  integer(kind=i_def), intent(in)    :: ndf_sc, undf_sc
  real(kind=r_solver), intent(inout) :: y(undf), x(undf)
  real(kind=r_solver), dimension(undf), intent(in) :: Helm_C,          &
                                                      Helm_N, Helm_E,  &
                                                      Helm_S, Helm_W,  &
                                                      Helm_U, Helm_UU, &
                                                      Helm_D, Helm_DD
  integer(kind=i_def), dimension(undf_lu), intent(in) :: lookup
  integer(kind=i_def), dimension(undf_sc), intent(in) :: set_counts
  integer(kind=l_def), intent(in) :: nindices
  logical(kind=l_def), intent(in) :: limited_area
  integer(kind=i_def), intent(in) :: map(ndf)
  integer(kind=i_def), intent(in) :: map_lu(ndf_lu)
  integer(kind=i_def), intent(in) :: map_sc(ndf_sc)

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def) :: ss, set_idx
  integer(kind=i_def) :: nsets
  integer(kind=i_def) :: neighbour
  integer(kind=i_def) :: recorded_face

  ! Note: global only (limited_area = .false.)
  ! This is checked for in the calling code

  do k = 2, nlayers - 1
    x(map(1) + k - 1) = x(map(1) + k - 1) + Helm_D(map(1) + k) * y(map(1) + k)
    x(map(1) + k - 2) = x(map(1) + k - 2) + Helm_DD(map(1) + k) * y(map(1) + k)
  end do

  k = 1
  x(map(1) + k - 1) = x(map(1) + k - 1) + Helm_D(map(1) + k) * y(map(1) + k)

  k = nlayers - 2
  x(map(1) + k + 1) = x(map(1) + k + 1) + Helm_U(map(1) + k) * y(map(1) + k)

  do k = 0, nlayers - 3
    x(map(1) + k + 1) = x(map(1) + k + 1) + Helm_U(map(1) + k) * y(map(1) + k)
    x(map(1) + k + 2) = x(map(1) + k + 2) + Helm_UU(map(1) + k) * y(map(1) + k)
  end do

  do k = 0, nlayers - 1
    x(map(1)+k) = x(map(1)+k) + Helm_C(map(1)+k)*y(map(1)+k)
  end do

  ! STENCIL operation
  ! Unpack the lookup table for this cell.
  nsets = set_counts( map_sc(1) )

  do ss = 1, nsets
    ! For each set of indices in the lookup table, apply the adjoint transformation.
    set_idx = map_lu(1) + (ss - 1)*nindices
    neighbour = lookup( set_idx )
    recorded_face = lookup( set_idx + 1 )  ! Face recorded as accessed by a neighbour.

    select case ( recorded_face )
      case ( W )
        do k = 0, nlayers - 1
          x(map(1)+k) = x(map(1)+k) + Helm_W(neighbour + k)*y(neighbour + k)
        end do
      case ( S )
        do k = 0, nlayers - 1
          x(map(1)+k) = x(map(1)+k) + Helm_S(neighbour + k)*y(neighbour + k)
        end do
      case ( E )
        do k = 0, nlayers - 1
          x(map(1)+k) = x(map(1)+k) + Helm_E(neighbour + k)*y(neighbour + k)
        end do
      case ( N )
        do k = 0, nlayers - 1
          x(map(1)+k) = x(map(1)+k) + Helm_N(neighbour + k)*y(neighbour + k)
        end do
    end select
  end do

end subroutine adj_apply_helmholtz_op_lookup_code

end module adj_apply_helmholtz_op_lookup_kernel_mod
