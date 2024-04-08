!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Samples a flux at nodal points: F = u*q.
!>
!> Samples a flux field F at nodal points where the flux is defined as the
!> product of a velocity u and a scalar q The quadrature rule is overloaded to
!> give the nodal points of the flux space.
!> As used in Multigrid, u = 1*1/2 = 1/2
!>
module mg_flux_kernel_mod

  use argument_mod,      only : arg_type, func_type,   &
                                GH_FIELD, GH_REAL,     &
                                GH_READ, GH_INC,       &
                                ANY_SPACE_1, GH_BASIS, &
                                CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only : r_def, r_solver, i_def
  use fs_continuity_mod, only : W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: mg_flux_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                    &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),         &
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_SPACE_1) &
         /)
    type(func_type) :: meta_funcs(1) = (/                  &
         func_type(ANY_SPACE_1, GH_BASIS)                  &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: mg_flux_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: mg_flux_code

contains

!> @brief Kernel to mg a flux at nodal points: F = u*q
!! @param[in] nlayers Number of layers
!! @param[in] ndf_f Number of degrees of freedom per cell for w2
!! @param[in] undf_f Number of unique degrees of freedom for w2
!! @param[in] map_f Dofmap for the cell at the base of the column for w2
!! @param[in,out] flux Field to contain the right hand side to be computed
!! @param[in] rmultiplicity Reciprocal of How many times the dof has been visited in total
!! @param[in] u Advecting wind
!! @param[in] ndf_q Number of degrees of freedom per cell for the field to be advected
!! @param[in] undf_q  Number of unique degrees of freedom for the advected field
!! @param[in] map_q Dofmap for the cell at the base of the column for the field to be advected
!! @param[in] basis_q Basis functions evaluated at gaussian quadrature points
!! @param[in] q Advected field
subroutine mg_flux_code(nlayers,                      &
                        flux, q ,                     &
                        ndf_f, undf_f, map_f,         &
                        ndf_q, undf_q, map_q, basis_q &
                        )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_f, ndf_q, undf_f, undf_q

  integer(kind=i_def), dimension(ndf_f), intent(in) :: map_f
  integer(kind=i_def), dimension(ndf_q), intent(in) :: map_q

  real(kind=r_def),    dimension(1,ndf_q,ndf_f), intent(in)    :: basis_q
  real(kind=r_solver), dimension(undf_f),        intent(inout) :: flux
  real(kind=r_solver), dimension(undf_q),        intent(in)    :: q

  ! Internal variables
  integer(kind=i_def)                :: df, df_q, k

  real(kind=r_solver), dimension(ndf_q) :: q_cell
  real(kind=r_solver)                   :: q_at_node

  do k = 0, nlayers-1
    do df_q = 1, ndf_q
      q_cell(df_q) = q( map_q(df_q) + k )
    end do
    do df = 1, ndf_f
      q_at_node = 0.0_r_solver
      do df_q = 1,ndf_q
        q_at_node = q_at_node + q_cell(df_q)*real(basis_q(1,df_q,df),r_solver)
      end do
      ! base and lid would have factor of 1 not 0.5 but field is zero there anyway!
      flux( map_f(df) + k ) = flux( map_f(df) + k ) + 0.5_r_solver*q_at_node
    end do
  end do

end subroutine mg_flux_code

end module mg_flux_kernel_mod
