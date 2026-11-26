!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of sampling of a flux at nodal points.
!>
module adj_sample_flux_kernel_mod

  use, intrinsic :: iso_fortran_env, only: real32, real64

  use argument_mod,      only : arg_type, func_type,     &
                                GH_FIELD, GH_REAL,       &
                                GH_READ, GH_INC,         &
                                ANY_SPACE_1,             &
                                GH_BASIS, GH_DIFF_BASIS, &
                                CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W0, W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: adj_sample_flux_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),         &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),         &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),         &
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_SPACE_1) &
         /)
    type(func_type) :: meta_funcs(1) = (/                  &
         func_type(ANY_SPACE_1, GH_BASIS)                  &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_sample_flux_code

  ! Generic interface for real32 and real64 types
  interface adj_sample_flux_code
    module procedure  &
      adj_sample_flux_code_real32, &
      adj_sample_flux_code_real64
  end interface
contains

!> @brief Adjoint of kernel that samples a flux at nodal points
!! @param[in] nlayers Number of layers
!! @param[in] flux Field to contain the right hand side to be computed
!! @param[in,out] u Advecting wind
!! @param[in] rmultiplicity Reciprocal of How many times the dof has been visited in total
!! @param[in] q Advected field
!! @param[in] ndf_f Number of degrees of freedom per cell for w2
!! @param[in] undf_f Number of unique degrees of freedom for w2
!! @param[in] map_f Dofmap for the cell at the base of the column for w2
!! @param[in] ndf_q Number of degrees of freedom per cell for the field to be advected
!! @param[in] undf_q  Number of unique degrees of freedom for the advected field
!! @param[in] map_q Dofmap for the cell at the base of the column for the field to be advected
!! @param[in] basis_q Basis functions evaluated at gaussian quadrature points

! REAL32 PRECISION
! ==================
subroutine adj_sample_flux_code_real32( nlayers,                      &
                                        flux, u, rmultiplicity, q,    &
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
  real(kind=real32),   dimension(undf_f),        intent(in)    :: flux, rmultiplicity
  real(kind=real32),   dimension(undf_f),        intent(inout) :: u
  real(kind=real32),   dimension(undf_q),        intent(in)    :: q

  ! Internal variables
  integer(kind=i_def)                   :: df, df_q, k, loc
  real(kind=real32),   dimension(ndf_q) :: q_cell
  real(kind=real32)                     :: q_at_node

  do k = nlayers-1, 0, -1
    do df_q = 1, ndf_q
      q_cell(df_q) = q( map_q(df_q) + k )
    end do
    do df = ndf_f, 1, -1
      q_at_node = 0.0_real32
      do df_q = 1, ndf_q
        q_at_node = q_at_node + q_cell(df_q)*real(basis_q(1,df_q,df), real32)
      end do
      loc = map_f(df) + k
      u(loc) = u(loc) + q_at_node * rmultiplicity(loc) * flux(loc)
    end do
  end do

end subroutine adj_sample_flux_code_real32

! REAL64 PRECISION
! ==================
subroutine adj_sample_flux_code_real64( nlayers,                      &
                                        flux, u, rmultiplicity, q,    &
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
  real(kind=real64),   dimension(undf_f),        intent(in)    :: flux, rmultiplicity
  real(kind=real64),   dimension(undf_f),        intent(inout) :: u
  real(kind=real64),   dimension(undf_q),        intent(in)    :: q

  ! Internal variables
  integer(kind=i_def)                   :: df, df_q, k, loc
  real(kind=real64),   dimension(ndf_q) :: q_cell
  real(kind=real64)                     :: q_at_node

  do k = nlayers-1, 0, -1
    do df_q = 1, ndf_q
      q_cell(df_q) = q( map_q(df_q) + k )
    end do
    do df = ndf_f, 1, -1
      q_at_node = 0.0_real64
      do df_q = 1, ndf_q
        q_at_node = q_at_node + q_cell(df_q)*real(basis_q(1,df_q,df), real64)
      end do
      loc = map_f(df) + k
      u(loc) = u(loc) + q_at_node * rmultiplicity(loc) * flux(loc)
    end do
  end do

end subroutine adj_sample_flux_code_real64

end module adj_sample_flux_kernel_mod
