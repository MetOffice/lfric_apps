!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes v.Jv on nodal points for normalising w2 fields.
!>
!> Compute the normalsation factor fro W2 fields as vJv on W2 nodes.
!> For a regular orthogonal grid J = diag(dx,dy,dz) so:
!> vJv = (dx,0,0) for u components
!> vJv = (0,dy,0) for v components
!> vJv = (0,0,dz) for w components
!>
!> @todo Create unit test for this kernel, see #2935
module w2_normalisation_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD, GH_REAL,         &
                                GH_INC, GH_READ,           &
                                ANY_SPACE_9,               &
                                ANY_DISCONTINUOUS_SPACE_3, &
                                GH_BASIS, GH_DIFF_BASIS,   &
                                CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only : r_def, r_solver, i_def
  use fs_continuity_mod, only : W2, Wchi
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: w2_normalisation_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                        &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, Wchi),                      &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2)                         &

         /)
    type(func_type) :: meta_funcs(2) = (/         &
         func_type(W2,   GH_BASIS),               &
         func_type(Wchi, GH_BASIS, GH_DIFF_BASIS) &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: w2_normalisation_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: w2_normalisation_code

contains

!> @brief Compute the normalisation factor for W2 fields as vJv
!! @param[in]     nlayers         Number of layers
!! @param[in,out] normalisation   Normalisation field to compute
!! @param[in]     chi_1           1st coordinate field in Wchi
!! @param[in]     chi_2           2nd coordinate field in Wchi
!! @param[in]     chi_3           3rd coordinate field in Wchi
!! @param[in]     panel_id        Field giving the ID for mesh panels
!! @param[in]     w2_multiplicity Dof multiplicity of W2 space
!! @param[in]     ndf             Number of degrees of freedom per cell
!! @param[in]     undf            Total number of degrees of freedom
!! @param[in]     map             Dofmap for the cell at the base of the column
!! @param[in]     basis           Basis functions evaluated at W2 nodal points
!! @param[in]     ndf_chi         Number of dofs per cell for the coordinate field
!! @param[in]     undf_chi        Total number of degrees of freedom
!! @param[in]     map_chi         Dofmap for the coordinate field
!! @param[in]     chi_basis       Wchi basis functions evaluated at W2 nodal points
!! @param[in]     chi_diff_basis  Wchi diff basis functions evaluated at W2 nodal points
!! @param[in]     ndf_pid         Number of degrees of freedom per cell for panel_id
!! @param[in]     undf_pid        Number of unique degrees of freedom for panel_id
!! @param[in]     map_pid         Dofmap for the cell at the base of the column for panel_id
subroutine w2_normalisation_code(nlayers,                 &
                                 normalisation,           &
                                 chi_1, chi_2, chi_3,     &
                                 panel_id,                &
                                 w2_multiplicity,         &
                                 ndf, undf,               &
                                 map, basis,              &
                                 ndf_chi, undf_chi,       &
                                 map_chi,                 &
                                 chi_basis,               &
                                 chi_diff_basis,          &
                                 ndf_pid, undf_pid,       &
                                 map_pid                  &
                                 )

  use sci_coordinate_jacobian_mod, only : coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_chi
  integer(kind=i_def), intent(in) :: undf, undf_chi
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid

  integer(kind=i_def), dimension(ndf),     intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(3,ndf_chi,ndf) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,ndf) :: chi_basis
  real(kind=r_def), intent(in), dimension(3,ndf,ndf)     :: basis

  real(kind=r_def), dimension(undf),     intent(inout) :: normalisation
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id
  real(kind=r_def), dimension(undf),     intent(in)    :: w2_multiplicity

  ! Internal variables
  integer(kind=i_def)                    :: df, k
  integer(kind=i_def)                    :: ipanel
  real(kind=r_def), dimension(ndf)       :: dj
  real(kind=r_def), dimension(3,3,ndf)   :: jacobian
  real(kind=r_def), dimension(ndf_chi)   :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def), dimension(3)         :: Jv
  real(kind=r_def), dimension(3,3)       :: JTJ
  real(kind=r_def)                       :: norm

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi_1_cell(df) = chi_1(map_chi(df) + k)
      chi_2_cell(df) = chi_2(map_chi(df) + k)
      chi_3_cell(df) = chi_3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, &
                             ndf, &
                             chi_1_cell, &
                             chi_2_cell, &
                             chi_3_cell, &
                             ipanel,     &
                             chi_basis,  &
                             chi_diff_basis, &
                             jacobian, &
                             dj)
    do df = 1,ndf
      JTJ =  matmul(transpose(jacobian(:,:,df)),jacobian(:,:,df))
      Jv = matmul(JTJ,basis(:,df,df))
      norm = sqrt(dot_product(Jv,basis(:,df,df)))/w2_multiplicity(map(df)+k)
      normalisation(map(df)+k) = normalisation(map(df)+k) + norm
    end do
  end do

end subroutine w2_normalisation_code

end module w2_normalisation_kernel_mod
