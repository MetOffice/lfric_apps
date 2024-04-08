!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes a mass increment for the bottom layer of each column
!> from a flux relating to a source or sink.
!>
!> @details The kernel computes the mass for each column from a flux,
!> the bottom area and the time-step.
!> \f[ M = F * A * dt \f]
!>
module compute_bottom_mass_from_flux_kernel_mod

  use argument_mod,            only : arg_type,               &
                                      GH_FIELD, GH_REAL,      &
                                      GH_WRITE, GH_READ,      &
                                      GH_SCALAR, CELL_COLUMN, &
                                      ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W2
  use kernel_mod,              only : kernel_type
  use reference_element_mod,   only : B

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: compute_bottom_mass_from_flux_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                     &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),                        &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: compute_bottom_mass_from_flux_code
  end type
  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_bottom_mass_from_flux_code
contains

!> @brief Compute the mass in a column corresponding to an incoming flux
!! @param[in] nlayers The number of layers in the mesh.
!! @param[in,out] mass  A 2D discontinuous field whose values in the bottom layer
!!                      will correspond to the flux
!! @param[in] flux    A 2D discontinuous field corresponding to a flux
!!                    on the bottom boundary
!! @param[in] area    A field in W2 whose values are the areas of each
!!                    face of each element
!! @param[in] dt      The model timestep length
!! @param[in] ndf_2d  The number of degrees of freedom per cell for 2D field
!! @param[in] undf_2d The number of unique degrees of freedom for 2D field
!! @param[in] map_2d  Dofmap for the cell at the base of the column
!!                    for the 2D field
!! @param[in] ndf_w2  The number of degrees of freedom per cell for W2
!! @param[in] undf_w2 The number of unique degrees of freedom for W2
!! @param[in] map_w2  Dofmap for the cell at the base of the column for w2
subroutine compute_bottom_mass_from_flux_code(                            &
                                               nlayers,                   &
                                               mass,                      &
                                               flux,                      &
                                               area,                      &
                                               dt,                        &
                                               ndf_2d, undf_2d, map_2d,   &
                                               ndf_w2, undf_w2, map_w2   &
                                             )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d, ndf_w2, undf_w2
  integer(kind=i_def), dimension(ndf_2d), intent(in)    :: map_2d
  integer(kind=i_def), dimension(ndf_w2), intent(in)    :: map_w2
  real(kind=r_def), dimension(undf_2d),   intent(inout) :: mass
  real(kind=r_def), dimension(undf_2d),   intent(in)    :: flux
  real(kind=r_def), dimension(undf_w2),   intent(in)    :: area
  real(kind=r_def),                       intent(in)    :: dt

  mass( map_2d(1) ) = flux( map_2d(1) ) * area( map_w2(B) ) * dt

end subroutine compute_bottom_mass_from_flux_code

end module compute_bottom_mass_from_flux_kernel_mod
