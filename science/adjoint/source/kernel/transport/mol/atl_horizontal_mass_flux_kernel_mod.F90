!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint horizontal mass flux of a density by a wind
!!        field.
module atl_horizontal_mass_flux_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD,                  &
                                GH_REAL,                   &
                                GH_INC, GH_READ,           &
                                CELL_COLUMN,               &
                                ANY_DISCONTINUOUS_SPACE_1, &
                                ANY_W2
  use constants_mod,     only : r_def, i_def
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: atl_horizontal_mass_flux_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, ANY_W2),                   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, ANY_W2),                   &
         arg_type(GH_FIELD,  GH_REAL, GH_INC,  ANY_W2),                   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1) &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: atl_horizontal_mass_flux_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: atl_horizontal_mass_flux_code

  contains

  !> @brief Computes the adjoint horizontal mass flux of a density by a wind
  !> @param[in]     nlayers        Number of layers
  !> @param[in]     mass_flux      Upwind mass flux to be computed
  !> @param[in]     ls_wind        Linear wind field
  !> @param[in,out] wind           Perturbation wind field
  !> @param[in]     reconstruction Density field computed on cell edges
  !> @param[in]     ndf_w2         Number of degrees of freedom per cell for the wind and flux fields
  !> @param[in]     undf_w2        Number of unique degrees of freedom for the wind and flux fields
  !> @param[in]     map_w2         Dofmap for the cell at the base of the column for the wind and flux fields
  !> @param[in]     ndf_md         Number of degrees of freedom per cell for the multidata density field
  !> @param[in]     undf_md        Number of unique degrees of freedom for the multidata density field
  !> @param[in]     map_md         Dofmap for the cell at the base of the column for the multidata density field
  subroutine atl_horizontal_mass_flux_code( nlayers,        &
                                            mass_flux,      &
                                            ls_wind,        &
                                            wind,           &
                                            reconstruction, &
                                            ndf_w2,         &
                                            undf_w2,        &
                                            map_w2,         &
                                            ndf_md,         &
                                            undf_md,        &
                                            map_md )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)                    :: nlayers
    integer(kind=i_def), intent(in)                    :: ndf_md
    integer(kind=i_def), intent(in)                    :: undf_md
    integer(kind=i_def), dimension(ndf_md), intent(in) :: map_md
    integer(kind=i_def), intent(in)                    :: ndf_w2
    integer(kind=i_def), intent(in)                    :: undf_w2
    integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

    real(kind=r_def), dimension(undf_md), intent(in)    :: reconstruction
    real(kind=r_def), dimension(undf_w2), intent(inout) :: wind
    real(kind=r_def), dimension(undf_w2), intent(in)    :: ls_wind
    real(kind=r_def), dimension(undf_w2), intent(in)    :: mass_flux

    ! Internal variables
    integer(kind=i_def)                    :: k, df, ijp
    integer(kind=i_def), parameter         :: nfaces = 4
    real(kind=r_def)                       :: direction
    real(kind=r_def), dimension(nfaces)    :: v_dot_n

    ! Implied direction of outward normals dotted with basis functions.
    ! If u*u_dot_n > 0 then this is the upwind cell
    v_dot_n = (/ -1.0_r_def, 1.0_r_def, 1.0_r_def, -1.0_r_def /)

    do df = nfaces, 1, -1
      do k = nlayers - 1, 0, -1
        direction = ls_wind(map_w2(df) + k)*v_dot_n(df)
        if ( direction > 0.0_r_def ) then
          ! Take value on edge from this column
          ijp = map_md(1) + (df-1)*nlayers
          wind(map_w2(df) + k) = wind(map_w2(df) + k) + &
                                 reconstruction(ijp + k) * mass_flux(map_w2(df) + k)
        end if
      end do
    end do

  end subroutine atl_horizontal_mass_flux_code

end module atl_horizontal_mass_flux_kernel_mod
