!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint vertical mass flux using an
!!        upwind reconstruction of a field on cell edges
module adj_vertical_mass_flux_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_REAL,           &
                                GH_READINC, GH_READWRITE,    &
                                GH_READ,                     &
                                CELL_COLUMN,                 &
                                ANY_DISCONTINUOUS_SPACE_1,   &
                                ANY_W2
  use constants_mod,     only : r_tran, i_def, l_def
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: adj_vertical_mass_flux_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                      &
         arg_type(GH_FIELD,  GH_REAL, GH_READINC,   ANY_W2), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_W2), &
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1) &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_vertical_mass_flux_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: adj_vertical_mass_flux_code

  contains

  !> @brief Computes the adjoint vertical mass flux.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] mass_flux      Vertical mass flux
  !> @param[in]     wind           Wind field
  !> @param[in,out] reconstruction Tracer field reconstructed on cell edges
  !> @param[in]     ndf_w2         Number of degrees of freedom per cell
  !> @param[in]     undf_w2        Number of unique degrees of freedom for the wind field
  !> @param[in]     map_w2         Dofmap for the cell at the base of the column
  !> @param[in]     ndf_md         Number of degrees of freedom per cell
  !> @param[in]     undf_md        Number of unique degrees of freedom for the
  !!                               reconstructed field
  !> @param[in]     map_md         Dofmap for the cell at the base of the column
  subroutine adj_vertical_mass_flux_code( nlayers,                &
                                          mass_flux,              &
                                          wind,                   &
                                          reconstruction,         &
                                          ndf_w2,                 &
                                          undf_w2,                &
                                          map_w2,                 &
                                          ndf_md,                 &
                                          undf_md,                &
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

    real(kind=r_tran), dimension(undf_md), intent(inout) :: reconstruction
    real(kind=r_tran), dimension(undf_w2), intent(in)    :: wind
    real(kind=r_tran), dimension(undf_w2), intent(inout) :: mass_flux

    ! Internal variables
    integer(kind=i_def) :: k, df, ijp1, ijp2, offset
    real(kind=r_tran)   :: wgt

    ! df of dof on face
    if  ( ndf_w2 == 2 ) then
      ! W2v space, reconstruction has ndata=2
      df = 1
      offset = 0
    else
      ! W2 space, reconstruction has ndata=6
      df = 5
      offset = 4*nlayers
    end if

    mass_flux( map_w2(df) + nlayers ) = 0.0_r_tran
    do k = nlayers - 1, 1, -1
      ! Value from bottom face of cell above edge
      ijp2 = map_md(1) + offset + nlayers + k-1
      ! Take value from top face of cell below edge
      ijp1 = map_md(1) + offset + k

      wgt = (0.5_r_tran - sign(0.5_r_tran, wind(map_w2(df)+k)))

      reconstruction(ijp1) = reconstruction(ijp1) &
                           + wind(map_w2(df)+k)*wgt*mass_flux(map_w2(df) + k)
      reconstruction(ijp2) = reconstruction(ijp2) &
                           + wind(map_w2(df)+k)*(1.0_r_tran-wgt)*mass_flux(map_w2(df) + k)
      mass_flux(map_w2(df) + k) = 0.0_r_tran

    end do
    mass_flux( map_w2(df) ) = 0.0_r_tran

  end subroutine adj_vertical_mass_flux_code

end module adj_vertical_mass_flux_kernel_mod
