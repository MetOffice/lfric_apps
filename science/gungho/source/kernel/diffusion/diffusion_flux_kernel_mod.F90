!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!
!> @brief Computes the horizontal flux for the tracer diffusion equation
!> @details Computes the flux at W2H points for the diffusion of a tracer,
!!          using a linear reconstruction. For tracer m, flux F, viscosity nu,
!!          and density rho, the flux is given by:
!!          F = <nu>*<rho>*grad(m)
!!          where < > indicates averaging to W2H points, and grad(m) is given
!!          by a difference divided by the grid spacing.
module diffusion_flux_kernel_mod

  use argument_mod,       only : arg_type,            &
                                 GH_FIELD, GH_REAL,   &
                                 GH_INC, GH_READ,     &
                                 CELL_COLUMN
  use constants_mod,      only : r_def, i_def, EPS
  use fs_continuity_mod,  only : W2H, W3, W2
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: diffusion_flux_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                        &
        arg_type(GH_FIELD, GH_REAL, GH_INC,  W2H),                             &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                              &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W2H),                             &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W2H),                             &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                              &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                              &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W2H)                              &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: diffusion_flux_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: diffusion_flux_code

contains

  !> @brief Computes the horizontal flux for the tracer diffusion equation
  !> @param[in]     nlayers    The number of layers (in the shifted mesh)
  !> @param[in,out] flux       Flux field in shifted W2H to compute
  !> @param[in]     tracer     Tracer field in shifted W3 space
  !> @param[in]     visc       Viscosity field in shifted W2H space
  !> @param[in]     rho        Density field in shifted W2H space
  !> @param[in]     dx_at_w2   Distance between cell centres at W2 points
  !> @param[in]     dA_at_w2   Area of faces at W2 points
  !> @param[in]     rmult      Reciprocal of the multiplicity at W2H points
  !> @param[in]     ndf_w2h    Number of degrees of freedom per cell for W2h
  !> @param[in]     undf_w2h   Number of partition degrees of freedom for W2h
  !> @param[in]     map_w2h    Dofmap for W2h
  !> @param[in]     ndf_w3     Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3    Number of partition degrees of freedom for W3
  !> @param[in]     map_w3     Dofmap for W3
  !> @param[in]     ndf_w2     Number of degrees of freedom per cell for W2
  !> @param[in]     undf_w2    Number of partition degrees of freedom for W2
  !> @param[in]     map_w2     Dofmap for W2
  subroutine diffusion_flux_code(                                              &
      nlayers,                                                                 &
      flux, tracer, visc, rho, dx_at_w2, dA_at_w2, rmult,                      &
      ndf_w2h, undf_w2h, map_w2h,                                              &
      ndf_w3, undf_w3, map_w3,                                                 &
      ndf_w2, undf_w2, map_w2                                                  &
  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h
    integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
    integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
    real(kind=r_def),    intent(inout) :: flux(undf_w2h)
    real(kind=r_def),    intent(in)    :: tracer(undf_w3)
    real(kind=r_def),    intent(in)    :: visc(undf_w2h)
    real(kind=r_def),    intent(in)    :: rho(undf_w2h)
    real(kind=r_def),    intent(in)    :: dx_at_w2(undf_w2)
    real(kind=r_def),    intent(in)    :: dA_at_w2(undf_w2)
    real(kind=r_def),    intent(in)    :: rmult(undf_w2h)

    ! Encodes the direction of the gradient at the W2H points (W, S, E, N)
    ! Note that the y-direction is flipped
    integer(kind=i_def), parameter :: grad_sign(4) = (/ 1, -1, -1, 1 /)

    integer(kind=i_def) :: df
    integer(kind=i_def) :: w2h_b_idx, w2h_t_idx
    integer(kind=i_def) :: w2_b_idx,  w2_t_idx
    integer(kind=i_def) :: w3_b_idx,  w3_t_idx

    w3_b_idx = map_w3(1)+1
    w3_t_idx = map_w3(1)+nlayers-2

    do df = 1, ndf_w2h
      ! Extract array indices for bottom and top values to use in column
      ! Don't compute any values for top and bottom layers
      w2h_b_idx = map_w2h(df)+1
      w2h_t_idx = map_w2h(df)+nlayers-2
      w2_b_idx = map_w2(df)+1
      w2_t_idx = map_w2(df)+nlayers-2

      ! Only compute calculation when not at domain boundary
      if (rmult(w2h_b_idx) < 1.0_r_def - EPS) then
        flux(w2h_b_idx:w2h_t_idx) = flux(w2h_b_idx:w2h_t_idx)                  &
          + visc(w2h_b_idx:w2h_t_idx) * rho(w2h_b_idx:w2h_t_idx)               &
          * dA_at_w2(w2_b_idx:w2_t_idx) / dx_at_w2(w2_b_idx:w2_t_idx)          &
          * grad_sign(df) * tracer(w3_b_idx:w3_t_idx)
      end if
    end do


  end subroutine diffusion_flux_code

end module diffusion_flux_kernel_mod
