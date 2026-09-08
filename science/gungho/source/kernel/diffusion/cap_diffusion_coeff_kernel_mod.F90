!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!
!> @brief Caps the diffusion coefficient
module cap_diffusion_coeff_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 GH_READWRITE, GH_READ, &
                                 CELL_COLUMN
  use constants_mod,      only : r_def, i_def
  use fs_continuity_mod,  only : Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: cap_diffusion_coeff_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                        &
        arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta),                     &
        arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta)                      &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: cap_diffusion_coeff_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: cap_diffusion_coeff_code

contains

  !> @brief Caps the diffusion coefficient to ensure stability
  !> @param[in]     nlayers    The number of layers
  !> @param[in,out] visc       Diffusion coefficient to be capped
  !> @param[in]     cap        Field in Wtheta to cap to
  !> @param[in]     ndf_wt     Number of degrees of freedom per cell for Wtheta
  !> @param[in]     undf_wt    Number of partition degrees of freedom for Wtheta
  !> @param[in]     map_wt     Dofmap for Wtheta
  subroutine cap_diffusion_coeff_code(                                         &
      nlayers,                                                                 &
      visc, cap,                                                               &
      ndf_wt, undf_wt, map_wt                                                  &
  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_wt, undf_wt
    integer(kind=i_def), intent(in)    :: map_wt(ndf_wt)
    real(kind=r_def),    intent(inout) :: visc(undf_wt)
    real(kind=r_def),    intent(in)    :: cap(undf_wt)

    integer(kind=i_def) :: b_idx,  t_idx

    b_idx = map_wt(1)
    t_idx = map_wt(1)+nlayers

    visc(b_idx:t_idx) = MIN(visc(b_idx:t_idx), cap(b_idx:t_idx))

  end subroutine cap_diffusion_coeff_code

end module cap_diffusion_coeff_kernel_mod
