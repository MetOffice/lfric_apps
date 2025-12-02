!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Applies horizontal Smagorinsky diffusion visc_h * (d2dx2 + d2dy2) to
!>        a tracer variable in the Wtheta space for lowest order elements.
!>
module tracer_smagorinsky_diff_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_WRITE,         &
                                    CELL_COLUMN, STENCIL, CROSS
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta, W2
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.

  type, public, extends(kernel_type) :: tracer_smagorinsky_diff_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta, STENCIL(CROSS)),   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W2)                        &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tracer_smagorinsky_diff_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: tracer_smagorinsky_diff_code

contains

!> @brief Calculates horizontal Smagorinsky diffusion for a tracer variable
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] theta_inc Diffusion increment for temperature field
!! @param[in] theta_n Input temperature field
!! @param[in] map_wt_stencil_size Number of cells in the stencil at the base
!!                                of the column for Wtheta
!! @param[in] map_wt_stencil Array holding the dofmap for the stencil at the
!!                           base of the column for Wtheta
!! @param[in] visc_h Diffusion coefficient for scalars on Wtheta points
!! @param[in] dx_at_w2 Grid length at cell faces
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom for theta space
!! @param[in] map_wt Cell dofmap for theta space
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2 space
!! @param[in] undf_w2  Number of unique degrees of freedom for w2 space
!! @param[in] map_w2 Cell dofmap for w2 space
subroutine tracer_smagorinsky_diff_code( nlayers,                              &
                                         theta_inc,                            &
                                         theta_n,                              &
                                         map_wt_stencil_size, map_wt_stencil,  &
                                         visc_h,                               &
                                         dx_at_w2,                             &
                                         ndf_wt, undf_wt, map_wt,              &
                                         ndf_w2, undf_w2, map_w2               &
                                        )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: map_wt_stencil_size
  integer(kind=i_def), dimension(ndf_wt,map_wt_stencil_size), intent(in)  :: map_wt_stencil
  integer(kind=i_def), dimension(ndf_wt),  intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2

  real(kind=r_def), dimension(undf_wt),  intent(inout) :: theta_inc
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: theta_n
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: visc_h
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: dx_at_w2

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_def)    :: d2dx, d2dy
  real(kind=r_def), dimension(0:nlayers-1) :: idx2, idy2

  ! Assumed direction for derivatives in this kernel is:
  !  y
  !  ^
  !  |_> x
  !

  ! The layout of the cells in the stencil is:
  !
  !          -----
  !          |   |
  !          | 5 |
  !     ---------------
  !     |    |   |    |
  !     |  2 | 1 |  4 |
  !     ---------------
  !          |   |
  !          | 3 |
  !          -----

  ! If the full stencil isn't available, we must be at the domain edge.
  ! Simply set the increment to 0 for now, and exit the routine.
  if (map_wt_stencil_size < 5_i_def) then
    do k = 0, nlayers
      theta_inc(map_wt(1) + k) = 0.0_r_def
    end do
    return
  end if

  ! Compute horizontal grid spacing - average of cell face values
  ! N.B. not accounting for difference between w3 and wtheta heights
  do k = 1, nlayers - 1
    idx2(k) = (2.0_r_def/(dx_at_w2(map_w2(1)+k)+dx_at_w2(map_w2(3)+k)))**2
    idy2(k) = (2.0_r_def/(dx_at_w2(map_w2(2)+k)+dx_at_w2(map_w2(4)+k)))**2
  end do

  ! Horizontal theta diffusion
  ! Set to zero at k=0 as shear(k=0) isn't defined
  k = 0
  theta_inc(map_wt(1) + k) = 0.0_r_def

  do k = 1, nlayers - 1
    d2dx = (theta_n(map_wt_stencil(1,2) + k)  - 2.0_r_def*theta_n(map_wt_stencil(1,1) + k) +      &
            theta_n(map_wt_stencil(1,4) + k) ) * idx2(k)
    d2dy = (theta_n(map_wt_stencil(1,3) + k)  - 2.0_r_def*theta_n(map_wt_stencil(1,1) + k) +      &
            theta_n(map_wt_stencil(1,5) + k) ) * idy2(k)
    theta_inc(map_wt(1) + k) = visc_h(map_wt(1) + k) * (d2dx + d2dy)
  end do

  k = nlayers
  theta_inc(map_wt(1) + k) = 0.0_r_def

end subroutine tracer_smagorinsky_diff_code

end module tracer_smagorinsky_diff_kernel_mod
