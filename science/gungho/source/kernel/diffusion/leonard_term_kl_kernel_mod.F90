!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculate the Leonard term parameter, kl, needed for calculating the
!>        Leonard term fluxes
!>
module leonard_term_kl_kernel_mod

  use argument_mod,          only : arg_type,                     &
                                    GH_FIELD, GH_SCALAR, GH_REAL, &
                                    GH_READ, GH_WRITE,            &
                                    CELL_COLUMN, STENCIL, CROSS,  &
                                    GH_INTEGER
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta, W3
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.

  type, public, extends(kernel_type) :: leonard_term_kl_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                    &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W3),                       &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta, STENCIL(CROSS)),   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: leonard_term_kl_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: leonard_term_kl_code

contains

!> @brief Calculates kl on w3 levels accounting for stability limit
!! @param[in] nlayers  Number  of layers in the mesh
!! @param[in,out] kl  Stability-corrected Leonard term parameter
!! @param[in] velocity_w2v   velocity normal to cell top/bottom
!! @param[in] map_wt_stencil_size  Number of cells in the stencil at the base
!!                                 of the column for Wtheta
!! @param[in] map_wt_stencil  Array holding the dofmap for the stencil at the
!!                            base of the column for Wtheta
!! @param[in] height_wth  Height of wth space levels above the surface
!! @param[in] leonard_kl  The user-specified Leonard term parameter
!! @param[in] dt  The model timestep length
!! @param[in] bl_levels   The number of boundary-layer levels
!! @param[in] ndf_w3  Number of degrees of freedom per cell for w3 space
!! @param[in] undf_w3  Number of unique degrees of freedom for w3 space
!! @param[in] map_w3  Cell dofmap for w3 space
!! @param[in] ndf_wt  Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom for theta space
!! @param[in] map_wt  Cell dofmap for theta space
subroutine leonard_term_kl_code( nlayers,                               &
                                 kl,                                    &
                                 velocity_w2v,                          &
                                 map_wt_stencil_size, map_wt_stencil,   &
                                 height_wth,                            &
                                 leonard_kl,                            &
                                 dt, bl_levels,                         &
                                 ndf_w3, undf_w3, map_w3,               &
                                 ndf_wt, undf_wt, map_wt                &
                                )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, bl_levels
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: map_wt_stencil_size
  integer(kind=i_def), dimension(ndf_wt,map_wt_stencil_size), intent(in)  :: map_wt_stencil
  integer(kind=i_def), dimension(ndf_wt),  intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_w3),  intent(in)  :: map_w3

  real(kind=r_def), intent(in)    :: leonard_kl
  real(kind=r_def), intent(in)    :: dt
  real(kind=r_def), dimension(undf_w3),  intent(inout) :: kl
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: height_wth
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: velocity_w2v

  ! Internal variables
  integer(kind=i_def) :: k, kp

  ! If the full stencil isn't available, we must be at the domain edge.
  ! Simply set the increment to 0 for now, and exit the routine.
  if (map_wt_stencil_size < 5_i_def) then
    do k = 0, bl_levels
      kl(map_w3(1) + k) = leonard_kl
    end do
    return
  end if

  ! Leonard term parameter is the min of the input leonard_kl
  ! and the max stable value   6 * dz / ( dt * dw )
  ! For dw we use the maximum horizontal finite difference that
  ! contributes to the flux at each w3 point.

  ! Set to zero at k=0 as Leonard term is zero at k=0
  k = 0
  kl(map_w3(1) + k) = 0.0_r_def

  do k = 1, bl_levels - 1
    kp = k + 1
    kl(map_w3(1) + k) = MIN( leonard_kl,                                 &
                        6.0_r_def * ( height_wth(map_wt(1) + kp) -       &
                                      height_wth(map_wt(1) + k) )        &
                        / (dt * MAX(                                     &
                        ! Difference on wth level above
                        ABS( velocity_w2v(map_wt_stencil(1,2) + kp) -    &
                             velocity_w2v(map_wt_stencil(1,1) + kp) ),   &
                        ABS( velocity_w2v(map_wt_stencil(1,3) + kp) -    &
                             velocity_w2v(map_wt_stencil(1,1) + kp) ),   &
                        ABS( velocity_w2v(map_wt_stencil(1,4) + kp) -    &
                             velocity_w2v(map_wt_stencil(1,1) + kp) ),   &
                        ABS( velocity_w2v(map_wt_stencil(1,5) + kp) -    &
                             velocity_w2v(map_wt_stencil(1,1) + kp) ),   &
                        ! Difference on wth level below
                        ABS( velocity_w2v(map_wt_stencil(1,2) + k) -     &
                             velocity_w2v(map_wt_stencil(1,1) + k) ),    &
                        ABS( velocity_w2v(map_wt_stencil(1,3) + k) -     &
                             velocity_w2v(map_wt_stencil(1,1) + k) ),    &
                        ABS( velocity_w2v(map_wt_stencil(1,4) + k) -     &
                             velocity_w2v(map_wt_stencil(1,1) + k) ),    &
                        ABS( velocity_w2v(map_wt_stencil(1,5) + k) -     &
                             velocity_w2v(map_wt_stencil(1,1) + k) ),    &
                        EPSILON( leonard_kl )                            &
                        ) ) )
  end do

  ! Set to zero at k=bl_levels as Leonard term isn't defined at k=bl_levels
  k = bl_levels
  kl(map_w3(1) + k) = 0.0_r_def

end subroutine leonard_term_kl_code

end module leonard_term_kl_kernel_mod
