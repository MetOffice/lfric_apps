!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculates Leonard term vertical fluxes and increments for any
!>        field defined at theta-points.
!>
module leonard_term_th_kernel_mod

  use argument_mod,          only : arg_type,                     &
                                    GH_FIELD, GH_SCALAR, GH_REAL, &
                                    GH_READ, GH_WRITE,            &
                                    CELL_COLUMN, STENCIL, CROSS,  &
                                    GH_INTEGER
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta, W3
  use kernel_mod,            only : kernel_type
  use mixing_config_mod,     only : leonard_kl

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.

  type, public, extends(kernel_type) :: leonard_term_th_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                                    &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta, STENCIL(CROSS)),   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta, STENCIL(CROSS)),   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                       &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                       &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                       &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: leonard_term_th_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: leonard_term_th_code

contains

!> @brief Calculates increment due to vertical Leonard term flux
!>        for variables held in wt space
!! @param[in] nlayers  Number  of layers in the mesh
!! @param[in,out] field_inc  Leonard term increment on wtheta levels
!! @param[in] field  Field on wtheta levels
!! @param[in] map_wth_stencil_size  Number of cells in the stencil at the base
!!                                 of the column for wtheta
!! @param[in] map_wth_stencil  Array holding the dofmap for the stencil at the
!!                            base of the column for wtheta
!! @param[in] velocity_w2v  velocity normal to cell top/bottom
!! @param[in] map_wt_stencil_size  Number of cells in the stencil at the base
!!                                 of the column for Wtheta
!! @param[in] map_wt_stencil  Array holding the dofmap for the stencil at the
!!                            base of the column for Wtheta
!! @param[in] dtrdz_tq_bl  Array of timestep / ( dz * rho * r^2 ),
!!                         precalculated in the BL scheme
!! @param[in] kl  Leonard term parameter on w3 space levels
!! @param[in] rho  Density on w3 space levels
!! @param[in] height_w3  Height of w3 space levels above the surface
!! @param[in] planet_radius The planet radius
!! @param[in] bl_levels   The number of boundary-layer levels
!! @param[in] ndf_wt  Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom for theta space
!! @param[in] map_wt  Cell dofmap for theta space
!! @param[in] ndf_w3  Number of degrees of freedom per cell for w3 space
!! @param[in] undf_w3  Number of unique degrees of freedom for w3 space
!! @param[in] map_w3  Cell dofmap for w3 space
subroutine leonard_term_th_code( nlayers,                               &
                                 field_inc,                             &
                                 field,                                 &
                                 map_wth_stencil_size, map_wth_stencil, &
                                 velocity_w2v,                          &
                                 map_wt_stencil_size, map_wt_stencil,   &
                                 dtrdz_tq_bl,                           &
                                 kl,                                    &
                                 rho,                                   &
                                 height_w3,                             &
                                 planet_radius, bl_levels,              &
                                 ndf_wt, undf_wt, map_wt,               &
                                 ndf_w3, undf_w3, map_w3                &
                                )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, bl_levels
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: map_wth_stencil_size
  integer(kind=i_def), dimension(ndf_wt,map_wth_stencil_size), intent(in)  :: map_wth_stencil
  integer(kind=i_def), intent(in) :: map_wt_stencil_size
  integer(kind=i_def), dimension(ndf_wt,map_wt_stencil_size), intent(in)  :: map_wt_stencil
  integer(kind=i_def), dimension(ndf_wt),  intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_w3),  intent(in)  :: map_w3

  real(kind=r_def), dimension(undf_wt),  intent(inout) :: field_inc
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: field
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: velocity_w2v
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: dtrdz_tq_bl
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: kl
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: rho
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: height_w3
  real(kind=r_def),                      intent(in)    :: planet_radius

  ! Internal variables
  integer(kind=i_def) :: k, kp

  ! Leonard term vertical flux on w3 levels
  real(kind=r_def), dimension(0:bl_levels) :: flux
  ! density * r^2 on w3 levels
  real(kind=r_def), dimension(1:bl_levels) :: rho_rsq

  ! If the full stencil isn't available, we must be at the domain edge.
  ! Simply set the increment to 0 for now, and exit the routine.
  if (map_wth_stencil_size < 5_i_def) then
    do k = 0, nlayers
      field_inc(map_wt(1) + k) = 0.0_r_def
    end do
    return
  end if

  ! Calculate rho * r^2
  do k = 1, bl_levels
     rho_rsq(k) = rho(map_w3(1) +k) *                             &
                  ( height_w3(map_w3(1) + k) + planet_radius ) *  &
                  ( height_w3(map_w3(1) + k) + planet_radius )
  end do

  ! Calculate vertical flux on w3 levels:
  ! Leonard term sub-grid vertical flux is proportional to the
  ! product of the grid-scale horizontal differences in field and w.

  ! Set flux and increment to zero at k=0
  k = 0
  flux(k) = 0.0_r_def
  field_inc(map_wt(1) + k) = 0.0_r_def

  do k = 1, bl_levels-1
    kp = k + 1
    flux(k) = ( kl(map_w3(1) + k) / 12.0_r_def )        &
            ! 4 terms contribute to each of x and y direction, so scale by 1/4
            * ( 1.0_r_def / 4.0_r_def ) * (                       &
            ! Terms from gradient in x-direction...
            ! from wth-point below:
              ( field(map_wth_stencil(1,1) + k) -                 &
                field(map_wth_stencil(1,2) + k) )                 &
            * ( velocity_w2v(map_wt_stencil(1,1) + k) -           &
                velocity_w2v(map_wt_stencil(1,2) + k) )           &
            + ( field(map_wth_stencil(1,4) + k) -                 &
                field(map_wth_stencil(1,1) + k) )                 &
            * ( velocity_w2v(map_wt_stencil(1,4) + k) -           &
                velocity_w2v(map_wt_stencil(1,1) + k) )           &
            ! from wth-point above:
            + ( field(map_wth_stencil(1,1) + kp) -                &
                field(map_wth_stencil(1,2) + kp) )                &
            * ( velocity_w2v(map_wt_stencil(1,1) + kp) -          &
                velocity_w2v(map_wt_stencil(1,2) + kp) )          &
            + ( field(map_wth_stencil(1,4) + kp) -                &
                field(map_wth_stencil(1,1) + kp) )                &
            * ( velocity_w2v(map_wt_stencil(1,4) + kp) -          &
                velocity_w2v(map_wt_stencil(1,1) + kp) )          &
            ! Terms from gradient in y-direction...
            ! from wth-point below:
            + ( field(map_wth_stencil(1,1) + k) -                 &
                field(map_wth_stencil(1,3) + k) )                 &
            * ( velocity_w2v(map_wt_stencil(1,1) + k) -           &
                velocity_w2v(map_wt_stencil(1,3) + k) )           &
            + ( field(map_wth_stencil(1,5) + k) -                 &
                field(map_wth_stencil(1,1) + k) )                 &
            * ( velocity_w2v(map_wt_stencil(1,5) + k) -           &
                velocity_w2v(map_wt_stencil(1,1) + k) )           &
            ! from wth-point above:
            + ( field(map_wth_stencil(1,1) + kp) -                &
                field(map_wth_stencil(1,3) + kp) )                &
            * ( velocity_w2v(map_wt_stencil(1,1) + kp) -          &
                velocity_w2v(map_wt_stencil(1,3) + kp) )          &
            + ( field(map_wth_stencil(1,5) + kp) -                &
                field(map_wth_stencil(1,1) + kp) )                &
            * ( velocity_w2v(map_wt_stencil(1,5) + kp) -          &
                velocity_w2v(map_wt_stencil(1,1) + kp) ) )        &
            * rho_rsq(k)
            ! Flux now includes density * r^2 factor
  end do

  ! Set flux to zero at k=nlayers-1
  k = bl_levels
  flux(k) = 0.0_r_def

! Difference the flux in the vertical to get increment on wth
  do k = 1, bl_levels
    field_inc(map_wt(1) + k) = -(flux(k) - flux(k-1))                     &
                               * dtrdz_tq_bl(map_wt(1) + k)
  end do

  ! Set increment to zero above bl_levels
  do k = bl_levels+1, nlayers
    field_inc(map_wt(1) + k) = 0.0_r_def
  end do

end subroutine leonard_term_th_code

end module leonard_term_th_kernel_mod
