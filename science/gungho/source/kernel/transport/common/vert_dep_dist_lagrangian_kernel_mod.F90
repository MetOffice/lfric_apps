!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the departure distances for cell faces in the
!!        vertical direction, using a Lagrangian method.
!> @details This code calculates the distance which is swept through a cell
!!          in the z direction during one timestep. The arrival point is the
!!          cell face and the departure point is calculated. This kernel
!!          performs Lagrangian calculations: options are the midpoint rule or
!!          the trapezoidal rule.
!!          This kernel returns the distance between the arrival and departure
!!          point for each cell face in the vertical and this value is positive
!!          if the w wind (radial wind) is positive, i.e. increasing in height.
!!          The integer part of the distance is the number of complete cells and
!!          the fractional part is the fraction of the cell in which the
!!          departure point resides.

module vert_dep_dist_lagrangian_kernel_mod

use argument_mod,                only : arg_type,              &
                                        GH_FIELD, GH_REAL,     &
                                        GH_WRITE, GH_READ,     &
                                        GH_SCALAR, GH_INTEGER, &
                                        CELL_COLUMN
use fs_continuity_mod,           only : W2, W2v
use constants_mod,               only : r_tran, i_def
use kernel_mod,                  only : kernel_type
use departure_points_config_mod, only : vertical_limit_boundary, &
                                        vertical_limit_exponential
implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vert_dep_dist_lagrangian_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/               &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2v), & ! dep_dist
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2v), & ! cfl
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2v), & ! u_n
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2v), & ! u_np1
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),  & ! heights_w2
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),    & ! iterations
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),    & ! method
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),    & ! limit
       arg_type(GH_SCALAR, GH_REAL, GH_READ)        & ! dt
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vert_dep_dist_lagrangian_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vert_dep_dist_lagrangian_code

contains

!> @brief Kernel which computes the departure distances for cell faces in the
!!        vertical direction, using a Lagrangian method.
!> @param[in]     nlayers             The number of layers
!> @param[in,out] dep_dist_z          The departure distances in the vertical
!> @param[in,out] cfl                 The vertical CFL calculated via the
!!                                    vertical departure points
!> @param[in]     u_n                 The W2v wind field at time level n
!> @param[in]     u_np1               The W2v wind field at time level n+1
!> @param[in]     height_w2           Physical height of W2 dofs
!> @param[in]     n_dep_pt_iterations The number of departure point iterations
!> @param[in]     vertical_method     Enumerator for the vertical method to be
!!                                    used for computing the departure points
!> @param[in]     dt                  The model timestep length
!> @param[in]     vertical_limit      Enumerator for the method to be used to
!!                                    limit the vertical departure points
!> @param[in]     ndf_w2v             The number of degrees of freedom per W2v cell
!> @param[in]     undf_w2v            The number of unique degrees of freedom for W2v
!> @param[in]     map_w2v             The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w2              The number of degrees of freedom per cell
!> @param[in]     undf_w2             The number of unique degrees of freedom
!> @param[in]     map_w2              The dofmap for the cell at the base of the column
subroutine vert_dep_dist_lagrangian_code( nlayers,             &
                                          dep_dist_z,          &
                                          cfl,                 &
                                          u_n,                 &
                                          u_np1,               &
                                          height_w2,           &
                                          n_dep_pt_iterations, &
                                          vertical_method,     &
                                          vertical_limit,      &
                                          dt,                  &
                                          ndf_w2v,             &
                                          undf_w2v,            &
                                          map_w2v,             &
                                          ndf_w2,              &
                                          undf_w2,             &
                                          map_w2 )

  use departure_points_support_mod, only : calc_vertical_dep_cfl, &
                                           vertical_increasing_check

  implicit none

  ! Arguments - DOFs
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2
  integer(kind=i_def), intent(in) :: undf_w2
  integer(kind=i_def), intent(in) :: ndf_w2v
  integer(kind=i_def), intent(in) :: undf_w2v

  ! Arguments - Maps
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w2v), intent(in) :: map_w2v

  ! Arguments
  real(kind=r_tran),   dimension(undf_w2v), intent(inout) :: dep_dist_z
  real(kind=r_tran),   dimension(undf_w2v), intent(inout) :: cfl
  real(kind=r_tran),   dimension(undf_w2v), intent(in)    :: u_n
  real(kind=r_tran),   dimension(undf_w2v), intent(in)    :: u_np1
  real(kind=r_tran),   dimension(undf_w2),  intent(in)    :: height_w2
  integer(kind=i_def),                      intent(in)    :: n_dep_pt_iterations
  integer(kind=i_def),                      intent(in)    :: vertical_method
  integer(kind=i_def),                      intent(in)    :: vertical_limit
  real(kind=r_tran),                        intent(in)    :: dt

  ! Indices
  integer(kind=i_def) :: k

  ! Local fields
  integer(kind=i_def) :: nCellEdges
  real(kind=r_tran)   :: xArrival_comp
  real(kind=r_tran)   :: xArrival_phys
  real(kind=r_tran)   :: u_n_local(1:nlayers+1)
  real(kind=r_tran)   :: u_np1_local(1:nlayers+1)
  real(kind=r_tran)   :: u_dep(1:nlayers+1)
  real(kind=r_tran)   :: height_local(1:nlayers+1)
  real(kind=r_tran)   :: dz_local_m, dz_local_p
  real(kind=r_tran)   :: upwind_dz_n
  real(kind=r_tran)   :: upwind_dz_np1
  real(kind=r_tran)   :: dep_local(1:nlayers-1)
  real(kind=r_tran)   :: cfl_local

  ! Number of cell edgs
  nCellEdges = nlayers+1

  ! Compute the physical velocity as departure_wind*upwind_dz
  ! Apply vertical boundary conditions
  height_local(1) = height_w2(map_w2(5))
  u_dep(1)        = 0.0_r_tran
  u_n_local(1)    = 0.0_r_tran
  u_np1_local(1)  = 0.0_r_tran
  ! Compute local values
  do k=1,nlayers-1
    dz_local_m = height_w2(map_w2(6)+k-1) - height_w2(map_w2(5)+k-1)
    dz_local_p = height_w2(map_w2(6)+k)   - height_w2(map_w2(5)+k)

    upwind_dz_n   = ( 0.5_r_tran + sign( 0.5_r_tran, u_n(map_w2v(1)+k) ) ) * dz_local_m + &
                    ( 0.5_r_tran - sign( 0.5_r_tran, u_n(map_w2v(1)+k) ) ) * dz_local_p
    upwind_dz_np1 = ( 0.5_r_tran + sign( 0.5_r_tran, u_np1(map_w2v(1)+k) ) ) * dz_local_m + &
                    ( 0.5_r_tran - sign( 0.5_r_tran, u_np1(map_w2v(1)+k) ) ) * dz_local_p
    u_dep(k+1)       = 0.5_r_tran * ( u_n(map_w2v(1)+k) + u_np1(map_w2v(1)+k) )
    u_n_local(k+1)   = u_n(map_w2v(1)+k) * upwind_dz_n
    u_np1_local(k+1) = u_np1(map_w2v(1)+k) * upwind_dz_np1
    height_local(k+1) = height_w2(map_w2(5)+k)
  end do
  ! Apply vertical boundary conditions
  height_local(nCellEdges) = height_w2(map_w2(6)+nlayers-1)
  u_dep(nCellEdges)        = 0.0_r_tran
  u_n_local(nCellEdges)    = 0.0_r_tran
  u_np1_local(nCellEdges)  = 0.0_r_tran

  ! Loop over all layers except the bottom layer.
  ! This code is hard-wired to work with 6 W2 dofs per cell where dof=5 is the
  ! vertical dof at the bottom of the cell, and 2 W2v dofs per cell
  do k=1,nlayers-1
    xArrival_comp = real(k,r_tran)
    xArrival_phys = height_w2(map_w2(5)+k)
    call calc_vertical_dep_cfl( xArrival_comp,        &
                                xArrival_phys,        &
                                nCellEdges,           &
                                u_n_local,            &
                                u_np1_local,          &
                                u_dep,                &
                                height_local,         &
                                dt,                   &
                                vertical_method,      &
                                n_dep_pt_iterations,  &
                                vertical_limit,       &
                                dep_local(k),         &
                                cfl_local )
    cfl( map_w2v(1) + k ) =  cfl_local
  end do

  if ( vertical_limit == vertical_limit_exponential ) then
    ! Ensure vertical departure points are monotonic
    call vertical_increasing_check(dep_local, nlayers)
  end if

  ! Compute departure points.
  dep_dist_z( map_w2v(1) ) =  0.0_r_tran
  do k=1,nlayers-1
    xArrival_comp = real(k,r_tran)
    dep_dist_z( map_w2v(1) + k ) =  xArrival_comp - dep_local(k)
  end do
  dep_dist_z( map_w2v(2)+nlayers-1 ) =  0.0_r_tran

end subroutine vert_dep_dist_lagrangian_code

end module vert_dep_dist_lagrangian_kernel_mod
