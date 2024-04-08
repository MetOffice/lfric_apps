!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel which computes vertical cell edges value based on the monotone
!!        Koren scheme.
!> @details The kernel computes advective increment at Wtheta-points, using
!!          edge-values at w3-points computed with the Koren scheme.
!!  Ref.   Koren, B. (1993), "A robust upwind discretisation method for advection,
!!         diffusion and source terms", in Vreugdenhil, C.B.; Koren, B. (eds.),
!!         Numerical Methods for Advection/Diffusion Problems, Braunschweig:
!!         Vieweg, p. 117, ISBN 3-528-07645-3.
module polyv_wtheta_koren_kernel_mod

use argument_mod,         only : arg_type, GH_FIELD,     &
                                 GH_SCALAR, GH_REAL,     &
                                 GH_INTEGER, GH_LOGICAL, &
                                 GH_READWRITE,           &
                                 GH_READ, CELL_COLUMN
use constants_mod,        only : r_tran, r_def, i_def, l_def, tiny_eps, EPS_R_TRAN
use fs_continuity_mod,    only : W2v, Wtheta
use kernel_mod,           only : kernel_type
use koren_support_mod,    only : interpolate_to_regular_grid

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: polyv_wtheta_koren_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                 &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),         & ! advective update
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2v),            & ! wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),         & ! tracer
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),         & ! theta_height
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                      & ! reversible
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                       & ! logspace
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: polyv_wtheta_koren_code
end type
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: polyv_wtheta_koren_code

contains
!> @brief Computes the vertical fluxes for a tracer density.
!> @param[in]     nlayers      Number of layers
!> @param[in,out] advective    Advective update to increment
!> @param[in]     wind         Wind field
!> @param[in]     tracer       Tracer field to advect
!> @param[in]     theta_height Vertical height of Wtheta DOFs
!> @param[in]     logspace     Perform interpolation in log space
!> @param[in]     ndf_wt       Number of degrees of freedom per cell
!> @param[in]     undf_wt      Number of unique degrees of freedom for the tracer field
!> @param[in]     map_wt       Cell dofmaps for the tracer space
!> @param[in]     ndf_w2v      Number of degrees of freedom per cell
!> @param[in]     undf_w2v     Number of unique degrees of freedom for the flux &
!!                             wind fields
!> @param[in]     map_w2v      Dofmap for the cell at the base of the column
subroutine polyv_wtheta_koren_code( nlayers,              &
                                    advective,            &
                                    wind,                 &
                                    tracer,               &
                                    theta_height,         &
                                    reversible,           &
                                    logspace,             &
                                    ndf_wt,               &
                                    undf_wt,              &
                                    map_wt,               &
                                    ndf_w2v,              &
                                    undf_w2v,             &
                                    map_w2v               )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt
  integer(kind=i_def), intent(in) :: undf_wt
  integer(kind=i_def), intent(in) :: ndf_w2v
  integer(kind=i_def), intent(in) :: undf_w2v

  integer(kind=i_def), dimension(ndf_w2v), intent(in) :: map_w2v
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt

  real(kind=r_tran), dimension(undf_wt),  intent(inout) :: advective
  real(kind=r_tran), dimension(undf_w2v), intent(in)    :: wind
  real(kind=r_tran), dimension(undf_wt),  intent(in)    :: tracer
  real(kind=r_def),  dimension(undf_wt),  intent(in)    :: theta_height

  logical(kind=l_def), intent(in)    :: reversible
  logical(kind=l_def), intent(in)    :: logspace

  ! Internal variables
  integer(kind=i_def) :: k, k1, k2, k3

  real(kind=r_tran) :: x, y, r, r1, r2, phi
  real(kind=r_tran) :: t1, t2, t3
  real(kind=r_tran), dimension(nlayers+1)   :: wind_1d
  real(kind=r_tran), dimension(nlayers+1)   :: tracer_edge_t, tracer_edge_b
  real(kind=r_tran), dimension(nlayers)     :: dtracerdz
  real(kind=r_tran), dimension(0:nlayers+2) :: tracer_1d, dz
  real(kind=r_tran), dimension(nlayers+1)   :: zl
  real(kind=r_tran), dimension(nlayers)     :: zm

  !Extract vertical 1d-arrays from global data
  do k=0,nlayers
    wind_1d(k+1)   = wind(map_w2v(1)+k)
    tracer_1d(k+1) = tracer(map_wt(1)+k)
    zl(k+1)        = real(theta_height(map_wt(1)+k),r_tran)
  end do
  do k = 1,nlayers
    zm(k) = 0.5_r_tran*(zl(k) + zl(k+1))
  end do
  do k = 2,nlayers
    dz(k) = zm(k) - zm(k-1)
  end do
  dz(0:1) = 2.0_r_tran*(zm(1)-zl(1))
  dz(nlayers+1:nlayers+2) = 2.0_r_tran*(zl(nlayers+1)-zm(nlayers))

  ! Add 2 extra points for tracer_1d array outside the boundaries
  ! to avoid treating the boundaries differently from the rest
  ! of the column.
  !
  ! These extra points are constructed by preserving the slope
  ! near the boundaries (the slope is one-sided slope).

  t1 = tracer_1d(2) - tracer_1d(1)
  t2 = tracer_1d(nlayers+1) - tracer_1d(nlayers)
  tracer_1d(0) = tracer_1d(1) - t1
  tracer_1d(nlayers+2) = tracer_1d(nlayers+1) + t2

  ! Apply log to tracer_1d if required
  ! If using the logspace option, the tracer is forced to be positive

  if (logspace) then
    do k=0,nlayers+2
      tracer_1d(k) = log(max(EPS_R_TRAN,abs(tracer_1d(k))))
    end do
  end if

  ! Compute tracers at edges of cells centred around wtheta-points
  ! For Reversible case use linear reconstruction
  if (reversible) then
    do k = 1, nlayers
      ! Top edge tracer(k+1/2) using the 2 cells [k, k+1]
      t1 = dz(k) + dz(k+1)
      t2 = dz(k)/t1
      t3 = 1.0_r_tran - t2
      tracer_edge_t(k) = t3*tracer_1d(k) + t2*tracer_1d(k+1)
    end do
    if (logspace) then
      do k=1,nlayers
        tracer_edge_t(k) = exp(tracer_edge_t(k))
      end do
    end if
    do k = 2, nlayers
      dtracerdz(k) = tracer_edge_t(k) - tracer_edge_t(k-1)
    end do
  else
    do k = 1, nlayers + 1
      ! Top edge tracer(k+1/2) using the 3 cells [k-1, k, k+1] (w>0)
      k1 = k - 1
      k2 = k
      k3 = k + 1

      ! Compute 3 upwind regular grid values (t1,t2,t3) from
      ! the non-regular ones, tracer_1d(k1:k3).
      ! The koren scheme will be applied to (t1,t2,t3) instead
      ! of the non-regular ones, tracer_1d(k1:k3).

      call interpolate_to_regular_grid(tracer_1d(k1),        &
                                       tracer_1d(k2),        &
                                       tracer_1d(k3),        &
                                       dz(k1),dz(k2),dz(k3), &
                                       t1,t2,t3)
      x = t2 - t1
      y = t3 - t2
      r = (y + tiny_eps)/(x + tiny_eps)
      r1 = 2.0_r_tran*r
      r2 = (1.0_r_tran + r1)/3.0_r_tran
      phi = max(0.0_r_tran, min(r1,r2,2.0_r_tran))
      tracer_edge_t(k) = t2 + 0.5_r_tran*phi*x

      ! Bottom edge tracer(k-1/2) using the 3 cells [k+1, k, k-1] (w<0)
      ! Note that we can apply the Koren scheme for the same values (t1,t2,t3)
      ! used for top edge but reversed (t3,t2,t1) [t3 --> t1 and t1 --> t3]

      x = t2 - t3
      y = t1 - t2
      r = (y + tiny_eps)/(x + tiny_eps)
      r1 = 2.0_r_tran*r
      r2 = (1.0_r_tran + r1)/3.0_r_tran
      phi = max(0.0_r_tran, min(r1,r2,2.0_r_tran))
      tracer_edge_b(k) = t2 + 0.5_r_tran*phi*x
    end do

    if (logspace) then
      do k=1,nlayers + 1
         tracer_edge_t(k) = exp(tracer_edge_t(k))
         tracer_edge_b(k) = exp(tracer_edge_b(k))
      end do
    end if

    do k = 2, nlayers
      if ( wind_1d(k) > 0.0_r_tran ) then
        dtracerdz(k) = tracer_edge_t(k) - tracer_edge_t(k-1)
      else
        dtracerdz(k) = tracer_edge_b(k+1) - tracer_edge_b(k)
      end if
    end do
  end if

  do k = 1, nlayers - 1
    advective(map_wt(1) + k ) = advective(map_wt(1) + k )   &
                                + wind(map_w2v(1) + k )     &
                                * dtracerdz(k+1)
  end do

end subroutine polyv_wtheta_koren_code

end module polyv_wtheta_koren_kernel_mod
