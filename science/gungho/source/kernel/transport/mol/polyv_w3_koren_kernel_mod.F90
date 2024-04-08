!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel which computes vertical cell egdes value based on the monotone
!!        Koren scheme.
!> @details The kernel computes vertical cell edges at W2v-points from
!!          W3-field using the Koren scheme.
!!  Ref.   Koren, B. (1993), "A robust upwind discretisation method for advection,
!!         diffusion and source terms", in Vreugdenhil, C.B.; Koren, B. (eds.),
!!         Numerical Methods for Advection/Diffusion Problems, Braunschweig:
!!         Vieweg, p. 117, ISBN 3-528-07645-3.
module polyv_w3_koren_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_LOGICAL, GH_WRITE,        &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              GH_READ, CELL_COLUMN
use fs_continuity_mod, only : W3, Wtheta
use constants_mod,     only : i_def, l_def, tiny_eps, EPS_R_TRAN, r_tran, r_def
use kernel_mod,        only : kernel_type
use koren_support_mod, only:  interpolate_to_regular_grid

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: polyv_w3_koren_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                         &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  & ! reconstruction
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                         & ! tracer
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                     & ! theta_height
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                              & ! reversible
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                               & ! logspace
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: polyv_w3_koren_code
end type
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: polyv_w3_koren_code
contains
!> @brief Computes the vertical reconstructions for a tracer.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] reconstruction Mass reconstruction field to compute
!> @param[in]     tracer         Tracer field
!> @param[in]     theta_height   Vertical height of Wtheta DOFs
!> @param[in]     reversible     Use the reversible scheme
!> @param[in]     logspace       Perform interpolation in log space;
!> @param[in]     ndf_md         Number of degrees of freedom per cell
!> @param[in]     undf_md        Number of unique degrees of freedom for the
!!                               reconstrution
!> @param[in]     map_md         Dofmap for the cell at the base of the column
!> @param[in]     ndf_w3         Number of degrees of freedom per cell
!> @param[in]     undf_w3        Number of unique degrees of freedom for the tracer field
!> @param[in]     map_w3         Cell dofmaps for the tracer space
!> @param[in]     ndf_wtheta     Number of degrees of cell-height (theta_height)
!> @param[in]     undf_wtheta    Number of unique degrees of freedom for theta_height
!> @param[in]     map_wtheta     Cell dofmaps for the theta_height
subroutine polyv_w3_koren_code( nlayers,                           &
                                reconstruction,                    &
                                tracer,                            &
                                theta_height,                      &
                                reversible,                        &
                                logspace,                          &
                                ndf_md,                            &
                                undf_md,                           &
                                map_md,                            &
                                ndf_w3,                            &
                                undf_w3,                           &
                                map_w3,                            &
                                ndf_wtheta,                        &
                                undf_wtheta,                       &
                                map_wtheta                         )
  implicit none
  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_md
  integer(kind=i_def), intent(in) :: undf_md
  integer(kind=i_def), intent(in) :: ndf_w3
  integer(kind=i_def), intent(in) :: undf_w3
  integer(kind=i_def), intent(in) :: ndf_wtheta
  integer(kind=i_def), intent(in) :: undf_wtheta

  integer(kind=i_def), dimension(ndf_md),     intent(in) :: map_md
  integer(kind=i_def), dimension(ndf_w3),     intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta

  real(kind=r_tran), dimension(undf_md),     intent(inout) :: reconstruction
  real(kind=r_tran), dimension(undf_w3),     intent(in)    :: tracer
  real(kind=r_def),  dimension(undf_wtheta), intent(in)    :: theta_height

  logical(kind=l_def), intent(in) :: reversible
  logical(kind=l_def), intent(in) :: logspace

  ! Local variables
  integer(kind=i_def), parameter                  :: ext = 2
  integer(kind=i_def)                             :: k, k1, k2, k3
  real(kind=r_tran)                               :: x, y, r, r1, r2, phi
  real(kind=r_tran)                               :: t1, t2, t3
  real(kind=r_tran), dimension(nlayers+1)         :: tracer_edge_t,tracer_edge_b
  real(kind=r_tran), dimension(1-ext:nlayers+ext) :: tracer_1d
  real(kind=r_tran), dimension(nlayers+1)         :: zl
  real(kind=r_tran), dimension(1-ext:nlayers+ext) :: dz

  ! Extract the global data into 1d-array
  do k = 0,nlayers
    zl(k+1) = real(theta_height(map_wtheta(1)+k),r_tran)
  end do
  do k = 0,nlayers - 1
    tracer_1d(k+1) = tracer(map_w3(1)+k)
  end do
  do k = 1,nlayers
    dz(k) = zl(k+1) - zl(k)
  end do

  ! Extend the tracer_1d array to go beyond the boundaries
  ! to avoid treating the boundaries differently from the rest
  ! of the column.
  !
  ! Ideally these extra points can be constructed by
  ! preserving the slope near the boundaries, but this
  ! approach can cause the data to go outside the initial bounds.
  ! Therefore we use zero slope assumption.

  tracer_1d(1-ext:0) = tracer_1d(1)
  tracer_1d(nlayers+1:nlayers+ext) = tracer_1d(nlayers)
  dz(1-ext:0) = dz(1)
  dz(nlayers+1:nlayers+ext) = dz(nlayers)

  ! Apply log to tracer_1d if required
  ! If using the logspace option, the tracer is forced to be positive

  if (logspace) then
      do k = 1-ext,nlayers+ext
         tracer_1d(k) = log(max(EPS_R_TRAN,abs(tracer_1d(k))))
      end do
  end if

  ! Compute tracers at edges of cells centred around w3-points
  ! For Reversible case use linear reconstruction
  if (reversible) then
    do k = 1, nlayers
      ! Top edge tracer(k+1/2) using the 2 cells [k, k+1]
      t1 = dz(k) + dz(k+1)
      t2 = dz(k)/t1
      t3 = 1.0_r_tran - t2
      tracer_edge_t(k) = t3*tracer_1d(k) + t2*tracer_1d(k+1)

      ! Bottom edge tracer(k-1/2) using the 2 cells [k-1, k]
      t1 = dz(k-1) + dz(k)
      t2 = dz(k)/t1
      t3 = 1.0_r_tran - t2
      tracer_edge_b(k) = t2*tracer_1d(k-1) + t3*tracer_1d(k)
    end do
  else
    do k = 1, nlayers
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
  end if

  if (logspace) then
    do k=1,nlayers
      tracer_edge_t(k) = exp(tracer_edge_t(k))
      tracer_edge_b(k) = exp(tracer_edge_b(k))
    end do
  end if

  ! Remap tracer_edge_b to reconstruction(map_md(1)+k[k=0:nlayers-1]) and
  !       tracer_edge_t to reconstruction(map_md(1)+nlyaers+k[k=0:nlayers-1])
  do k = 0, nlayers - 1
    reconstruction(map_md(1) + k )           =  tracer_edge_b(k+1)
    reconstruction(map_md(1) + nlayers + k ) =  tracer_edge_t(k+1)
  end do
end subroutine polyv_w3_koren_code

end module polyv_w3_koren_kernel_mod
