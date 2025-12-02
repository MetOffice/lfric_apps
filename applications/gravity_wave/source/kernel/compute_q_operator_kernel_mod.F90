!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

module compute_q_operator_kernel_mod

use argument_mod,      only: arg_type, func_type,   &
                             GH_OPERATOR, GH_REAL,  &
                             GH_WRITE, ANY_SPACE_1, &
                             GH_BASIS, CELL_COLUMN, &
                             GH_QUADRATURE_XYoZ
use constants_mod,     only: r_def, i_def, r_single, r_double
use fs_continuity_mod, only: W2
use kernel_mod,        only: kernel_type

implicit none

private

! Precomputed arrays
real(kind=r_def), allocatable, private :: delta_z(:)

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_q_operator_type
  private
  type(arg_type) :: meta_args(1) = (/                            &
       arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, ANY_SPACE_1) &
       /)
  type(func_type) :: meta_funcs(2) = (/                          &
       func_type(W2,          GH_BASIS),                         &
       func_type(ANY_SPACE_1, GH_BASIS)                          &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: compute_q_operator_init
public :: compute_q_operator_code
public :: compute_q_operator_final

  ! Generic interface for real32 and real64 types
  interface compute_q_operator_code
    module procedure  &
      compute_q_operator_code_r_single, &
      compute_q_operator_code_r_double
  end interface

contains

!> @brief Computes the q operator which is the projection of the vertical
!!        bouyancy term
!! @param[in] cell Cell id
!! @param[in] nlayers Number of layers
!! @param[in] ncell_3d ncell*nlayers
!! @param[in] ndf_w2 Number of degrees of freedom per cell
!! @param[in] basis_w2 Basis functions evaluated at quadrature points
!! @param[in] ndf_wt Number of degrees of freedom per cell
!! @param[in] basis_wt Basis functions evaluated at quadrature points
!! @param[in,out] q Local stencil of the q operator
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights

! R_DOUBLE PRECISION
! ==================
subroutine compute_q_operator_code_r_double(cell, nlayers, ncell_3d,     &
                                            q,                           &
                                            ndf_w2, basis_w2,            &
                                            ndf_wt, basis_wt,            &
                                            nqp_h, nqp_v, wqp_h, wqp_v )

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nqp_h, nqp_v
  integer(kind=i_def),                    intent(in) :: nlayers
  integer(kind=i_def),                    intent(in) :: ncell_3d
  integer(kind=i_def),                    intent(in) :: ndf_wt, ndf_w2

  real(kind=r_def), intent(in) :: basis_wt(1,ndf_wt,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: basis_w2(3,ndf_w2,nqp_h,nqp_v)

  real(kind=r_double), dimension(ncell_3d,ndf_w2,ndf_wt), intent(inout) :: q
  real(kind=r_def),    dimension(nqp_h),                  intent(in)    :: wqp_h
  real(kind=r_def),    dimension(nqp_v),                  intent(in)    :: wqp_v

  ! Internal variables
  integer(kind=i_def)            :: df2, dft, k, ik
  integer(kind=i_def)            :: qp1, qp2
  real(kind=r_double)            :: integrand
  real(kind=r_def), dimension(3) :: z_hat(3)


  do k = 0, nlayers - 1
    z_hat = (/ 0.0_r_def, 0.0_r_def, delta_z(k+1) /)
    ik = k + 1 + (cell-1)*nlayers
    do dft = 1, ndf_wt
      do df2 = 1, ndf_w2
        q(ik,df2,dft) = 0.0_r_double
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h
            integrand = real(wqp_h(qp1)*wqp_v(qp2)                &
                      *dot_product(basis_w2(:,df2,qp1,qp2),z_hat) &
                      *basis_wt(1,dft,qp1,qp2), r_double)
            q(ik,df2,dft) = q(ik,df2,dft) + integrand
          end do
        end do
      end do
    end do
  end do
end subroutine compute_q_operator_code_r_double

! R_SINGLE PRECISION
! ==================
subroutine compute_q_operator_code_r_single(cell, nlayers, ncell_3d,     &
                                            q,                           &
                                            ndf_w2, basis_w2,            &
                                            ndf_wt, basis_wt,            &
                                            nqp_h, nqp_v, wqp_h, wqp_v )

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nqp_h, nqp_v
  integer(kind=i_def),                    intent(in) :: nlayers
  integer(kind=i_def),                    intent(in) :: ncell_3d
  integer(kind=i_def),                    intent(in) :: ndf_wt, ndf_w2

  real(kind=r_def), intent(in) :: basis_wt(1,ndf_wt,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: basis_w2(3,ndf_w2,nqp_h,nqp_v)

  real(kind=r_single), dimension(ncell_3d,ndf_w2,ndf_wt), intent(inout) :: q
  real(kind=r_def),    dimension(nqp_h),                  intent(in)    :: wqp_h
  real(kind=r_def),    dimension(nqp_v),                  intent(in)    :: wqp_v

  ! Internal variables
  integer(kind=i_def)            :: df2, dft, k, ik
  integer(kind=i_def)            :: qp1, qp2
  real(kind=r_single)            :: integrand
  real(kind=r_def), dimension(3) :: z_hat(3)


  do k = 0, nlayers - 1
    z_hat = (/ 0.0_r_def, 0.0_r_def, delta_z(k+1) /)
    ik = k + 1 + (cell-1)*nlayers
    do dft = 1, ndf_wt
      do df2 = 1, ndf_w2
        q(ik,df2,dft) = 0.0_r_single
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h
            integrand = real(wqp_h(qp1)*wqp_v(qp2)                &
                      *dot_product(basis_w2(:,df2,qp1,qp2),z_hat) &
                      *basis_wt(1,dft,qp1,qp2), r_single)
            q(ik,df2,dft) = q(ik,df2,dft) + integrand
          end do
        end do
      end do
    end do
  end do
end subroutine compute_q_operator_code_r_single

!=============================================================================!
!> @brief Initialise the q computation kernel, copies the dz array into a local
!!        copy
!! @param[in] dz Layer thickness array
!! @param[in] nlayers The number of layers
subroutine compute_q_operator_init(dz, nlayers)

  implicit none

  integer(kind=i_def),                  intent(in) :: nlayers
  real(kind=r_def), dimension(nlayers), intent(in) :: dz

  allocate( delta_z(nlayers) )
  delta_z = dz
end subroutine compute_q_operator_init

!=============================================================================!
!> @brief Reclaims memory from private allocatable arrays
subroutine compute_q_operator_final()

  implicit none

  if (allocated(delta_z)) deallocate (delta_z)

end subroutine compute_q_operator_final

end module compute_q_operator_kernel_mod
