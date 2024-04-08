!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to compute the vertical quintic interpolation coefficients.
!> @Details Computes the quintic vertical interpolation coefficients
!!          used by the semi-Lagrangian transport scheme.

module compute_vertical_quintic_coef_kernel_mod

  use argument_mod,          only : arg_type,                     &
                                    GH_FIELD, GH_WRITE,           &
                                    GH_REAL, GH_READ,             &
                                    GH_INTEGER, CELL_COLUMN,      &
                                    ANY_DISCONTINUOUS_SPACE_1
  use fs_continuity_mod,     only : W2v, Wtheta
  use constants_mod,         only : r_tran, i_def, l_def, EPS_R_TRAN
  use kernel_mod,            only : kernel_type
  use sl_support_mod,        only : compute_quintic_coeffs

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed
  !>                                      by the PSy layer.
  type, public, extends(kernel_type) :: compute_vertical_quintic_coef_kernel_type
    private
    type(arg_type) :: meta_args(16) = (/                                       &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v),                       & ! dep_dist_z
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! theta_height
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_coef
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! quintic_indices
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! linear_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! linear_coef
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: compute_vertical_quintic_coef_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: compute_vertical_quintic_coef_code

  contains

  !-------------------------------------------------------------------------------
  !> @details This kernel calculates the departure point of w/theta-points using
  !!          only w (i.e., vertical motion only), then interpolate theta at the
  !!          departure point using 1d-Cubic-Lagrange interpolation.
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     dep_dist_z    The vertical departure distance used for SL advection
  !> @param[in]     theta_height  The height of theta-points
  !> @param[inout]  quintic_coef  The cubic interpolation coefficients (1-6)
  !> @param[inout]  quintic_indices
  !!                              The cubic interpolation indices (1-6)
  !> @param[inout]  linear_coef   The linear interpolation coefficients (1-2, used for monotonicity)
  !> @param[in]     ndf_w2v       The number of degrees of freedom per cell
  !!                              on W2v space
  !> @param[in]     undf_w2v      The number of unique degrees of freedom
  !!                              on W2v space
  !> @param[in]     map_w2v       The dofmap for the cell at the base of the column
  !!                              on W2v space
  !> @param[in]     ndf_wt        The number of degrees of freedom per cell
  !!                              on Wtheta space
  !> @param[in]     undf_wt       The number of unique degrees of freedom
  !!                              on Wtheta space
  !> @param[in]     map_wt        The dofmap for the cell at the base of the column
  !!                              on Wtheta space
  !> @param[in]     ndf_wq        The number of degrees of freedom per cell
  !!                              for the coefficients space
  !> @param[in]     undf_wq       The number of unique degrees of freedom
  !!                              for the coefficients space
  !> @param[in]     map_wq        The dofmap for the cell at the base of the column
  !!                              for the coefficients space
  !-------------------------------------------------------------------------------

  subroutine compute_vertical_quintic_coef_code( nlayers,                              &
                                                 dep_dist_z,                           &
                                                 theta_height,                         &
                                                 quintic_coef_1, quintic_coef_2,       &
                                                 quintic_coef_3, quintic_coef_4,       &
                                                 quintic_coef_5, quintic_coef_6,       &
                                                 quintic_indices_1, quintic_indices_2, &
                                                 quintic_indices_3, quintic_indices_4, &
                                                 quintic_indices_5, quintic_indices_6, &
                                                 linear_coef_1, linear_coef_2,         &
                                                 ndf_w2v, undf_w2v, map_w2v,           &
                                                 ndf_wt, undf_wt, map_wt,              &
                                                 ndf_wq, undf_wq, map_wq )

    implicit none

    ! Arguments
    integer(kind=i_def),                      intent(in) :: nlayers
    integer(kind=i_def),                      intent(in) :: ndf_w2v
    integer(kind=i_def),                      intent(in) :: undf_w2v
    integer(kind=i_def), dimension(ndf_w2v),  intent(in) :: map_w2v
    integer(kind=i_def),                      intent(in) :: ndf_wt
    integer(kind=i_def),                      intent(in) :: undf_wt
    integer(kind=i_def), dimension(ndf_wt),   intent(in) :: map_wt
    integer(kind=i_def),                      intent(in) :: ndf_wq
    integer(kind=i_def),                      intent(in) :: undf_wq
    integer(kind=i_def), dimension(ndf_wq),   intent(in) :: map_wq
    real(kind=r_tran),   dimension(undf_w2v), intent(in) :: dep_dist_z
    real(kind=r_tran),   dimension(undf_wt),  intent(in) :: theta_height
    real(kind=r_tran),   dimension(undf_wq),  intent(inout) :: quintic_coef_1
    real(kind=r_tran),   dimension(undf_wq),  intent(inout) :: quintic_coef_2
    real(kind=r_tran),   dimension(undf_wq),  intent(inout) :: quintic_coef_3
    real(kind=r_tran),   dimension(undf_wq),  intent(inout) :: quintic_coef_4
    real(kind=r_tran),   dimension(undf_wq),  intent(inout) :: quintic_coef_5
    real(kind=r_tran),   dimension(undf_wq),  intent(inout) :: quintic_coef_6
    integer(kind=i_def), dimension(undf_wq),  intent(inout) :: quintic_indices_1
    integer(kind=i_def), dimension(undf_wq),  intent(inout) :: quintic_indices_2
    integer(kind=i_def), dimension(undf_wq),  intent(inout) :: quintic_indices_3
    integer(kind=i_def), dimension(undf_wq),  intent(inout) :: quintic_indices_4
    integer(kind=i_def), dimension(undf_wq),  intent(inout) :: quintic_indices_5
    integer(kind=i_def), dimension(undf_wq),  intent(inout) :: quintic_indices_6
    real(kind=r_tran),   dimension(undf_wq),  intent(inout) :: linear_coef_1
    real(kind=r_tran),   dimension(undf_wq),  intent(inout) :: linear_coef_2

    ! Local arrays and indices
    real(kind=r_tran),  dimension(nlayers+1)    :: zm, zmd
    real(kind=r_tran),  dimension(nlayers+2)    :: dist, zl, zld
    integer(kind=i_def)                         :: k, nz, nl, nzl, km1, km2, si, nlp
    real(kind=r_tran)                           :: d, r, sr
    real(kind=r_tran),   dimension(6,nlayers+1) :: qc
    integer(kind=i_def), dimension(6,nlayers+1) :: sq
    real(kind=r_tran),   dimension(2,nlayers+1) :: cl
    real(kind=r_tran),   dimension(nlayers+1)   :: dz

    ! nl = nlayers-1  for W3
    !    = nlayers    for Wtheta
    nl  = nlayers - 1 + (ndf_wq - 1)
    nlp = nl+1
    nz  = nlayers
    nzl = nlayers + 1

    ! Extract and fill local column from global data
    ! Map departure points into 1d-array dist
    ! Map theta-height to 1d-zl (zl is cell-edges in the vertical)
    do k = 0, nlayers
      dist(k+1) = dep_dist_z(map_w2v(1)+k)
      zl(k+1)   = theta_height(map_wt(1)+k)
    end do

    ! Recover the physical departure points of cell edges zld
    do k = 1, nzl
       d     = abs(dist(k))
       sr    = sign(1.0_r_tran,dist(k))
       si    = int(sr,i_def)
       km1   = int( d,i_def)
       r     = d - real(km1,r_tran)
       km1   = k - km1*si
       km2   = km1 - si
       km1   = max(1_i_def, min(km1,nzl))
       km2   = max(1_i_def, min(km2,nzl))
       zld(k) = zl(km1) - sr*r*abs(zl(km2)-zl(km1))
       zld(k) = min(zl(nzl),max(zl(1),zld(k)))
    end do

    if ( ndf_wq == 1 ) then
      ! W3 field:
      ! Create a local 1d SL problem for cell centres
      ! Compute zm from zl and zmd from zld of cells where,
      ! (zm,zmd) are averaged from (zl,zld)
      !
      do k = 1, nz
        zm(k) = 0.5_r_tran*(zl(k) +  zl(k+1))
      end do
      do k = 1, nz
        zmd(k) = 0.5_r_tran*( zld(k) + zld(k+1) )
        zmd(k) = min(zm(nz),max(zm(1),zmd(k)))
      end do
      ! Define the spacing dz between zm-points
      do k = 1, nz - 1
        dz(k) = zm(k+1) - zm(k)
      end do
      dz(nz) = dz(nz - 1)
      ! Compute the quintic-interpolation coefficients
      call compute_quintic_coeffs(zmd,zm,dz,sq,qc,cl,nz,nz)
    else
      ! Wtheta field:
      ! Define the spacing dz between zl-points
      do k = 1, nzl - 1
        dz(k) = zl(k+1) - zl(k)
      end do
      dz(nzl) = dz(nzl-1)
      ! Compute the quintic-interpolation coefficients
      call compute_quintic_coeffs(zld,zl,dz,sq,qc,cl,nzl,nzl)
    end if

    ! Coeffs are W3/Wtheta fields
    do k = 0, nl
      quintic_coef_1(map_wq(1)+k) = qc(1,k+1)
      quintic_coef_2(map_wq(1)+k) = qc(2,k+1)
      quintic_coef_3(map_wq(1)+k) = qc(3,k+1)
      quintic_coef_4(map_wq(1)+k) = qc(4,k+1)
      quintic_coef_5(map_wq(1)+k) = qc(5,k+1)
      quintic_coef_6(map_wq(1)+k) = qc(6,k+1)
      quintic_indices_1(map_wq(1)+k) = sq(1,k+1)
      quintic_indices_2(map_wq(1)+k) = sq(2,k+1)
      quintic_indices_3(map_wq(1)+k) = sq(3,k+1)
      quintic_indices_4(map_wq(1)+k) = sq(4,k+1)
      quintic_indices_5(map_wq(1)+k) = sq(5,k+1)
      quintic_indices_6(map_wq(1)+k) = sq(6,k+1)
      linear_coef_1(map_wq(1)+k) = cl(1,k+1)
      linear_coef_2(map_wq(1)+k) = cl(2,k+1)
    end do

  end subroutine compute_vertical_quintic_coef_code

end module compute_vertical_quintic_coef_kernel_mod
