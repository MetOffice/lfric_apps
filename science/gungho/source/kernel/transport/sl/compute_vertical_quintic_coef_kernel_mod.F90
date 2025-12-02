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
    type(arg_type) :: meta_args(14) = (/                                       &
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
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! quintic_indices
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
  !!          departure point using 1d-Quintic-Lagrange interpolation.
  !> @param[in]     nlayers          The number of layers
  !> @param[in]     dep_dist_z       The vertical departure distance
  !> @param[in]     theta_height     The height of theta-points
  !> @param[in,out] quintic_coef     The cubic interpolation coefficients (1-6)
  !> @param[in,out] quintic_indices  The cubic interpolation indices (1-6)
  !> @param[in]     ndf_w2v          Num DoFs per cell for W2v
  !> @param[in]     undf_w2v         Num DoFs for this rank for W2v
  !> @param[in]     map_w2v          DoF map for W2v
  !> @param[in]     ndf_wt           Num DoFs per cell for Wtheta
  !> @param[in]     undf_wt          Num DoFs for this rank for Wtheta
  !> @param[in]     map_wt           DoF map for Wtheta
  !> @param[in]     ndf_wc           Num DoFs per cell for the coefficients
  !> @param[in]     undf_wc          Num DoFs for this rank for the coefficients
  !> @param[in]     map_wc           DoF map for the coefficients
  !-------------------------------------------------------------------------------
  subroutine compute_vertical_quintic_coef_code( nlayers,                      &
                                                 dep_dist_z,                   &
                                                 theta_height,                 &
                                                 quintic_coef_1,               &
                                                 quintic_coef_2,               &
                                                 quintic_coef_3,               &
                                                 quintic_coef_4,               &
                                                 quintic_coef_5,               &
                                                 quintic_coef_6,               &
                                                 quintic_indices_1,            &
                                                 quintic_indices_2,            &
                                                 quintic_indices_3,            &
                                                 quintic_indices_4,            &
                                                 quintic_indices_5,            &
                                                 quintic_indices_6,            &
                                                 ndf_w2v, undf_w2v, map_w2v,   &
                                                 ndf_wt, undf_wt, map_wt,      &
                                                 ndf_wc, undf_wc, map_wc )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_w2v
    integer(kind=i_def), intent(in)    :: undf_w2v
    integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
    integer(kind=i_def), intent(in)    :: ndf_wt
    integer(kind=i_def), intent(in)    :: undf_wt
    integer(kind=i_def), intent(in)    :: map_wt(ndf_wt)
    integer(kind=i_def), intent(in)    :: ndf_wc
    integer(kind=i_def), intent(in)    :: undf_wc
    integer(kind=i_def), intent(in)    :: map_wc(ndf_wc)
    real(kind=r_tran),   intent(in)    :: dep_dist_z(undf_w2v)
    real(kind=r_tran),   intent(in)    :: theta_height(undf_wt)
    real(kind=r_tran),   intent(inout) :: quintic_coef_1(undf_wc)
    real(kind=r_tran),   intent(inout) :: quintic_coef_2(undf_wc)
    real(kind=r_tran),   intent(inout) :: quintic_coef_3(undf_wc)
    real(kind=r_tran),   intent(inout) :: quintic_coef_4(undf_wc)
    real(kind=r_tran),   intent(inout) :: quintic_coef_5(undf_wc)
    real(kind=r_tran),   intent(inout) :: quintic_coef_6(undf_wc)
    integer(kind=i_def), intent(inout) :: quintic_indices_1(undf_wc)
    integer(kind=i_def), intent(inout) :: quintic_indices_2(undf_wc)
    integer(kind=i_def), intent(inout) :: quintic_indices_3(undf_wc)
    integer(kind=i_def), intent(inout) :: quintic_indices_4(undf_wc)
    integer(kind=i_def), intent(inout) :: quintic_indices_5(undf_wc)
    integer(kind=i_def), intent(inout) :: quintic_indices_6(undf_wc)

    ! Local arrays
    real(kind=r_tran)   :: displacement(nlayers+1)
    real(kind=r_tran)   :: frac_dist(nlayers+1)
    real(kind=r_tran)   :: sign_offset(nlayers+1)
    real(kind=r_tran)   :: z_arr(nlayers+1)
    real(kind=r_tran)   :: z_dep(nlayers+1)
    real(kind=r_tran)   :: z_dep_w3(nlayers)
    real(kind=r_tran)   :: z_arr_w3(nlayers)
    real(kind=r_tran)   :: dz(nlayers+ndf_wc-1)
    real(kind=r_tran)   :: coeff_local(nlayers+ndf_wc-1, 6)
    integer(kind=i_def) :: index_local(nlayers+ndf_wc-1, 6)
    integer(kind=i_def) :: int_disp(nlayers+1)
    integer(kind=i_def) :: k_dep(nlayers+1)
    integer(kind=i_def) :: k_dep_w3(nlayers)

    ! Local scalars
    integer(kind=i_def) :: wc_idx, k, j, nl, j_max, j_min, j_dep

    ! nl = nlayers for W3
    !    = nlayers+1 for Wtheta
    nl  = nlayers - 2 + ndf_wc

    wc_idx = map_wc(1)

    ! Extract departure distances and physical heights
    displacement(:) = dep_dist_z(map_w2v(1) : map_w2v(1)+nlayers)
    int_disp(:) = INT(displacement(:), i_def)
    frac_dist(:) = ABS(displacement(:) - REAL(int_disp(:), r_tran))
    sign_offset(:) = 0.5_r_tran*(1.0_r_tran + SIGN(1.0_r_tran, displacement(:)))
    z_arr(:) = theta_height(map_wt(1) : map_wt(1)+nlayers)

    ! Extract departure distances and physical heights
    displacement(:) = dep_dist_z(map_w2v(1) : map_w2v(1)+nlayers)
    int_disp(:) = INT(displacement(:), i_def)
    frac_dist(:) = ABS(displacement(:) - REAL(int_disp(:), r_tran))
    sign_offset(:) = 0.5_r_tran*(1.0_r_tran + SIGN(1.0_r_tran, displacement(:)))
    z_arr(:) = theta_height(map_wt(1) : map_wt(1)+nlayers)

    ! Wtheta departure heights and indices -------------------------------------
    ! Force bottom departure point to be zero
    k_dep(1) = 1
    z_dep(1) = z_arr(1)
    ! Calculate the index of the level below the departure distance, and the
    ! height of the departure point
    do k = 2, nlayers
      k_dep(k) = MAX(k - int_disp(k) - INT(sign_offset(k), i_def), 1)
      z_dep(k) = (                                                             &
        z_arr(k_dep(k)) * (                                                    &
          frac_dist(k)*sign_offset(k)                                          &
          + (1.0_r_tran - frac_dist(k))*(1.0_r_tran - sign_offset(k))          &
        )                                                                      &
      + z_arr(k_dep(k)+1) * (                                                  &
          (1.0_r_tran - frac_dist(k))*sign_offset(k)                           &
          + frac_dist(k)*(1.0_r_tran - sign_offset(k))                         &
        )                                                                      &
      )
    end do
    ! Force top departure point to be zero
    k_dep(nlayers+1) = nlayers
    z_dep(nlayers+1) = z_arr(nlayers+1)

    ! W3 departure heights and indices -----------------------------------------
    if (ndf_wc == 1) then
      ! W3 departure heights are the average of the Wtheta departure heights
      z_dep_w3(:) = 0.5_r_tran*(z_dep(1:nlayers) + z_dep(2:nlayers+1))
      z_arr_w3(:) = 0.5_r_tran*(z_arr(1:nlayers) + z_arr(2:nlayers+1))
      dz(1:nlayers-1) = z_arr_w3(2:nlayers) - z_arr_w3(1:nlayers-1)
      dz(nlayers) = dz(nlayers-1)  ! Copy top layer dz

      ! Bound W3 departure distances
      z_dep_w3(:) = MIN(z_arr_w3(nlayers), MAX(z_arr_w3(1), z_dep_w3(:)))

      ! Need to back out the indices of the corresponding levels
      ! Note that these aren't simply the average of the Wtheta indices
      k_dep_w3(1) = 1
      do k = 2, nlayers-1
        ! As first guesses, take the indices of the corresponding Wtheta dep pts
        j_max = MIN(MAX(k_dep(k), k_dep(k+1)), nlayers)
        j_min = MIN(k_dep(k), k_dep(k+1))
        j_dep = j_min
        ! Step downwards from upper guess to lower guess, to find the first
        ! model level that is below the W3 departure height. This gives the
        ! index to use for the interpolation
        do j = j_max, j_min, -1
          if (z_dep_w3(k) > z_arr_w3(j)) then
            j_dep = j
            EXIT
          end if
        end do
        k_dep_w3(k) = j_dep
      end do
      k_dep_w3(nlayers) = nlayers-1

    else
      ! Just need to compute layer depths
      dz(1:nlayers) = z_arr(2:nlayers+1) - z_arr(1:nlayers)
      dz(nlayers+1) = dz(nlayers)  ! Copy top layer dz
    end if

    ! Compute coefficients -----------------------------------------------------
    if ( ndf_wc == 1 ) then
      ! W3
      call compute_quintic_coeffs(                                             &
          k_dep_w3, z_dep_w3, z_arr_w3, dz, index_local, coeff_local, nlayers  &
      )
    else
      ! Wtheta
      call compute_quintic_coeffs(                                             &
          k_dep, z_dep, z_arr, dz, index_local, coeff_local, nlayers+1         &
      )
    end if

    ! Populate the global coefficient and index fields -------------------------
    quintic_coef_1(wc_idx : wc_idx+nl) = coeff_local(:,1)
    quintic_coef_2(wc_idx : wc_idx+nl) = coeff_local(:,2)
    quintic_coef_3(wc_idx : wc_idx+nl) = coeff_local(:,3)
    quintic_coef_4(wc_idx : wc_idx+nl) = coeff_local(:,4)
    quintic_coef_5(wc_idx : wc_idx+nl) = coeff_local(:,5)
    quintic_coef_6(wc_idx : wc_idx+nl) = coeff_local(:,6)
    quintic_indices_1(wc_idx : wc_idx+nl) = index_local(:,1)
    quintic_indices_2(wc_idx : wc_idx+nl) = index_local(:,2)
    quintic_indices_3(wc_idx : wc_idx+nl) = index_local(:,3)
    quintic_indices_4(wc_idx : wc_idx+nl) = index_local(:,4)
    quintic_indices_5(wc_idx : wc_idx+nl) = index_local(:,5)
    quintic_indices_6(wc_idx : wc_idx+nl) = index_local(:,6)

  end subroutine compute_vertical_quintic_coef_code

end module compute_vertical_quintic_coef_kernel_mod
