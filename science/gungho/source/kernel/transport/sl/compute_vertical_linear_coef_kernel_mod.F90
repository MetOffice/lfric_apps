!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to compute the vertical linear interpolation coefficients.
!> @details Computes the linear vertical interpolation coefficients used by the
!!          semi-Lagrangian transport scheme.

module compute_vertical_linear_coef_kernel_mod

  use argument_mod,          only : arg_type,                      &
                                    GH_FIELD, GH_WRITE, GH_SCALAR, &
                                    GH_REAL, GH_READ, GH_LOGICAL,  &
                                    GH_INTEGER, CELL_COLUMN,       &
                                    ANY_DISCONTINUOUS_SPACE_1,     &
                                    STENCIL, CROSS2D
  use fs_continuity_mod,     only : W2v, Wtheta
  use constants_mod,         only : r_tran, i_def, l_def, EPS_R_TRAN
  use kernel_mod,            only : kernel_type
  use sl_support_mod,        only : compute_linear_coeffs
  use log_mod,               only : log_event, LOG_LEVEL_ERROR

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed
  !>                                      by the PSy layer.
  type, public, extends(kernel_type) :: compute_vertical_linear_coef_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v,    STENCIL(CROSS2D)),  & ! dep_dist_z
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta, STENCIL(CROSS2D)),  & ! theta_height
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! linear_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! linear_coef
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: compute_vertical_linear_coef_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: compute_vertical_linear_coef_code

  contains

  !-------------------------------------------------------------------------------
  !> @details This kernel calculates the departure point of w/theta-points using
  !!          only w (i.e., vertical motion only), then interpolate theta at the
  !!          departure point using 1d-linear-Lagrange interpolation.
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     dep_dist_z    The vertical departure distance
  !> @param[in]     d_stencil_sizes Sizes of branches of the dep_dist_z
  !!                                 cross stencil
  !> @param[in]     d_stencil_max  Maximum size of a dep_dist_z stencil branch
  !> @param[in]     d_stencil_map  Dofmap for the dep_dist_z stencil
  !> @param[in]     theta_height  The height of theta-points
  !> @param[in]     t_stencil_sizes Sizes of branches of the theta_height
  !!                                 cross stencil
  !> @param[in]     t_stencil_max  Maximum size of a theta_height stencil branch
  !> @param[in]     t_stencil_map  Dofmap for the theta_height stencil
  !> @param[in,out] linear_coef   The linear interpolation coefficients (1-2)
  !> @param[in]     ndf_w2v       Num DoFs per cell for W2v
  !> @param[in]     undf_w2v      Num DoFs for this rank for W2v
  !> @param[in]     map_w2v       DoF map for W2v
  !> @param[in]     ndf_wt        Num DoFs per cell for Wtheta
  !> @param[in]     undf_wt       Num DoFs for this rank for Wtheta
  !> @param[in]     map_wt        DoF map for Wtheta
  !> @param[in]     ndf_wc        Num DoFs per cell for the coefficients
  !> @param[in]     undf_wc       Num DoFs for this rank for the coefficients
  !> @param[in]     map_wc        DoF map for the coefficients
  !-------------------------------------------------------------------------------
  subroutine compute_vertical_linear_coef_code( nlayers,                       &
                                                dep_dist_z,                    &
                                                d_stencil_sizes,               &
                                                d_stencil_max,                 &
                                                d_stencil_map,                 &
                                                theta_height,                  &
                                                t_stencil_sizes,               &
                                                t_stencil_max,                 &
                                                t_stencil_map,                 &
                                                linear_coef_1, linear_coef_2,  &
                                                ndf_w2v, undf_w2v, map_w2v,    &
                                                ndf_wt, undf_wt, map_wt,       &
                                                ndf_wc, undf_wc, map_wc )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_w2v
    integer(kind=i_def), intent(in)    :: undf_w2v
    integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
    integer(kind=i_def), intent(in)    :: d_stencil_max
    integer(kind=i_def), intent(in)    :: d_stencil_sizes(4)
    integer(kind=i_def), intent(in)    :: d_stencil_map(ndf_w2v,d_stencil_max,4)
    integer(kind=i_def), intent(in)    :: ndf_wt
    integer(kind=i_def), intent(in)    :: undf_wt
    integer(kind=i_def), intent(in)    :: map_wt(ndf_wt)
    integer(kind=i_def), intent(in)    :: t_stencil_max
    integer(kind=i_def), intent(in)    :: t_stencil_sizes(4)
    integer(kind=i_def), intent(in)    :: t_stencil_map(ndf_wt,t_stencil_max,4)
    integer(kind=i_def), intent(in)    :: ndf_wc
    integer(kind=i_def), intent(in)    :: undf_wc
    integer(kind=i_def), intent(in)    :: map_wc(ndf_wc)
    real(kind=r_tran),   intent(in)    :: dep_dist_z(undf_w2v)
    real(kind=r_tran),   intent(in)    :: theta_height(undf_wt)
    real(kind=r_tran),   intent(inout) :: linear_coef_1(undf_wc)
    real(kind=r_tran),   intent(inout) :: linear_coef_2(undf_wc)

    ! Local arrays
    real(kind=r_tran)   :: displacement(nlayers+1)
    real(kind=r_tran)   :: frac_dist(nlayers+1)
    real(kind=r_tran)   :: sign_offset(nlayers+1)
    real(kind=r_tran)   :: z_arr(nlayers+1)
    real(kind=r_tran)   :: z_dep(nlayers+1)
    real(kind=r_tran)   :: z_dep_w3(nlayers)
    real(kind=r_tran)   :: z_arr_w3(nlayers)
    real(kind=r_tran), allocatable :: dz(:)
    real(kind=r_tran), allocatable :: coeff_local(:,:)
    integer(kind=i_def) :: int_disp(nlayers+1)
    integer(kind=i_def) :: k_dep(nlayers+1)
    integer(kind=i_def) :: k_dep_w3(nlayers)

    ! Local scalars
    integer(kind=i_def) :: wc_idx, k, j, nl, j_max, j_min, j_dep, d, ndofs, df

    select case ( ndf_wc )
      case ( 1 ) ! W3
        nl = nlayers - 1
        ndofs = 1
      case ( 2 ) ! Wtheta
        nl = nlayers
        ndofs = 1
      case ( 4 ) ! W2h
        nl = nlayers - 1
        ndofs = 4
      case default
        call log_event('Invlaid space in compute_vertical_linear_coef', LOG_LEVEL_ERROR)
    end select

    allocate( dz(nl+1), coeff_local(nl+1, 2) )

    do df = 1, ndofs
      wc_idx = map_wc(df)

      ! Extract departure distances and physical heights
      if ( ndofs == 4 ) then
        ! Average dep_dist to edge midpoint corresponding to a u or v point
        ! but still on the wtheta layers, i.e. in 2d wt is wtheta point,
        ! u is w2 point and d is new dep_dist point
        ! d --- wt --- d --- wt --- d
        ! |            |            |
        ! u            u            u
        ! |            |            |
        ! d --- wt --- d --- wt --- d
        d = min(2,d_stencil_sizes(df))
        displacement(:) = 0.5_r_tran * (                                      &
            dep_dist_z(d_stencil_map(1,1,1)  : d_stencil_map(1,1,1)+nlayers)  &
          + dep_dist_z(d_stencil_map(1,d,df) : d_stencil_map(1,d,df)+nlayers))
        z_arr(:) = 0.5_r_tran* (                                                &
            theta_height(t_stencil_map(1,1,1)  : t_stencil_map(1,1,1) +nlayers) &
          + theta_height(t_stencil_map(1,d,df) : t_stencil_map(1,d,df)+nlayers) )
      else
        ! Transporting coordinates in W3 or Wtheta
        displacement(:) = dep_dist_z(map_w2v(1) : map_w2v(1)+nlayers)
        z_arr(:) = theta_height(map_wt(1) : map_wt(1)+nlayers)
      end if

      int_disp(:) = INT(displacement(:), i_def)
      frac_dist(:) = ABS(displacement(:) - REAL(int_disp(:), r_tran))
      sign_offset(:) = 0.5_r_tran*(1.0_r_tran + SIGN(1.0_r_tran, displacement(:)))

      ! Wtheta departure heights and indices -------------------------------------
      ! Force bottom departure point to be zero
      k_dep(1) = 1
      z_dep(1) = z_arr(1)
      ! Calculate the index of the level below the departure distance, and the
      ! height of the departure point
      do k = 2, nlayers
        k_dep(k) = MAX(k - int_disp(k) - INT(sign_offset(k), i_def), 1)
        z_dep(k) = (                                                           &
          z_arr(k_dep(k)) * (                                                  &
            frac_dist(k)*sign_offset(k)                                        &
            + (1.0_r_tran - frac_dist(k))*(1.0_r_tran - sign_offset(k))        &
          )                                                                    &
        + z_arr(k_dep(k)+1) * (                                                &
            (1.0_r_tran - frac_dist(k))*sign_offset(k)                         &
            + frac_dist(k)*(1.0_r_tran - sign_offset(k))                       &
          )                                                                    &
        )
      end do
      ! Force top departure point to be zero
      k_dep(nlayers+1) = nlayers
      z_dep(nlayers+1) = z_arr(nlayers+1)

      if (ndf_wc == 2) then
        ! Wtheta departure heights and indices ---------------------------------
        ! Just need to compute layer depths
        dz(1:nlayers) = z_arr(2:nlayers+1) - z_arr(1:nlayers)
        dz(nlayers+1) = dz(nlayers)  ! Copy top layer dz

        call compute_linear_coeffs(                                            &
            k_dep, z_dep, z_arr, dz, coeff_local, nlayers+1                    &
        )
      else
        ! W3 & W2h departure heights and indices -------------------------------
        ! Departure heights are the average of the Wtheta departure heights
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

        call compute_linear_coeffs(                                            &
            k_dep_w3, z_dep_w3, z_arr_w3, dz, coeff_local, nlayers             &
        )

      end if

      ! Populate the global coefficient fields -----------------------------------
      linear_coef_1(wc_idx : wc_idx+nl) = coeff_local(:,1)
      linear_coef_2(wc_idx : wc_idx+nl) = coeff_local(:,2)
    end do

    deallocate( dz, coeff_local )

  end subroutine compute_vertical_linear_coef_code

end module compute_vertical_linear_coef_kernel_mod
