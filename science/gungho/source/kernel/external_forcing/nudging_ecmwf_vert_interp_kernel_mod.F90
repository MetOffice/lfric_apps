!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Performs vertical remapping of nudging reference fields
module nudging_ecmwf_vert_interp_kernel_mod

  use argument_mod,      only: arg_type,                                       &
                               GH_FIELD, GH_SCALAR,                            &
                               GH_REAL, GH_INTEGER,                            &
                               GH_WRITE, GH_READ,                              &
                               ANY_DISCONTINUOUS_SPACE_1,                      &
                               ANY_DISCONTINUOUS_SPACE_2,                      &
                               CELL_COLUMN

  use constants_mod,     only: r_def, i_def, l_def
  use fs_continuity_mod, only: Wtheta, W3
  use kernel_mod,        only: kernel_type

  implicit none

  private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
  type, public, extends(kernel_type) :: nudging_ecmwf_vert_interp_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, Wtheta),                     & ! Interpolated theta field
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3),                         & ! Interpolated U field
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3),                         & ! Interpolated V field
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  & ! Input T/theta field
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  & ! Input U field
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  & ! Input V field
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                     & ! Exner on wtheta
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                         & ! Exner on w3
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2),  & ! Nudging pressure field
        arg_type(GH_SCALAR, GH_REAL,    GH_READ),                              & ! p_zero
        arg_type(GH_SCALAR, GH_REAL,    GH_READ),                              & ! kappa
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                               & ! nudge_levels
    /)
    integer :: operates_on = CELL_COLUMN
contains
    procedure, nopass :: nudging_ecmwf_vert_interp_code
end type nudging_ecmwf_vert_interp_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: nudging_ecmwf_vert_interp_code

contains
  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  !> @brief Performs vertical remapping of nudging reference fields
  !> @param[in]     nlayers         Number of layers in model mesh
  !> @param[in,out] theta_mod_levs  theta field in Wtheta
  !> @param[in,out] u_mod_levs      zonal wind at W3 on model levels
  !> @param[in,out] v_mod_levs      meridional wind at W3 on model levels
  !> @param[in]     T_nudge_levs    input nudging T or theta field
  !> @param[in]     u_nudge_levs    input nudging zonal wind
  !> @param[in]     v_nudge_levs    input nudging v wind
  !> @param[in]     exner_in_wth    exner pressure in Wtheta
  !> @param[in]     exner_in_w3     exner pressure in W3
  !> @param[in]     ln_surfp_nudge  surface pressure from nudging data
  !> @param[in]     p_zero          Reference pressure for Exner calculation
  !> @param[in]     kappa           Rd/Cp, used in definition of Exner
  !> @param[in]     nudge_levels    number of levels in reference data
  !> @param[in]     ndf_wth         number of degrees of freedom per cell
  !> @param[in]     undf_wth        num of DoFs in this partition for Wtheta
  !> @param[in]     map_wth         index of DoFs for lowest cell for Wtheta
  !> @param[in]     ndf_w3          number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3         num of DoFs in this partition for W3
  !> @param[in]     map_w3          index of DoFs for lowest cell for W3
  !> @param[in]     ndf_ec          num DoFs per cell for input nudging fields
  !> @param[in]     undf_ec         num DoFs in this partition for nudging field
  !> @param[in]     map_ec          dofmap for nudging fields
  !> @param[in]     ndf_ec_2d       num DoFs per cell for surface pressure
  !> @param[in]     undf_ec_2d      num DoFs in this partition for surf pressure
  !> @param[in]     map_ec_2d       dofmap for surf pressure field
  subroutine nudging_ecmwf_vert_interp_code(                                   &
          nlayers,                                                             &
          theta_mod_levs, u_mod_levs, v_mod_levs,                              &
          T_nudge_levs, u_nudge_levs, v_nudge_levs,                            &
          exner_in_wth, exner_in_w3, ln_surfp_nudge,                           &
          p_zero, kappa, nudge_levels,                                         &
          ndf_wth, undf_wth, map_wth,                                          &
          ndf_w3, undf_w3, map_w3,                                             &
          ndf_ec, undf_ec, map_ec,                                             &
          ndf_ec2d, undf_ec2d, map_ec2d                                        &
  )

    use ecmwf_level_support_mod, only: ak_ec_137, bk_ec_137,                   &
                                       ak_ec_137_50, bk_ec_137_50

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in)    :: map_wth(ndf_wth)
    integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in)    :: ndf_ec, undf_ec
    integer(kind=i_def), intent(in)    :: map_ec(ndf_ec)
    integer(kind=i_def), intent(in)    :: ndf_ec2d, undf_ec2d
    integer(kind=i_def), intent(in)    :: map_ec2d(ndf_ec2d)
    integer(kind=i_def), intent(in)    :: nudge_levels

    real(kind=r_def),    intent(inout) :: theta_mod_levs(undf_wth)
    real(kind=r_def),    intent(inout) :: u_mod_levs(undf_w3)
    real(kind=r_def),    intent(inout) :: v_mod_levs(undf_w3)
    real(kind=r_def),    intent(inout) :: T_nudge_levs(undf_ec)
    real(kind=r_def),    intent(in)    :: u_nudge_levs(undf_ec)
    real(kind=r_def),    intent(in)    :: v_nudge_levs(undf_ec)
    real(kind=r_def),    intent(in)    :: exner_in_wth(undf_wth)
    real(kind=r_def),    intent(in)    :: exner_in_w3(undf_w3)
    real(kind=r_def),    intent(in)    :: ln_surfp_nudge(undf_ec2d)
    real(kind=r_def),    intent(in)    :: p_zero, kappa

    ! Local variables
    ! Pointers to height coefficients for the current dataset
    real(kind=r_def) :: ak_ec(0:nudge_levels)
    real(kind=r_def) :: bk_ec(0:nudge_levels)

    ! Coefficients interpolated on full levels
    real(kind=r_def), dimension(nudge_levels) :: ak_int, bk_int

    ! Pressure profiles
    real(kind=r_def), dimension(nudge_levels) :: p_nudge_levels
    real(kind=r_def), dimension(nlayers)      :: press_in_wth
    real(kind=r_def), dimension(0:nlayers-1)  :: press_in_w3

    ! Interim copies of nudging variables - for inverting data
    real(kind=r_def), dimension(nudge_levels) :: theta_nudge_int
    real(kind=r_def), dimension(nudge_levels) :: u_nudge_int
    real(kind=r_def), dimension(nudge_levels) :: v_nudge_int
    real(kind=r_def) :: step_wth, step_w3, delta

    ! Loop counters
    integer(kind=i_def) :: ec_idx, wt_idx, w3_idx
    integer(kind=i_def) :: j, k, l, m
    integer(kind=i_def) :: kth, lth, mth

    ! Set height parameters, based on the number of nudging levels
    select case (nudge_levels)
    case (137)
      ak_ec(:) = ak_ec_137(:)
      bk_ec(:) = bk_ec_137(:)
    case (88)
      ak_ec(:) = ak_ec_137_50(:)
      bk_ec(:) = bk_ec_137_50(:)
    end select

    ec_idx = map_ec(1)
    w3_idx = map_w3(1)
    wt_idx = map_wth(1)

    ! Contruct height parameters on full levels
    do k = 1, nudge_levels

      ! Convert the levels to the same order as model
      ! (i.e  start at the bottom of the atmosphere)
      j = nudge_levels + 1 - k

      ak_int(k) = 0.5_r_def*(ak_ec(j-1) + ak_ec(j))
      bk_int(k) = 0.5_r_def*(bk_ec(j-1) + bk_ec(j))

      ! construct pressure profile for nudging data
      p_nudge_levels(k) = ak_int(k)  &
                          + bk_int(k) * exp(ln_surfp_nudge( map_ec2d(1)))
    end do

    ! Construct pressure profiles on model levels
    do k = 0, nlayers - 1
      kth = k + 1
      press_in_wth(kth) = p_zero * exner_in_wth(wt_idx+kth)**(1.0_r_def/kappa)
      press_in_w3(k) = p_zero * exner_in_w3(w3_idx+k)**(1.0_r_def/kappa)
    end do

    ! Copy to interim arrays
    ! Convert T to theta using the pressure profile generated here.
    ! Input field on W3 mesh goes from 0 -> nudge_levels-1
    do k = 0, nudge_levels - 1
      theta_nudge_int(k+1) =                                                   &
          T_nudge_levs(ec_idx+k) * (p_zero / p_nudge_levels(k+1))**kappa
      u_nudge_int(k+1) = u_nudge_levs(ec_idx+k)
      v_nudge_int(k+1) = v_nudge_levs(ec_idx+k)
    end do

    ! Now interpolate from nudging levels to model levels
    ! Initialise the data loop
    l = 1
    m = l + 1
    lth = 1
    mth = lth + 1

    ! Increment L until Nudging Pressures straddle 1st model Pressure level
    ! Separately for wth and w3 grids
    do while(press_in_w3(1) < p_nudge_levels(m) .and. m < nudge_levels)
      l = l + 1
      m = l + 1
    end do

    do while (press_in_wth(1) < p_nudge_levels(mth) .and. mth < nudge_levels)
      lth = lth + 1
      mth = lth + 1
    end do

    ! Loop over model levels interpolating on to these levels
    do k = 0, nlayers - 1
      kth = k + 1

      ! Calculate interpolation weights in log-pressure space
      step_w3 = log(press_in_w3(k)) - log(p_nudge_levels(l))
      step_w3 = step_w3/(log(p_nudge_levels(l+1)) - log(p_nudge_levels(l)))

      step_wth = log(press_in_wth(kth)) - log(p_nudge_levels(lth))
      step_wth = step_wth/(log(p_nudge_levels(lth+1)) - log(p_nudge_levels(lth)))

      ! Limit large extrapolations - keep step between -1.0 <-> 2.0
      step_w3 = max( min(step_w3, 2.0_r_def), -1.0_r_def )
      step_wth = max( min(step_wth, 2.0_r_def), -1.0_r_def )

      ! Interpolate the analyses on to the model levels
      delta = u_nudge_int(m) - u_nudge_int(l)
      u_mod_levs(w3_idx+k) = u_nudge_int( l ) + delta * step_w3

      delta = v_nudge_int(m) - v_nudge_int(l)
      v_mod_levs(w3_idx+k) = v_nudge_int(l) + delta * step_w3

      delta = theta_nudge_int(mth) - theta_nudge_int(lth)
      theta_mod_levs(wt_idx+kth) = theta_nudge_int(lth) + delta * step_wth

      ! If the pressure of the next nudge_level is greater than the next
      ! model level then increment to the nudge level
      ! Also check that we have not exceeded the number of nudge or model
      ! levels
      if (kth < nlayers) then
        mth = lth + 1
        do while (press_in_wth(kth+1) < p_nudge_levels(mth) .and. mth < nudge_levels)
          lth = lth + 1
          mth = lth + 1
        end do
      end if

      if (k < nlayers - 1) then
        m = l + 1
        do while (press_in_w3(k+1) < p_nudge_levels(m) .and. m < nudge_levels)
          l = l + 1
          m = l + 1
        end do
      end if

    end do      ! Loop over levels

    ! Set bottom value of theta to the value of the next level up
    theta_mod_levs(wt_idx) = theta_mod_levs(wt_idx+1)

  end subroutine nudging_ecmwf_vert_interp_code

end module nudging_ecmwf_vert_interp_kernel_mod
