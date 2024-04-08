!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the UM Lambert-Lewis convection scheme.
!>
module conv_ll_kernel_mod

  use argument_mod,           only : arg_type,          &
                                     GH_FIELD, GH_REAL, &
                                     GH_READ, GH_WRITE, &
                                     CELL_COLUMN,       &
                                     ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,          only : i_def, i_um, r_def, r_um
  use fs_continuity_mod,      only : W3, Wtheta
  use kernel_mod,             only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: conv_ll_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                  &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                   & ! dt_conv
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                   & ! dmv_conv
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                   & ! dmcl_conv
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                   & ! dcfl_conv
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! theta_star
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! m_v
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! m_cl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                       & ! exner_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1) & ! conv_rain
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: conv_ll_code
  end type conv_ll_kernel_type

  public :: conv_ll_code

contains

  !> @brief Interface to the Lambert-Lewis convection scheme
  !> @details The Lambert-Lewis convection scheme is a simple
  !>           convection parametrization that mixes theta and moisture
  !>           as documented in UMDP41
  !> @param[in]     nlayers      Number of layers
  !> @param[in,out] dt_conv      Convection temperature increment
  !> @param[in,out] dmv_conv     Convection vapour increment
  !> @param[in,out] dmcl_conv    Convection liquid increment
  !> @param[in,out] dcfl_conv    Convection liquid cloud fraction increment
  !> @param[in]     theta_star   Potential temperature predictor after advection
  !> @param[in]     m_v          Vapour mixing ratio after advection
  !> @param[in]     m_cl         Cloud liquid mixing ratio after advection
  !> @param[in]     exner_in_wth Exner pressure field in wth space
  !> @param[in]     exner_in_w3  Exner pressure field in density space
  !> @param[in,out] conv_rain_2d Convective rain from twod fields
  !> @param[in]     ndf_wth      Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth     Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth      Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3       Number of DOFs per cell for density space
  !> @param[in]     undf_w3      Number of unique DOFs for density space
  !> @param[in]     map_w3       Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d       Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d      Number of unique DOFs for 2D fields
  !> @param[in]     map_2d       Dofmap for the cell at the base of the column for 2D fields
  subroutine conv_ll_code(nlayers,      &
                          dt_conv,      &
                          dmv_conv,     &
                          dmcl_conv,    &
                          dcfl_conv,    &
                          theta_star,   &
                          m_v,          &
                          m_cl,         &
                          exner_in_wth, &
                          exner_in_w3,  &
                          conv_rain_2d, &
                          ndf_wth,      &
                          undf_wth,     &
                          map_wth,      &
                          ndf_w3,       &
                          undf_w3,      &
                          map_w3,       &
                          ndf_2d,       &
                          undf_2d,      &
                          map_2d)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use llcs, only: llcs_control
    use nlsizes_namelist_mod, only: row_length, rows
    use planet_constants_mod, only: p_zero, kappa

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3, ndf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_2d

    integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_2d),  intent(in) :: map_2d

    real(kind=r_def), dimension(undf_wth), intent(inout)  :: dt_conv, dmv_conv, &
                                                             dmcl_conv, dcfl_conv

    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_star, &
                                                           m_v, m_cl,  &
                                                           exner_in_wth

    real(kind=r_def), dimension(undf_w3),  intent(in)   :: exner_in_w3

    real(kind=r_def), dimension(undf_2d),  intent(inout)  :: conv_rain_2d

    ! Local variables for the kernel
    integer(kind=i_def) :: k

    real(r_um), dimension(row_length,rows,nlayers) :: theta_conv, q_conv, &
         theta_inc, q_inc, qcl_inc, cf_liquid_inc, p_theta_levels, qcl_conv
    real(r_um), dimension(row_length,rows,0:nlayers) :: p_rho_minus_one
    real(r_um), dimension(row_length,rows) ::  conv_rain, p_star

    !-----------------------------------------------------------------------
    ! Initialise variables required from input fields
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      theta_conv(1,1,k) = theta_star(map_wth(1) + k)
      q_conv(1,1,k) = m_v(map_wth(1) + k)
      qcl_conv(1,1,k) = m_cl(map_wth(1) + k)
      p_theta_levels(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
    end do
    ! Set p_rho_minus_one at lowest level with surface pressure
    p_rho_minus_one(1,1,0) = p_zero*(exner_in_wth(map_wth(1)))**(1.0_r_def/kappa)
    do k = 1, nlayers - 1
      ! Pressure on rho levels, without level 1
      p_rho_minus_one(1,1,k) = p_zero*(exner_in_w3(map_w3(1) + k))**(1.0_r_def/kappa)
    end do
    ! Initialised p_rho_minus_one at top of atmosphere
    p_rho_minus_one(1,1,nlayers) = 0.0_r_um
    ! Initialise increments and output fields to zero
    theta_inc(:,:,:) = 0.0_r_um
    q_inc(:,:,:) = 0.0_r_um
    qcl_inc(:,:,:) = 0.0_r_um
    cf_liquid_inc(:,:,:) = 0.0_r_um
    conv_rain(:,:) = 0.0_r_um

    !-----------------------------------------------------------------------
    ! Call the convection scheme
    !-----------------------------------------------------------------------
    call llcs_control(row_length, rows, theta_conv, q_conv, qcl_conv, &
                      p_rho_minus_one, p_theta_levels, theta_inc,     &
                      q_inc, qcl_inc, cf_liquid_inc, conv_rain)

    !-----------------------------------------------------------------------
    ! Update fields to pass out
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      ! "Increments" passed out are actually updated full fields, so
      ! calculate the real increment
      theta_inc(1,1,k) = theta_inc(1,1,k) - theta_conv(1,1,k)
      q_inc(1,1,k) = q_inc(1,1,k) - q_conv(1,1,k)
      ! Increments to pass out
      dt_conv(map_wth(1) + k) = theta_inc(1,1,k)*exner_in_wth(map_wth(1) + k)
      dmv_conv(map_wth(1) + k) = q_inc(1,1,k)
      dmcl_conv(map_wth(1) + k) = qcl_inc(1,1,k)
      dcfl_conv(map_wth(1) + k) = cf_liquid_inc(1,1,k)
     end do
    ! Set level 0 increment such that theta increment will equal level 1
    dt_conv(map_wth(1)) = theta_inc(1,1,1)*exner_in_wth(map_wth(1))
    dmv_conv(map_wth(1)) = dmv_conv(map_wth(1) + 1)
    dmcl_conv(map_wth(1)) = dmcl_conv(map_wth(1) + 1)
    dcfl_conv(map_wth(1)) = dcfl_conv(map_wth(1) + 1)

    ! Copy conv_rain
    conv_rain_2d(map_2d(1))  = conv_rain(1,1)

  end subroutine conv_ll_code

end module conv_ll_kernel_mod
