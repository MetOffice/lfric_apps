!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to cloud scheme.

module pc2_conv_coupling_kernel_mod

use argument_mod,      only: arg_type,              &
                             GH_FIELD, GH_REAL,     &
                             GH_READ, GH_READWRITE, &
                             GH_SCALAR, DOMAIN
use fs_continuity_mod, only: WTHETA
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: pc2_conv_coupling_kernel_type
  private
  type(arg_type) :: meta_args(15) = (/                         &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! theta_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! mv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! ml_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! cfl_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! cff_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! bcf_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! exner_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! pressure_inc_env
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dt_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dmv_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dmcl_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dcfl_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dcff_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dbcf_conv_wth
       arg_type(GH_SCALAR, GH_REAL, GH_READ)                    & ! dt
       /)
   integer :: operates_on = DOMAIN
contains
  procedure, nopass :: pc2_conv_coupling_code
end type

public :: pc2_conv_coupling_code

contains

!> @brief Interface to pc2 conv_coupling
!> @details Calculation of additional changes in cloud as a result of
!>          increments from a convection scheme. The convection scheme may have detrained cloud fraction
!>          and condensate and advected them down due to subsidence advection. Here the temperature,
!>          and humidity changes from convection scheme are used as forcing to calculate further cloud
!>          changes in the environment and we also calculate changes due to erosion.
!>          See UMDP 30 for more details.
!> @param[in]     nlayers       Number of layers
!> @param[in]     theta_wth     Potential temperature field
!> @param[in]     mv_wth        Vapour mass mixing ratio
!> @param[in]     ml_wth        Liquid cloud mass mixing ratio
!> @param[in]     cfl_wth       Liquid cloud fraction
!> @param[in]     cff_wth       Ice cloud fraction
!> @param[in]     bcf_wth       Bulk cloud fraction
!> @param[in]     exner_wth     Exner pressure in theta space
!> @param[in]     pressure_inc_env Environmental pressure inc from convection
!> @param[in,out] dt_conv_wth   Increment to theta in theta space
!> @param[in,out] dmv_conv_wth  Increment to vapour from convection in theta space
!> @param[in,out] dmcl_conv_wth Increment to liquid water content from convection in theta space
!> @param[in,out] dcfl_conv_wth Increment to liquid cloud fraction from convection in theta space
!> @param[in,out] dcff_conv_wth Increment to ice cloud fraction from convection in theta space
!> @param[in,out] dbcf_conv_wth Increment to bulk cloud fraction from convection in theta space
!> @param[in]     dt            The model timestep length
!> @param[in]     ndf_wth       Number of degrees of freedom per cell for theta space
!> @param[in]     undf_wth      Number of unique degrees of freedom in segment for theta space
!> @param[in]     map_wth       Dofmap for a segment of the theta space

subroutine pc2_conv_coupling_code( nlayers, seg_len,                           &
                                   ! Atmospheric fields
                                   theta_wth,                                  &
                                   mv_wth,                                     &
                                   ml_wth,                                     &
                                   cfl_wth,                                    &
                                   cff_wth,                                    &
                                   bcf_wth,                                    &
                                   exner_wth,                                  &
                                   pressure_inc_env,                           &
                                   ! Forcings in and combined responses out
                                   dt_conv_wth,                                &
                                   dmv_conv_wth,                               &
                                   dmcl_conv_wth,                              &
                                   dcfl_conv_wth,                              &
                                   dcff_conv_wth,                              &
                                   dbcf_conv_wth,                              &
                                   ! Other
                                   dt, ndf_wth, undf_wth, map_wth )

    use constants_mod, only: r_def, i_def, r_um, i_um

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use nlsizes_namelist_mod,       only: model_levels
    use pc2_hom_conv_mod,           only: pc2_hom_conv
    use cloud_inputs_mod,           only: dbsdtbs_turb_0, l_pc2_homog_conv_pressure
    use planet_constants_mod,       only: p_zero, kappa

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth, seg_len) :: map_wth

    ! Variables
    real(kind=r_def), intent(in),    dimension(undf_wth) :: theta_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: bcf_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: exner_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: pressure_inc_env

    ! The forcings (in) and the updated increments (out) as a result
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dt_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmcl_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcff_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dbcf_conv_wth

    ! The model timestep length
    real(kind=r_def), intent(in) :: dt

    ! Local variables
    real(r_um), dimension(seg_len,1) ::                                        &
                p_work,                                                        &
                ! Work arrays
                qv_work,  qcl_work,                                            &
                bcf_work, cfl_work, cff_work, t_work,                          &
                ! Forcings
                t_forcing, qv_forcing, cfl_forcing, p_forcing,                 &
                ! Increments
                qv_incr,  qcl_incr,                                            &
                bcf_incr, cfl_incr, t_incr,                                    &
                ! Other
                zeros

    integer(i_um) :: i, k
    integer(i_um), parameter :: j = 1

    do i = 1, seg_len
      zeros(i,j)=0.0_r_um
    end do

    do k = 1, model_levels-1

      if (l_pc2_homog_conv_pressure) then

        do i = 1, seg_len
          ! Pressure forcing - change in pressure from departure point
          p_forcing(i,j)  = pressure_inc_env(map_wth(1,i) + k)
          ! Pressure at departure point
          p_work(i,j) = p_zero*(exner_wth(map_wth(1,i) + k))                   &
                                        **(1.0_r_def/kappa) - p_forcing(i,j)

          ! Temperature forcing = theta after conv * exner (arrival - depart)
          t_forcing(i,j) = ( theta_wth(map_wth(1,i) + k) +                     &
                             dt_conv_wth(map_wth(1,i) + k) /                   &
                             exner_wth(map_wth(1,i) + k) ) *                   &
                           ( exner_wth(map_wth(1,i) + k) -                     &
                             (p_work(i,j)/p_zero)**kappa )
          ! Temperature at departure point
          t_work(i,j) = theta_wth(map_wth(1,i) + k) *                          &
                        exner_wth(map_wth(1,i) + k) +                          &
                        dt_conv_wth(map_wth(1,i) + k) - t_forcing(i,j)

          ! No forcing of vapour or cloud
          qv_forcing(i,j) = 0.0_r_um
          cfl_forcing(i,j)= 0.0_r_um

          ! Vapour after convection
          qv_work(i,j)    = mv_wth(map_wth(1,i) + k) + dmv_conv_wth(map_wth(1,i) + k)

          ! Cloud condensate and fraction after convection
          qcl_work(i,j)   = ml_wth(map_wth(1,i) + k) + dmcl_conv_wth(map_wth(1,i) + k)
          cfl_work(i,j)   = cfl_wth(map_wth(1,i) + k) + dcfl_conv_wth(map_wth(1,i) + k)
          cff_work(i,j)   = cff_wth(map_wth(1,i) + k) + dcff_conv_wth(map_wth(1,i) + k)
          bcf_work(i,j)   = bcf_wth(map_wth(1,i) + k) + dbcf_conv_wth(map_wth(1,i) + k)

          ! Output Increments from PC2
          t_incr(i,j)     = 0.0_r_um
          qv_incr(i,j)    = 0.0_r_um
          qcl_incr(i,j)   = 0.0_r_um
          cfl_incr(i,j)   = 0.0_r_um
          bcf_incr(i,j)   = 0.0_r_um
        end do

      else

        do i = 1, seg_len
          ! Pressure at centre of theta levels
          p_work(i,j) = p_zero*(exner_wth(map_wth(1,i) + k))                   &
                                        **(1.0_r_def/kappa)

          ! Temperature before convection
          t_work(i,j)    = theta_wth(map_wth(1,i) + k) *                       &
                           exner_wth(map_wth(1,i) + k)

          ! Vapour before convection
          qv_work(i,j)    = mv_wth(map_wth(1,i) + k)

          ! Cloud condensate and fraction after convection
          qcl_work(i,j)   = ml_wth(map_wth(1,i) + k) + dmcl_conv_wth(map_wth(1,i) + k)
          cfl_work(i,j)   = cfl_wth(map_wth(1,i) + k) + dcfl_conv_wth(map_wth(1,i) + k)
          cff_work(i,j)   = cff_wth(map_wth(1,i) + k) + dcff_conv_wth(map_wth(1,i) + k)
          bcf_work(i,j)   = bcf_wth(map_wth(1,i) + k) + dbcf_conv_wth(map_wth(1,i) + k)

          ! Forcings - convection increments except qcl
          t_forcing(i,j)  = dt_conv_wth(map_wth(1,i) + k)
          qv_forcing(i,j) = dmv_conv_wth(map_wth(1,i) + k)
          cfl_forcing(i,j)= dcfl_conv_wth(map_wth(1,i) + k)
          p_forcing(i,j)  = 0.0_r_um

          ! Output Increments from PC2
          t_incr(i,j)     = 0.0_r_um
          qv_incr(i,j)    = 0.0_r_um
          qcl_incr(i,j)   = 0.0_r_um
          cfl_incr(i,j)   = 0.0_r_um
          bcf_incr(i,j)   = 0.0_r_um
        end do

      end if

      call pc2_hom_conv(p_work,           & ! Pressure
                        real(dt,r_um),    & ! Model timestep in seconds
                        ! Variables
                        t_work,           & ! Temperature
                        qv_work,          & ! Water vapour
                        qcl_work,         & ! Liquid water content
                        bcf_work,         & ! Bulk cloud fraction
                        cfl_work,         & ! Liquid cloud fraction
                        cff_work,         & ! Ice cloud fraction
                        ! Forcings
                        t_forcing,        & ! Forcing of temperature
                        qv_forcing,       & ! Forcing of water vapour
                        zeros,            & ! Dummy dqclin forcing LWC
                        p_forcing,        & ! dpdt   forcing pressure
                        cfl_forcing,      & ! dcflin forcing lid cloud frac
                        ! Output variables
                        t_incr,           & ! Response to temperature
                        qv_incr,          & ! Response to water vapour
                        qcl_incr,         & ! Response to liquid water content
                        bcf_incr,         & ! Response to bulk cloud fraction
                        cfl_incr,         & ! Response liquid cloud fraction
                        ! Input variables (other quantities)
                        dbsdtbs_turb_0,   & ! pc2mixingrate dbsdtbs_turb_0
                        0.0_r_um          ) ! dbsdtbs1      dbsdtbs_turb_1

      ! Recast back to LFRic space
      do i = 1, seg_len
        dt_conv_wth  (map_wth(1,i) + k) = dt_conv_wth (map_wth(1,i) + k)  + t_incr(i,j)
        dmv_conv_wth (map_wth(1,i) + k) = dmv_conv_wth(map_wth(1,i) + k)  + qv_incr   (i,j)
        dmcl_conv_wth(map_wth(1,i) + k) = dmcl_conv_wth(map_wth(1,i) + k) + qcl_incr  (i,j)
        dcfl_conv_wth(map_wth(1,i) + k) = dcfl_conv_wth(map_wth(1,i) + k) + cfl_incr  (i,j)
        dbcf_conv_wth(map_wth(1,i) + k) = dbcf_conv_wth(map_wth(1,i) + k) + bcf_incr  (i,j)
      end do
    end do

    ! Set level 0 increment such that theta increment will equal level 1
    do i = 1, seg_len
      dt_conv_wth  (map_wth(1,i) + 0) = dt_conv_wth  (map_wth(1,i) + 1) &
                                      * exner_wth(map_wth(1,i) + 0)     &
                                      / exner_wth(map_wth(1,i) + 1)
      dmv_conv_wth (map_wth(1,i) + 0) = dmv_conv_wth (map_wth(1,i) + 1)
      dmcl_conv_wth(map_wth(1,i) + 0) = dmcl_conv_wth(map_wth(1,i) + 1)
      dcfl_conv_wth(map_wth(1,i) + 0) = dcfl_conv_wth(map_wth(1,i) + 1)
      dbcf_conv_wth(map_wth(1,i) + 0) = dbcf_conv_wth(map_wth(1,i) + 1)
    end do

end subroutine pc2_conv_coupling_code

end module pc2_conv_coupling_kernel_mod
