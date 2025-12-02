!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to cloud scheme.

module pc2_homog_iau_kernel_mod

use argument_mod,      only: arg_type,          &
                             GH_FIELD, GH_REAL, &
                             GH_READ, GH_WRITE, &
                             GH_SCALAR, CELL_COLUMN
use fs_continuity_mod, only: WTHETA
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: pc2_homog_iau_kernel_type
  private
  type(arg_type) :: meta_args(17) = (/                 &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! mv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! ml_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! cfl_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! cff_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! bcf_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! theta_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! exner_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! exner_earlier_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! dtheta_forcing_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! dqv_forcing_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! dqcl_forcing_wth
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! dtheta_inc
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! dqv_inc_wth
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! dqcl_inc_wth
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! dcfl_inc_wth
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! dbcf_inc_wth
       arg_type(GH_SCALAR, GH_REAL, GH_READ)           & ! dt
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: pc2_homog_iau_code
end type

public :: pc2_homog_iau_code

contains

!> @brief Interface to the pc2 homog_iau forcing code
!> @param[in]     nlayers            Number of layers
!> @param[in]     mv_wth             Vapour mass mixing ratio
!> @param[in]     ml_wth             Liquid cloud mass mixing ratio
!> @param[in]     cfl_wth            Liquid cloud fraction
!> @param[in]     cff_wth            Ice cloud fraction
!> @param[in]     bcf_wth            Bulk cloud fraction
!> @param[in]     theta_wth          Potential temperature field
!> @param[in]     exner_wth          Exner pressure in potential temperature space
!> @param[in]     exner_earlier_wth  Exner pressure in potential temperature space
!> @param[in]     dtheta_forcing_wth Potential temperature forcing
!> @param[in]     dqv_forcing_wth    Water vapour forcing
!> @param[in]     dqcl_forcing_wth   Liquid water content forcing
!> @param[in,out] dtheta_inc_wth     Increment to potential temperature
!> @param[in,out] dqv_inc_wth        Increment to water vapour
!> @param[in,out] dqcl_inc_wth       Increment to liquid water content
!> @param[in,out] dcfl_inc_wth       Increment to liquid cloud fraction
!> @param[in,out] dbcf_inc_wth       Increment to bulk cloud fraction
!> @param[in]     dt                 The model timestep length
!> @param[in]     ndf_wth            Number of degrees of freedom per cell for
!!                                    potential temperature space
!> @param[in]     undf_wth           Number of unique of degrees of freedom
!!                                    for potential temperature space
!> @param[in]     map_wth            Dofmap for the cell at the base of the column
!!                                    for potential temperature space

subroutine pc2_homog_iau_code( nlayers,                    &
                               mv_wth,                     &
                               ml_wth,                     &
                               cfl_wth,                    &
                               cff_wth,                    &
                               bcf_wth,                    &
                               theta_wth,                  &
                               exner_wth,                  &
                               exner_earlier_wth,          &
                               ! Forcings
                               dtheta_forcing_wth,         &
                               dqv_forcing_wth,            &
                               dqcl_forcing_wth,           &
                               ! Response increments
                               dtheta_inc_wth,             &
                               dqv_inc_wth,                &
                               dqcl_inc_wth,               &
                               dcfl_inc_wth,               &
                               dbcf_inc_wth,               &
                               ! Other
                               dt,                         &
                               ndf_wth, undf_wth, map_wth  )

    use constants_mod, only: r_def, i_def, r_um, i_um, r_second

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use nlsizes_namelist_mod,       only: row_length, rows, model_levels
    use pc2_homog_plus_turb_mod,    only: pc2_homog_plus_turb
    use planet_constants_mod,       only: p_zero, kappa
    use gen_phys_inputs_mod,        only: l_mr_physics

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    real(kind=r_second),    intent(in) :: dt

    real(kind=r_def), intent(in),  dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: bcf_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: theta_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_earlier_wth

    ! The forcings coming IN
    real(kind=r_def), intent(in),  dimension(undf_wth) :: dtheta_forcing_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: dqv_forcing_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: dqcl_forcing_wth

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    ! The changes to the fields as a result
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dtheta_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dqv_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dqcl_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dbcf_inc_wth

    real(r_um), dimension(row_length,rows,model_levels) :: &
                  qv_work, qcl_work,                       &
                  cfl_work, cff_work, bcf_work,            &
                  t_work, theta_work, pressure

    real(r_um), dimension(row_length,rows,model_levels) :: &
         dtdt, dqdt, dldt, dpdt

    real(r_um), dimension(row_length,rows,model_levels) :: &
      t_earlier, p_earlier

    integer(i_um) :: k

    ! Hardwired 0.0, to be used for dbsdtbs0 and dbsdtbs1
    real(r_um), parameter :: zero = 0.0_r_um

    !-----------------------------------------------------------------------
    ! Initialisation of prognostic variables and arrays
    !-----------------------------------------------------------------------
    do k = 1, model_levels
      ! Calculate temperature from theta.
      t_earlier(1,1,k)   = theta_wth(map_wth(1) + k) *            &
                           exner_earlier_wth(map_wth(1) + k)

      ! Heating/cooling forcing has come in as potential temperature (theta)
      ! forcing, so calculate a provisional new temperature...
      t_work(1,1,k) = ( theta_wth(map_wth(1) + k) + dtheta_forcing_wth(map_wth(1) + k) ) * &
                        exner_wth(map_wth(1) + k)
      ! ... and calculate change in temperature (T) from that.
      ! This will be a forcing for PC2.
      dtdt(1,1,k) = t_work(1,1,k) - t_earlier(1,1,k)

      ! Calculate pressures
      p_earlier(1,1,k)  = p_zero * ( exner_earlier_wth(map_wth(1) + k) ) &
               **(1.0_r_def/kappa)

      ! Current pressure
      pressure(1,1,k) = p_zero*(exner_wth(map_wth(1) + k))     &
                        **(1.0_r_def/kappa)

      ! The pressure forcing used for PC2
      dpdt(1,1,k) = pressure(1,1,k) - p_earlier(1,1,k)

      ! Other forcings used by PC2
      dqdt(1,1,k) = dqv_forcing_wth (map_wth(1) + k)
      dldt(1,1,k) = dqcl_forcing_wth(map_wth(1) + k)

      ! Moist prognostics
      qv_work(1,1,k)   = mv_wth(map_wth(1) + k)
      qcl_work(1,1,k)  = ml_wth(map_wth(1) + k)

      ! Cast LFRic cloud fractions onto work arrays.
      cfl_work(1,1,k) = cfl_wth(map_wth(1) + k)
      cff_work(1,1,k) = cff_wth(map_wth(1) + k)
      bcf_work(1,1,k) = bcf_wth(map_wth(1) + k)
    end do

    ! Call to pc2_homog_plus_turb routine.
    call pc2_homog_plus_turb(p_earlier,          & ! Pressure related fields
                             model_levels,       & ! levels in the vertical
                             dt,                 & ! Timestep
                                                   ! Fields that get updated:
                             t_earlier,          & !   temperature
                             bcf_work,           & !   bulk cloud fraction
                             cfl_work,           & !   liquid cloud fraction
                             cff_work,           & !   ice cloud fraction
                             qv_work,            & !   vapour
                             qcl_work,           & !   liquid water content
                                                   ! Forcings quantities for PC2:
                             dtdt,               & !   temperature forcing
                             dqdt,               & !   vapour forcing
                             dldt,               & !   liquid water content forcing
                             dpdt,               & !   pressure forcing
                                                   ! Quantities for PC2 erosion:
                             zero,               & !    dbsdtbs0
                             zero,               & !    dbsdtbs1
                                                   ! Model switches
                             l_mr_physics )        !   mixing ratio

    do k = 1, model_levels

      ! New theta found from new temperature.
      theta_work(1,1,k) = t_earlier(1,1,k) / exner_wth(map_wth(1) + k)

      ! All increments found from difference between the updated *_work values
      ! and the values that were intent in.
      dtheta_inc_wth(map_wth(1) + k)  = (theta_work(1,1,k) - theta_wth(map_wth(1) + k)) &
                                            ! Above is the change as a result of the forcing
                                        - (dtdt(1,1,k) / exner_wth(map_wth(1) + k) )
                                            ! so take forcing off to get the response to it.

      dqv_inc_wth    (map_wth(1) + k) = (qv_work (1,1,k) - mv_wth(map_wth(1) + k)) &
                                            ! Above is the change as a result of the forcing
                                            - dqdt(1,1,k)
                                            ! so take forcing off to get the response to it.

      dqcl_inc_wth   (map_wth(1) + k) = (qcl_work(1,1,k) - ml_wth(map_wth(1) + k)) &
                                            ! Above is the change as a result of the forcing
                                            - dldt(1,1,k)
                                            ! so take forcing off to get the response to it.

      dcfl_inc_wth   (map_wth(1) + k) = cfl_work(1,1,k) - cfl_wth( map_wth(1) + k )

      dbcf_inc_wth   (map_wth(1) + k) = bcf_work(1,1,k) - bcf_wth( map_wth(1) + k )

      ! N.B. we are not updating the variables that came in, just
      !      providing increments that need to be added on later.

    end do

    ! Copy all level 1 info into level 0.
    dtheta_inc_wth (map_wth(1) + 0) = dtheta_inc_wth (map_wth(1) + 1)
    dqv_inc_wth    (map_wth(1) + 0) = dqv_inc_wth    (map_wth(1) + 1)
    dqcl_inc_wth   (map_wth(1) + 0) = dqcl_inc_wth   (map_wth(1) + 1)
    dcfl_inc_wth   (map_wth(1) + 0) = dcfl_inc_wth   (map_wth(1) + 1)
    dbcf_inc_wth   (map_wth(1) + 0) = dbcf_inc_wth   (map_wth(1) + 1)

end subroutine pc2_homog_iau_code

end module pc2_homog_iau_kernel_mod
