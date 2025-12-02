!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Remove SPT moisture increments if saturation is reached
module spt_saturation_cap_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD,        &
                               GH_WRITE, GH_REAL,         &
                               GH_READ, CELL_COLUMN

  use fs_continuity_mod, only: Wtheta
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: spt_saturation_cap_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & !dtheta
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & !dmv
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),  & !T_latest
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),  & !mv
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA)   & !pressure
         /)
         integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass ::  spt_saturation_cap_code
  end type spt_saturation_cap_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public  spt_saturation_cap_code
contains

  !> @brief Set Remove points where SPT perturbations reach saturation
  !> @param[in]       nlayers   The number of layers
  !> @param[in,out]   dtheta    Field with theta perturbations
  !> @param[in,out]   dmv       Field with mv perturbations
  !> @param[in]       T_latest  Latest Temperature (after solver)
  !> @param[in]       mv        Water vapour mixing ratios
  !> @param[in]       pressure  Pressure
  !> @param[in]       ndf_wth   Number of degrees of freedom per cell for wtheta
  !> @param[in]       undf_wth  Number of total degrees of freedom for wtheta
  !> @param[in]       map_wth   Dofmap for the cell at the base of the column

  subroutine  spt_saturation_cap_code(nlayers,  &
                                      dtheta,   &
                                      dmv,      &
                                      T_latest, &
                                      mv,       &
                                      pressure, &
                                      ndf_wth,  &
                                      undf_wth, &
                                      map_wth   &
                                      )

    use stochastic_physics_config_mod,    only: spt_mse_conservation, &
                                                spt_level_bottom, &
                                                spt_level_top

    use qsat_mod,                   only: qsat_wat_mix

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in),  dimension(ndf_wth)  :: map_wth

    ! Fields perturbations + tendencies
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: dtheta
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: dmv
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: T_latest
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: mv
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: pressure


    !UM definition for qsat capping
    real(kind=r_def) :: qsat

    integer(kind=i_def) :: k

    do k= spt_level_bottom, spt_level_top

      ! Call qsat from UM routine
      call qsat_wat_mix(qsat,  T_latest(map_wth(1) + k), pressure(map_wth(1) + k))

      ! ! Remove points where perturbed mv is above saturation
      if ( mv(map_wth(1) + k) + dmv(map_wth(1) + k) > qsat ) then
        if (spt_mse_conservation) then
            dtheta(map_wth(1) + k) = 0.0_r_def
        end if
        dmv(map_wth(1) + k) = 0.0_r_def
      end if
    end do

  end subroutine  spt_saturation_cap_code

end module spt_saturation_cap_kernel_mod
