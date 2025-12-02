!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to wet bulb freezing level calculation.

module freeze_lev_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD,                  &
                                  GH_READ, GH_WRITE,         &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,        only: r_def, i_def, l_def, r_um
  use fs_continuity_mod,    only: Wtheta
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: freeze_lev_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    &!theta
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    &!mv
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    &!total_mass
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    &!exner
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    &!height
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &!freeze_lev
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: freeze_lev_code
  end type freeze_lev_kernel_type

  public :: freeze_lev_code

contains

  !> @details Calculates wet bulb freezing level using UM routine thetaw.
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     theta         Potential temperature (K)
  !> @param[in]     mv            Vapour mixing ratio (kg/kg)
  !> @param[in]     total_mass    Total of all mixing ratios (kg/kg)
  !> @param[in]     exner         Exner pressure (Pa)
  !> @param[in]     height        Height above planet_radius (m)
  !> @param[out]    freeze_lev    Wet bulb Freezing level height (m)
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell
  !> @param[in]     undf_wth      Number of total degrees of freedom
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell
  !> @param[in]     undf_2d       Number of total degrees of freedom
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column
  subroutine freeze_lev_code(nlayers,                         &
                             theta,                           &
                             mv,                              &
                             total_mass,                      &
                             exner,                           &
                             height,                          &
                             freeze_lev,                      &
                             ndf_wth, undf_wth, map_wth,      &
                             ndf_2d, undf_2d, map_2d)

    use thetaw_mod, only: thetaw
    use planet_constants_mod, only: p_zero, kappa, lapse
    use conversions_mod, only: zerodegc

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_wth) :: mv
    real(kind=r_def),    intent(in), dimension(undf_wth) :: theta
    real(kind=r_def),    intent(in), dimension(undf_wth) :: exner
    real(kind=r_def),    intent(in), dimension(undf_wth) :: total_mass
    real(kind=r_def),    intent(in), dimension(undf_wth) :: height

    real(kind=r_def),    intent(inout), dimension(undf_2d) :: freeze_lev

    ! Internal variables
    integer(kind=i_def) :: k
    integer(kind=i_def), parameter :: seg_len = 1
    logical(kind=l_def), parameter :: l_potential = .false.
    logical(kind=l_def) :: freeze_lev_found
    real(kind=r_um) :: temperature(seg_len), pressure(seg_len), qv(seg_len), &
                        temperature_wb(nlayers), dz, frac

    k = 1
    ! Save input temperature, specific humidity and pressure
    temperature = theta(map_wth(1)+k) * exner(map_wth(1)+k)
    qv = mv(map_wth(1)+k) / total_mass(map_wth(1)+k)
    pressure = p_zero*(exner(map_wth(1)+k))**(1.0_r_def/kappa)
    ! Call thetaw code and update thetaw array
    call thetaw(seg_len, temperature, qv, pressure, l_potential, &
                temperature_wb(k))
    if (temperature_wb(k) <= zerodegc) then
      dz = (temperature_wb(k) - zerodegc) / lapse
      freeze_lev(map_2d(1)) = max(height(map_wth(1)+k) + dz, 0.0_r_def)
      freeze_lev_found = .true.
    else
      freeze_lev_found = .false.
    end if

    if (.not. freeze_lev_found) then
      do k = 2, nlayers

        ! Save input temperature, specific humidity and pressure
        temperature = theta(map_wth(1)+k) * exner(map_wth(1)+k)
        qv = mv(map_wth(1)+k) / total_mass(map_wth(1)+k)
        pressure = p_zero*(exner(map_wth(1)+k))**(1.0_r_def/kappa)
        ! Call thetaw code and update thetaw array
        call thetaw(seg_len, temperature, qv, pressure, l_potential, &
                    temperature_wb(k))

        if (temperature_wb(k) <= zerodegc) then
          frac = (zerodegc - temperature_wb(k-1)) / &
                 (temperature_wb(k) - temperature_wb(k-1))
          freeze_lev(map_2d(1)) = height(map_wth(1)+k) * frac + &
                                  height(map_wth(1)+k-1) * (1.0_r_def - frac)
          exit
        end if

      end do
    end if

  end subroutine freeze_lev_code

end module freeze_lev_kernel_mod
