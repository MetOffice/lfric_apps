!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code (Claude Opus 5.5), 2026-10-20.
!-----------------------------------------------------------------------------
!> @brief Magnitude of the orographic drag stress on one pressure level.

module mountain_wave_stress_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_SCALAR,       &
                           GH_READ, GH_READWRITE,     &
                           GH_REAL, GH_INTEGER,       &
                           CELL_COLUMN,               &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def
  use fs_continuity_mod, only: Wtheta
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: mountain_wave_stress_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,                         &
                                         ANY_DISCONTINUOUS_SPACE_1),           &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: mountain_wave_stress_code
  end type mountain_wave_stress_kernel_type

  public :: mountain_wave_stress_code

contains

  !> @brief   Magnitude of the orographic drag stress on one pressure level.
  !> @details The eastward and northward stress components are interpolated
  !>          linearly in pressure between the model levels either side of
  !>          the pressure level (using the lowest or highest two levels if
  !>          the level is outside the column), and the magnitude of the
  !>          interpolated stress is returned.
  !> @param[in]     nlayers     Number of layers
  !> @param[in]     taux_orog   Eastward orographic drag stress (Pa)
  !> @param[in]     tauy_orog   Northward orographic drag stress (Pa)
  !> @param[in]     exner_wth   Exner pressure on model levels
  !> @param[in,out] stress_plev Magnitude of the stress on pressure
  !>                            levels (Pa)
  !> @param[in]     press_lev   Pressure of this level (Pa)
  !> @param[in]     kp          Index of this level in the output field
  !> @param[in]     p_zero      Reference pressure (Pa)
  !> @param[in]     kappa       Ratio of gas constant to specific heat
  !> @param[in]     ndf_wth     Number of degrees of freedom per cell for
  !>                            Wtheta
  !> @param[in]     undf_wth    Number of unique degrees of freedom for
  !>                            Wtheta
  !> @param[in]     map_wth     Dofmap for the cell at the base of the column
  !> @param[in]     ndf_pl      Number of degrees of freedom per cell for the
  !>                            pressure-level field
  !> @param[in]     undf_pl     Number of unique degrees of freedom for the
  !>                            pressure-level field
  !> @param[in]     map_pl      Dofmap for the pressure-level field
  subroutine mountain_wave_stress_code(nlayers,                       &
                                       taux_orog, tauy_orog,          &
                                       exner_wth, stress_plev,        &
                                       press_lev, kp, p_zero, kappa,  &
                                       ndf_wth, undf_wth, map_wth,    &
                                       ndf_pl, undf_pl, map_pl)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in) :: ndf_pl, undf_pl
    integer(kind=i_def), intent(in), dimension(ndf_pl) :: map_pl

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(in),    dimension(undf_wth) :: taux_orog
    real(kind=r_def), intent(in),    dimension(undf_wth) :: tauy_orog
    real(kind=r_def), intent(in),    dimension(undf_wth) :: exner_wth
    real(kind=r_def), intent(inout), dimension(undf_pl)  :: stress_plev
    real(kind=r_def), intent(in)    :: press_lev
    integer(kind=i_def), intent(in) :: kp
    real(kind=r_def), intent(in)    :: p_zero, kappa

    ! Internal variables
    integer(kind=i_def) :: k, k_lower, dof_lower, dof_upper
    real(kind=r_def)    :: p_lower, p_upper, weight, taux, tauy

    ! Model level k_lower is the level below the pressure level. Levels
    ! 1 to nlayers are used; level 0 is the surface.
    if ( pressure(exner_wth(map_wth(1) + 1)) <= press_lev ) then
      k_lower = 1_i_def
    else if ( pressure(exner_wth(map_wth(1) + nlayers)) > press_lev ) then
      k_lower = nlayers - 1_i_def
    else
      k_lower = 1_i_def
      do k = 1, nlayers - 1
        if ( pressure(exner_wth(map_wth(1) + k + 1)) <= press_lev ) then
          k_lower = k
          exit
        end if
      end do
    end if

    dof_lower = map_wth(1) + k_lower
    dof_upper = dof_lower + 1
    p_lower = pressure(exner_wth(dof_lower))
    p_upper = pressure(exner_wth(dof_upper))

    weight = ( press_lev - p_lower ) / ( p_upper - p_lower )
    taux = taux_orog(dof_lower)                                             &
           + weight * ( taux_orog(dof_upper) - taux_orog(dof_lower) )
    tauy = tauy_orog(dof_lower)                                             &
           + weight * ( tauy_orog(dof_upper) - tauy_orog(dof_lower) )

    stress_plev(map_pl(1) + kp - 1) = sqrt( taux**2 + tauy**2 )

  contains

    !> @brief  Pressure from Exner pressure.
    !> @param[in] exner Exner pressure
    !> @return    Pressure (Pa)
    function pressure(exner) result(p)

      implicit none

      real(kind=r_def), intent(in) :: exner
      real(kind=r_def) :: p

      p = p_zero * exner**( 1.0_r_def / kappa )

    end function pressure

  end subroutine mountain_wave_stress_code

end module mountain_wave_stress_kernel_mod
