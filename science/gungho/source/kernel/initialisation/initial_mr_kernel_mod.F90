!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial mixing ratio fields

!> @details The kernel computes initial mixing ratio fields for mr in the same
!>          space as that of theta

module initial_mr_kernel_mod

    use argument_mod,                  only: arg_type, func_type,       &
                                             GH_FIELD, GH_REAL,         &
                                             GH_SCALAR, GH_BASIS,       &
                                             GH_WRITE, GH_READ,         &
                                             ANY_SPACE_9,               &
                                             ANY_DISCONTINUOUS_SPACE_3, &
                                             GH_EVALUATOR, CELL_COLUMN
    use fs_continuity_mod,             only: W3, Wtheta
    use constants_mod,                 only: r_def, i_def
    use kernel_mod,                    only: kernel_type
    use section_choice_config_mod,     only: cloud, cloud_um
    use idealised_config_mod,          only: test, test_bryan_fritsch, &
                                             test_grabowski_clark, &
                                             test_isot_dry_atm
    use initial_pressure_config_mod,   only: method, method_balanced

    implicit none
    private

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_mr_kernel_type
        private
        type(arg_type) :: meta_args(9) = (/                                      &
             arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                    &
             arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        &
             arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        &
             arg_type(GH_FIELD*6, GH_REAL, GH_WRITE, Wtheta),                    &
             arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
             arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
             arg_type(GH_SCALAR,  GH_REAL, GH_READ),                             &
             arg_type(GH_SCALAR,  GH_REAL, GH_READ),                             &
             arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
             /)
        type(func_type) :: meta_funcs(1) = (/                                    &
             func_type(ANY_SPACE_9, GH_BASIS)                                    &
             /)
        integer :: operates_on = CELL_COLUMN
        integer :: gh_shape = GH_EVALUATOR
    contains
        procedure, nopass :: initial_mr_code
    end type

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public :: initial_mr_code
contains

  !> @brief The subroutine which is called directly by the Psy layer
  !! @param[in] nlayers Integer the number of layers
  !! @param[in] theta Potential temperature
  !! @param[in] exner Exner pressure variable
  !! @param[in] rho Density of dry air
  !! @param[in,out] mr_v Water vapour mixing ratio
  !! @param[in,out] mr_cl Liquid cloud mixing ratio
  !! @param[in,out] mr_r Rain mixing ratio
  !! @param[in,out] mr_ci Ice cloud mixing ratio
  !! @param[in,out] mr_s Snow mixing ratio
  !! @param[in,out] mr_g Graupel mixing ratio
  !! @param[in] chi_1 First component of the chi coordinate field
  !! @param[in] chi_2 Second component of the chi coordinate field
  !! @param[in] chi_3 Third component of the chi coordinate field
  !! @param[in] panel_id A field giving the ID for mesh panels
  !! @param[in] p_zero Reference surface pressure
  !! @param[in] Rd Gas constant for dry air
  !! @param[in] kappa Ratio of Rd and cp
  !! @param[in] ndf_wtheta The number of degrees of freedom per cell for wtheta
  !! @param[in] undf_wtheta The number of total degrees of freedom for wtheta
  !! @param[in] map_wtheta Integer array holding the dofmap for the cell at the base of the column
  !! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
  !! @param[in] undf_w3 Number of unique degrees of freedom  for w3
  !! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
  !! @param[in] ndf_chi Number of degrees of freedom per cell for chi
  !! @param[in] undf_chi Number of total degrees of freedom for chi
  !! @param[in] map_chi Dofmap for the cell at the base of the column
  !! @param[in] chi_basis Basis functions evaluated at Wtheta points
  !! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
  !! @param[in] undf_pid Number of unique degrees of freedom for panel_id
  !! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
  subroutine initial_mr_code(nlayers, theta, exner, rho,            &
                             mr_v, mr_cl, mr_r, mr_ci, mr_s, mr_g,  &
                             chi_1, chi_2, chi_3,                   &
                             panel_id,                              &
                             p_zero, Rd, kappa,                     &
                             ndf_wtheta, undf_wtheta, map_wtheta,   &
                             ndf_w3, undf_w3, map_w3,               &
                             ndf_chi, undf_chi, map_chi, chi_basis, &
                             ndf_pid, undf_pid, map_pid)

    use analytic_moisture_profiles_mod, only : analytic_moisture
    use sci_chi_transform_mod,          only : chi2xyz

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, undf_wtheta
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
    integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
    integer(kind=i_def), dimension(ndf_wtheta),  intent(in) :: map_wtheta
    integer(kind=i_def), dimension(ndf_w3),      intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_chi),     intent(in) :: map_chi
    integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid

    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_v, mr_cl, mr_r
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_ci, mr_s, mr_g
    real(kind=r_def), dimension(undf_wtheta), intent(in)    :: theta
    real(kind=r_def), dimension(undf_w3),     intent(in)    :: exner
    real(kind=r_def), dimension(undf_w3),     intent(in)    :: rho
    real(kind=r_def), dimension(undf_chi),    intent(in)    :: chi_1, chi_2, chi_3
    real(kind=r_def), dimension(undf_pid),    intent(in)    :: panel_id
    real(kind=r_def),                         intent(in)    :: p_zero
    real(kind=r_def),                         intent(in)    :: Rd
    real(kind=r_def),                         intent(in)    :: kappa

    real(kind=r_def), dimension(1,ndf_chi,ndf_wtheta), intent(in) :: chi_basis

    ! Internal variables
    integer(kind=i_def) :: k, df, dfc, kp1, ipanel
    real(kind=r_def)    :: theta_at_dof, rho_at_dof, pressure_at_dof
    real(kind=r_def)    :: exner_at_dof, temperature_at_dof
    real(kind=r_def)    :: coords(3), xyz(3)
    real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e

    ipanel = int(panel_id(map_pid(1)), i_def)

    ! Two separate methods for pointwise computing mixing ratio fields
    ! Either we get the pressure from Exner or from rho
    !--------------------------------------------------------------------------!
    ! Method 1: Get pressure from Exner
    !--------------------------------------------------------------------------!
    if (method == method_balanced) then

      do k = 0, nlayers - 1
        !----------------------------------------------------------------------!
        ! Find coordinates at top DoF of cell
        !----------------------------------------------------------------------!
        ! Only visit top Wtheta dof
        df = 2

        do dfc = 1, ndf_chi
          chi_1_e(dfc) = chi_1( map_chi(dfc) + k)
          chi_2_e(dfc) = chi_2( map_chi(dfc) + k)
          chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
        end do

        coords(:) = 0.0_r_def
        do dfc = 1, ndf_chi
          coords(1) = coords(1) + chi_1_e(dfc)*chi_basis(1,dfc,df)
          coords(2) = coords(2) + chi_2_e(dfc)*chi_basis(1,dfc,df)
          coords(3) = coords(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
        end do

        call chi2xyz(coords(1), coords(2), coords(3), &
                     ipanel, xyz(1), xyz(2), xyz(3))

        !----------------------------------------------------------------------!
        ! Get thermodynamic variables at DoF
        !----------------------------------------------------------------------!
        ! Extrapolate exner if at top boundary
        if (k == nlayers - 1) then
          exner_at_dof = exner(map_w3(1) + k) * sqrt(exner(map_w3(1) + k) /    &
                                                     exner(map_w3(1) + k - 1))
        else
          exner_at_dof = 0.5_r_def * (exner(map_w3(1) + k) + exner(map_w3(1) + k + 1))
        end if

        theta_at_dof = theta(map_wtheta(df) + k)
        temperature_at_dof = theta_at_dof * exner_at_dof
        pressure_at_dof = p_zero * exner_at_dof ** (1.0_r_def/kappa)

        !----------------------------------------------------------------------!
        ! Set mixing ratios at DoF
        !----------------------------------------------------------------------!
        mr_cl(map_wtheta(df) + k) = 0.0_r_def
        mr_r(map_wtheta(df) + k) = 0.0_r_def
        mr_ci(map_wtheta(df) + k) = 0.0_r_def
        mr_s(map_wtheta(df) + k) = 0.0_r_def
        mr_g(map_wtheta(df) + k) = 0.0_r_def
        mr_v(map_wtheta(df) + k) =  &
          analytic_moisture(xyz, temperature_at_dof, pressure_at_dof, test)

      end do

      ! TODO: This should either be removed or moved to end to deal with both
      ! methods. This will be done by #2877
      ! Reduce humidity at top of model for cloudy cases.
      if ( (cloud == cloud_um) .and. &
           (test /= test_bryan_fritsch) .and. (test /= test_grabowski_clark) &
           .and. ( test /= test_isot_dry_atm ) ) then
        mr_v(map_wtheta(1) + nlayers) = 1.0e-8_r_def
        mr_v(map_wtheta(1) + nlayers-1) = 1.0e-8_r_def
      end if

    !--------------------------------------------------------------------------!
    ! Method 2: Get pressure from equation of state
    !--------------------------------------------------------------------------!
    else
      do k = 0, nlayers-1
        !----------------------------------------------------------------------!
        ! Find coordinates at top DoF of cell
        !----------------------------------------------------------------------!
        ! Only visit top Wtheta dof
        df = 2

        do dfc = 1, ndf_chi
          chi_1_e(dfc) = chi_1( map_chi(dfc) + k)
          chi_2_e(dfc) = chi_2( map_chi(dfc) + k)
          chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
        end do

        coords(:) = 0.0_r_def
        do dfc = 1, ndf_chi
          coords(1) = coords(1) + chi_1_e(dfc)*chi_basis(1,dfc,df)
          coords(2) = coords(2) + chi_2_e(dfc)*chi_basis(1,dfc,df)
          coords(3) = coords(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
        end do

        call chi2xyz(coords(1), coords(2), coords(3), &
                     ipanel, xyz(1), xyz(2), xyz(3))

        !----------------------------------------------------------------------!
        ! Get thermodynamic variables at DoF
        !----------------------------------------------------------------------!

        kp1 = min(k+1,nlayers-1)
        theta_at_dof = theta(map_wtheta(df) + k)
        rho_at_dof = 0.5_r_def*(rho(map_w3(1) + k) + rho(map_w3(1) + kp1))
        pressure_at_dof = p_zero * &
           (rho_at_dof*Rd/p_zero*theta_at_dof)**(1.0_r_def/(1.0_r_def-kappa))
        exner_at_dof = (pressure_at_dof / p_zero) ** kappa
        temperature_at_dof = theta_at_dof * exner_at_dof

        !----------------------------------------------------------------------!
        ! Set mixing ratios at DoF
        !----------------------------------------------------------------------!
        mr_cl(map_wtheta(df) + k) = 0.0_r_def
        mr_r(map_wtheta(df) + k) = 0.0_r_def
        mr_ci(map_wtheta(df) + k) = 0.0_r_def
        mr_s(map_wtheta(df) + k) = 0.0_r_def
        mr_g(map_wtheta(df) + k) = 0.0_r_def
        mr_v(map_wtheta(df) + k) =  &
          analytic_moisture(xyz, temperature_at_dof, pressure_at_dof, test)
      end do

    end if

    ! Set bottom value
    k = 0
    df = 1
    mr_v(map_wtheta(df) + k) = mr_v(map_wtheta(df) + k + 1)
    mr_cl(map_wtheta(df) + k) = 0.0_r_def
    mr_r(map_wtheta(df) + k) = 0.0_r_def
    mr_ci(map_wtheta(df) + k) = 0.0_r_def
    mr_s(map_wtheta(df) + k) = 0.0_r_def
    mr_g(map_wtheta(df) + k) = 0.0_r_def

  end subroutine initial_mr_code

end module initial_mr_kernel_mod
