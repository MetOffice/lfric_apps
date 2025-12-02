!-------------------------------------------------------------------------------
!(c) Crown copyright 2020 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Holds methane oxidation code
module methox_kernel_mod

  use argument_mod,      only: arg_type,            &
                               GH_FIELD, GH_SCALAR, &
                               GH_READ, GH_WRITE,   &
                               GH_REAL, CELL_COLUMN
  use fs_continuity_mod, only: Wtheta
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  ! Maximum of 2CH4+H2O used to imply methane amount
  real(kind=r_def), parameter :: max_2ch4_h2o = 3.75e-6_r_def

  ! Coefficients used in water vapour change calculation
  real(kind=r_def), allocatable, save :: methox_coeff(:)
  real(kind=r_def), allocatable, save :: photol_coeff(:)

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: methox_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ         )  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: methox_code
  end type methox_kernel_type

  public :: methox_code
  public :: methox_coeff, photol_coeff

contains

  !> @brief Calculate chemical changes to water vapour due to methane oxidation
  !>        and hydrogen photolysis, using the method used at ECMWF
  !> @details Follows (Untch et al (ECMWF Newsletter No 82, winter
  !>          1998/99, pp 2-8) and Simmons (pers. comm.)). The model methane
  !>          mixing ratio is implicit, and derived from the assumption that
  !>          2 [CH4] + [H2O] = 3.75  ppmm throughout the stratosphere.
  !>          The methane oxidation and hydrogen photolysis rate coefficients
  !>          vary only with pressure, which is calculated for a standard
  !>          atmosphere assuming a surface pressure of p_zero.
  !> @param[in]     nlayers     The number of layers
  !> @param[in,out] dmv_methox  m_v increment from methane oxidation
  !> @param[in]     m_v         vapour mixing ratio
  !> @param[in]     dt          Timestep length (s)
  !> @param[in]     ndf_wth     Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth    Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth     Dofmap for the cell at the base of the column
  subroutine methox_code(nlayers,       &
                         dmv_methox,    &
                         m_v,           &
                         dt,            &
                         ndf_wth,       &
                         undf_wth,      &
                         map_wth)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    real(kind=r_def),    intent(in) :: dt
    real(kind=r_def),    intent(in), dimension(undf_wth) :: m_v
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: dmv_methox

    ! Internal variables
    integer(kind=i_def) :: k

    do k = 0, nlayers

      ! Only do anything if mixing ratio in the required range
      if ( m_v(map_wth(1)+k) > 0.0_r_def .and. &
           m_v(map_wth(1)+k) < max_2ch4_h2o ) then

        ! Calculate vapour increment
        dmv_methox(map_wth(1)+k) = ( methox_coeff(k) *                     &
                                     (max_2ch4_h2o - m_v(map_wth(1)+k))    &
                                   - photol_coeff(k) * m_v(map_wth(1)+k) ) &
                                   * dt

      else
        dmv_methox(map_wth(1)+k) = 0.0_r_def
      end if

    end do

  end subroutine methox_code

end module methox_kernel_mod
