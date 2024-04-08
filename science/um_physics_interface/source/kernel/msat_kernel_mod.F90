!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Return the saturation mixing ratio from the UM function

module msat_kernel_mod

  use argument_mod,         only: arg_type,            &
                                  GH_FIELD, GH_SCALAR, &
                                  GH_READ, GH_WRITE, GH_LOGICAL, &
                                  GH_REAL, CELL_COLUMN, &
                                  ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,        only: r_def, i_def, l_def
  use kernel_mod,           only: kernel_type
  use planet_constants_mod, only: recip_kappa, p_zero

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: msat_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                     &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                           &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: msat_code
  end type msat_kernel_type

  public :: msat_code

contains

  !> @details Call the UM msat function, with logical detemining whether the
  !>          water or ice version is required, and return the saturation
  !>          mixing ratio
  !> @param[in]     nlayers     The number of layers
  !> @param[out]    msat        saturation mixing ratio
  !> @param[in]     exner       exner pressure in theta space
  !> @param[in]     theta       potential temperature
  !> @param[in]     water_only  Water or ice below 0C
  !> @param[in]     ndf         Number of degrees of freedom per cell
  !> @param[in]     undf        Number of total degrees of freedom
  !> @param[in]     map         Dofmap for the cell at the base of the column

  subroutine msat_code(nlayers,            &
                       msat,               &
                       exner,              &
                       theta,              &
                       water_only,         &
                       ndf,                &
                       undf,               &
                       map)

    use qsat_mod, only: qsat_mix, qsat_wat_mix

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf
    integer(kind=i_def), intent(in), dimension(ndf)  :: map

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf) :: exner
    real(kind=r_def),    intent(in), dimension(undf) :: theta
    logical(kind=l_def), intent(in) :: water_only
    real(kind=r_def),    intent(inout), dimension(undf) :: msat

    ! Internal variables
    real(kind=r_def) :: pressure(0:nlayers), temperature(0:nlayers)
    integer(kind=i_def) :: k

    ! Calculate pressure and temperature
    do k = 0, nlayers
      pressure(k) = p_zero * (exner(map(1)+k)) ** recip_kappa
      temperature(k) = theta(map(1)+k) * exner(map(1)+k)
    end do

    ! Calculate msat
    if (water_only) then
      do k = 0, nlayers
        call qsat_wat_mix(msat(map(1)+k),temperature(k),pressure(k))
      end do
    else
      do k = 0, nlayers
        call qsat_mix(msat(map(1)+k),temperature(k),pressure(k))
      end do
    end if

  end subroutine msat_code

end module msat_kernel_mod
