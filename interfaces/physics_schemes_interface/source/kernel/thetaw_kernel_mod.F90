!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to wet bulb potential temperature code.

module thetaw_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_SCALAR,       &
                           GH_READ, GH_READWRITE,     &
                           GH_INTEGER,                &
                           GH_REAL, CELL_COLUMN,      &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def, l_def, r_um
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: thetaw_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR,GH_INTEGER, GH_READ)                               &
! PSyclone currently doesn't support passing scalar arrays, see PSyclone #1312.
! The following code can be uncommented when this is issue is resolved.
!        arg_type(GH_SCALAR_ARRAY,GH_REAL, GH_READ, 1)                         &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: thetaw_code
  end type thetaw_kernel_type

  public :: thetaw_code

contains

  !> @details Calculates wet bulb potential temperature on pressure levels
  !>           using UM routine thetaw.
  !> @param[in]     nlayers       The number of layers
  !> @param[in,out] thetaw        Wet bulb potential temperature (K)
  !> @param[in]     temperature   Temperature on pressure levels (K)
  !> @param[in]     qv            Specific humidity on pressure levels (kg/kg)
  !> @param[in]     nplev         Number of pressure levels
  !> @param[in]     plevs         Pressure level values (Pa)
  !> @param[in]     ndf           Number of degrees of freedom per cell
  !> @param[in]     undf          Number of total degrees of freedom
  !> @param[in]     map           Dofmap for the cell at the base of the column
  subroutine thetaw_code(nlayers,        &
                         theta_w,        &
                         temperature,    &
                         qv,             &
                         nplev,          &
                         plevs,          &
                         ndf, undf, map)

    use thetaw_mod, only: thetaw

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers, nplev
    integer(kind=i_def), intent(in) :: ndf, undf
    integer(kind=i_def), intent(in), dimension(ndf)  :: map

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf) :: temperature, qv
    real(kind=r_def),    intent(inout), dimension(undf) :: theta_w

    ! Constants passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(nplev) :: plevs

    ! Internal variables
    integer(kind=i_def) :: kp
    integer(kind=i_def), parameter :: segment_length = 1
    ! Logical to output potential temperature (rather than temperature) from
    ! thetaw subroutine
    logical(kind=l_def), parameter :: l_potential = .true.
    real(r_um), dimension(segment_length) :: temperature_um, qv_um, &
         pressure_um, theta_w_um

    do kp = 1, nplev

      temperature_um = real(temperature(map(1)+kp-1),r_um)
      qv_um = real(qv(map(1)+kp-1),r_um)
      pressure_um = real(plevs(kp),r_um)
      ! Call thetaw code and update thetaw array
      call thetaw(segment_length, temperature_um, qv_um, &
                  pressure_um, l_potential, theta_w_um )
      theta_w(map(1)+kp-1) = theta_w_um(1)

    end do

  end subroutine thetaw_code

end module thetaw_kernel_mod
