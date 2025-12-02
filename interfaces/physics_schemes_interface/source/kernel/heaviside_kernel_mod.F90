!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to Heaviside function calculation.

module heaviside_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_READWRITE,     &
                                  GH_INTEGER,                &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: heaviside_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),     &
         arg_type(GH_SCALAR,GH_INTEGER, GH_READ),                              &
!        arg_type(GH_SCALAR_ARRAY,GH_REAL, GH_READ, 1), see PSyclone issue #1312
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                                 &
         arg_type(GH_SCALAR,GH_REAL, GH_READ)                                  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: heaviside_code
  end type heaviside_kernel_type

  public :: heaviside_code

contains

  !> @details Calculates the Heaviside function, i.e. 1 if the pressure level
  !!          is above the surface and 0 if the pressure level is below it.
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     exner         Exner pressure at theta levels
  !> @param[in]     nplev         Number of pressure levels
  !> @param[in]     plevs         Pressure level values
  !> @param[in,out] plev_data     Data field on pressure levels
  !> @param[in]     p_zero        Reference pressure
  !> @param[in]     kappa         Rd / cp
  !> @param[in]     ndf_in        Number of degrees of freedom per cell for in fields
  !> @param[in]     undf_in       Number of total degrees of freedom for in fields
  !> @param[in]     map_in        Dofmap for the cell at the base of the column for in fields
  !> @param[in]     ndf_out       Number of degrees of freedom per cell for out fields
  !> @param[in]     undf_out      Number of total degrees of freedom for out fields
  !> @param[in]     map_out       Dofmap for the cell at the base of the column for out fields
  subroutine heaviside_code(nlayers,                    &
                            exner,                      &
                            nplev,                      &
                            plevs,                      &
                            plev_data,                  &
                            p_zero,                     &
                            kappa,                      &
                            ndf_in, undf_in, map_in,    &
                            ndf_out, undf_out, map_out)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers, nplev
    integer(kind=i_def), intent(in) :: ndf_out, undf_out
    integer(kind=i_def), intent(in), dimension(ndf_out)  :: map_out
    integer(kind=i_def), intent(in) :: ndf_in, undf_in
    integer(kind=i_def), intent(in), dimension(ndf_in)  :: map_in

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_in) :: exner
    real(kind=r_def),    intent(inout), dimension(undf_out) :: plev_data

    ! Constants passed explicitly from algorithm
    real(kind=r_def),    intent(in) :: p_zero, kappa
    real(kind=r_def),    intent(in), dimension(nplev) :: plevs

    ! Internal variables
    integer(kind=i_def) :: kp
    real(kind=r_def) :: desired_ex

    do kp = 1, nplev

      ! Level we want
      desired_ex = (plevs(kp)/p_zero)**kappa

      if ((exner(map_in(1)) < desired_ex) .or. &
          (exner(map_in(1)+nlayers) > desired_ex)) then

        plev_data(map_out(1)+kp-1) = 0.0_r_def

      end if

    end do

  end subroutine heaviside_code

end module heaviside_kernel_mod
