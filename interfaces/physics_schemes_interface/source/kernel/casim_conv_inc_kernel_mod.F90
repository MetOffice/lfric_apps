!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Distribute convection increment between cloud species
!>
module casim_conv_inc_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 GH_READWRITE,          &
                                 CELL_COLUMN

  use constants_mod,      only : r_def, i_def
  use fs_continuity_mod,  only : Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: casim_conv_inc_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                        &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),   & ! m_ci
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA)    & ! m_s
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: casim_conv_inc_code
  end type

  public :: casim_conv_inc_code

contains

  !> @details The GR convection scheme works on the total mass of ice cloud
  !>          m_ci + m_s, providing a single increment to this.
  !>          When the increment is positive, it's added to m_ci, since it
  !>          is assumed this is small ice. However, when the total increment
  !>          is negative, it needs to be removed from both species to avoid
  !>          negative values. We remove it from m_ci first, followed by m_s.
  !> @param[in]     nlayers       Number of layers
  !> @param[in,out] m_ci          Cloud ice mixing ratio in wth
  !> @param[in,out] m_s           Cloud snow mixing ratio in wth
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom for potential temperature space
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for potential temperature space
  subroutine casim_conv_inc_code(nlayers,      &
                                 m_ci,         &
                                 m_s,          &
                                 ndf_wth,      &
                                 undf_wth,     &
                                 map_wth)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_wth
    integer(kind=i_def), intent(in)     :: undf_wth

    integer(kind=i_def), intent(in),    dimension(ndf_wth)  :: map_wth

    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_ci, m_s

    ! Local variables for the kernel
    integer(i_def) :: k

    do k = 0, nlayers

      ! Increment has already been added to ice by builtin function,
      ! just need to redistribute any negatives to the snow
      ! Convection shouldn't be able to reduce the total (ice+snow) to below
      ! zero, therefore if the ice is negative there will always be enough
      ! snow to borrow from
      if (m_ci(map_wth(1)+k) < 0.0_r_def) then
        m_s(map_wth(1)+k) = m_s(map_wth(1)+k) + m_ci(map_wth(1)+k)
        m_ci(map_wth(1)+k) = 0.0_r_def
      end if

    end do

  end subroutine casim_conv_inc_code

end module casim_conv_inc_kernel_mod
