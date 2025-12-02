!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief q_tests

module saturation_limit_kernel_mod

  use argument_mod,      only : arg_type,     &
                                GH_FIELD,     &
                                GH_REAL,      &
                                GH_READWRITE, &
                                GH_READ,      &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : WTHETA
  use kernel_mod,        only : kernel_type
  use qsat_mod,          only : qsat_wat_mix

  implicit none

  private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

  type, public, extends(kernel_type) :: saturation_limit_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                  &
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! q_v
        arg_type(GH_FIELD,  GH_REAL, GH_READ,       WTHETA), & ! temperature
        arg_type(GH_FIELD,  GH_REAL, GH_READ,       WTHETA)  & ! pressure
       /)
      integer :: operates_on = CELL_COLUMN
    contains
      procedure, nopass :: saturation_limit_kernel_code
    end type

!-------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------
    public :: saturation_limit_kernel_code

contains

!> @brief Interface to the pc2 homogeneous forcing code
!> @param[in]        nlayers             number of layers
!> @param[in,out]    m_v                 vapour mixing ratio after iau increment
!> @param[in]        temperature         temperature on theta levels after iau inc
!> @param[in]        pressure            pressure on theta levels after iau inc
!> @param[in]        ndf_wth             the number of degrees of freedom for a field in wtheta space.
!> @param[in]        undf_wth            the number of unique degrees of freedom for a field in wtheta
!> @param[in]        map_wth             dofmap for the cell at the base of the column

  subroutine saturation_limit_kernel_code( nlayers,                    &
                                           m_v, temperature, pressure, &
                                           ndf_wth, undf_wth, map_wth  &
                                          )

      implicit none

      ! Arguments
      integer( kind=i_def ), intent(in) :: nlayers
      integer( kind=i_def ), intent(in) :: ndf_wth
      integer( kind=i_def ), intent(in) :: undf_wth

      real( kind=r_def ), intent(inout),  dimension(undf_wth) :: m_v
      real( kind=r_def ), intent(in),     dimension(undf_wth) :: temperature
      real( kind=r_def ), intent(in),     dimension(undf_wth) :: pressure
      integer( kind=i_def ), intent(in),  dimension(ndf_wth)  :: map_wth

      integer( kind=i_def ) :: k

      !UM definition for qsat capping
      real( kind=r_def ) :: msat

      do k = 0, nlayers

        ! Call qsat_wat_mix from UM routine
        call qsat_wat_mix( msat, temperature(map_wth(1) + k), pressure(map_wth(1) + k) )

        if ( m_v( map_wth(1) + k ) > msat ) then

          m_v( map_wth(1) + k ) = msat

        end if

      end do

  end subroutine saturation_limit_kernel_code

end module saturation_limit_kernel_mod
