!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief q_tests

module tropopause_limits_kernel_mod

  use argument_mod,      only : arg_type,     &
                                GH_FIELD,     &
                                GH_REAL,      &
                                GH_INTEGER,   &
                                GH_READWRITE, &
                                GH_READ,      &
                                CELL_COLUMN,  &
                                ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : WTHETA
  use kernel_mod,        only : kernel_type
  use qsat_mod,          only : qsat_mix

  implicit none

  private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

  type, public, extends(kernel_type) :: tropopause_limits_kernel_type
    private
    type(arg_type) :: meta_args(14) = (/                  &
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! m_v
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! m_cl
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! m_ci
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! m_r
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! m_s
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! m_g
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! liquid_fraction
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! frozen_fraction
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! area_fraction
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! bulk_fraction
        arg_type(GH_FIELD,  GH_REAL, GH_READ,       WTHETA), & ! mv_saved
        arg_type(GH_FIELD,  GH_REAL, GH_READ,       WTHETA), & ! temperature
        arg_type(GH_FIELD,  GH_REAL, GH_READ,       WTHETA), & ! pressure
        arg_type(GH_FIELD, GH_INTEGER,GH_READ,ANY_DISCONTINUOUS_SPACE_1) & ! trop_level
       /)
      integer :: operates_on = CELL_COLUMN
    contains
      procedure, nopass :: tropopause_limits_kernel_code
    end type

!-------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------
    public :: tropopause_limits_kernel_code

contains

!> @brief Interface to the pc2 homogeneous forcing code
!> @param[in]        nlayers             number of layers
!> @param[in,out]    m_v                 vapour mixing ratio after iau increment
!> @param[in,out]    m_cl                liquid mixing ratio after iau increment
!> @param[in,out]    m_ci                frozen mixing ratio after iau increment
!> @param[in,out]    m_r                 rain mixing ratio after iau increment
!> @param[in,out]    m_s                 snow mixing ratio after iau increment
!> @param[in,out]    m_g                 graupel mixing ratio after iau increment
!> @param[in]        temperature         temperature on theta levels after iau inc
!> @param[in]        pressure            pressure on theta levels after iau inc
!> @param[in]        ndf_wth             the number of degrees of freedom for a field in wtheta space.
!> @param[in]        undf_wth            the number of unique degrees of freedom for a field in wtheta
!> @param[in]        map_wth             dofmap for the cell at the base of the column

  subroutine tropopause_limits_kernel_code( nlayers,                           &
                                            m_v, m_cl, m_ci, m_r, m_s, m_g,    &
                                            liquid_fraction,                   &
                                            frozen_fraction,                   &
                                            area_fraction,                     &
                                            bulk_fraction,                     &
                                            mv_saved,                          &
                                            temperature,                       &
                                            pressure, trop_level,              &
                                            ndf_wth, undf_wth, map_wth,        &
                                            ndf_2d, undf_2d, map_2d )

      implicit none

      ! Arguments
      integer(kind=i_def), intent(in) :: nlayers
      integer(kind=i_def), intent(in) :: ndf_wth, ndf_2d
      integer(kind=i_def), intent(in) :: undf_wth, undf_2d

      real(kind=r_def), intent(inout),  dimension(undf_wth) :: m_v
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: mv_saved
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: m_cl
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: m_ci
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: m_r
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: m_s
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: m_g
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: liquid_fraction
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: frozen_fraction
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: area_fraction
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: bulk_fraction
      real(kind=r_def), intent(in),     dimension(undf_wth) :: temperature
      real(kind=r_def), intent(in),     dimension(undf_wth) :: pressure
      integer(kind=i_def), intent(in),  dimension(undf_2d)  :: trop_level
      integer(kind=i_def), intent(in),  dimension(ndf_wth)  :: map_wth
      integer(kind=i_def), intent(in),  dimension(ndf_2d)   :: map_2d

      integer(kind=i_def) :: k
      integer(kind=i_def) :: i_trop
      real(kind=r_def), parameter :: trop_min_rh = 0.01_r_def

      !UM definition for qsat capping
      real(kind=r_def) :: msat_out

      ! Level of tropopause
      i_trop = trop_level(map_2d(1))

      do k = 0, nlayers

        ! Apply limits depending on whether in troposphere
        if ( k <= i_trop ) then
        ! in troposphere

          ! Call qsat_mix from UM routine
          call qsat_mix( msat_out, temperature(map_wth(1) + k), pressure(map_wth(1) + k) )

          m_v(map_wth(1) + k) = MAX( m_v(map_wth(1) + k),msat_out * trop_min_rh )

        else
        ! above troposphere, so remove cloud and limit q_v

          m_v(map_wth(1) + k)  = mv_saved(map_wth(1) + k)
          m_cl(map_wth(1) + k) = 0.0_r_def
          m_ci(map_wth(1) + k) = 0.0_r_def
          m_r(map_wth(1) + k)  = 0.0_r_def
          m_s(map_wth(1) + k)  = 0.0_r_def
          m_g(map_wth(1) + k)  = 0.0_r_def
          liquid_fraction(map_wth(1) + k) = 0.0_r_def
          frozen_fraction(map_wth(1) + k) = 0.0_r_def
          area_fraction(map_wth(1) + k)   = 0.0_r_def
          bulk_fraction(map_wth(1) + k)   = 0.0_r_def

        end if

      end do

  end subroutine tropopause_limits_kernel_code

end module tropopause_limits_kernel_mod
