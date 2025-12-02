!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief q_tests

module q_limit_test_kernel_mod

  use argument_mod,      only : arg_type,     &
                                GH_FIELD,     &
                                GH_REAL,      &
                                GH_READWRITE, &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : WTHETA
  use kernel_mod,        only : kernel_type

  implicit none

  private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

  type, public, extends(kernel_type) :: q_limit_test_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                  &
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! q_wth
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! qcl_wth
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! qcf_wth
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! area_frac_wth
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! bulk_frac_wth
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! liquid_frac_wth
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA)  & ! frozen_frac_wth
       /)
      integer :: operates_on = CELL_COLUMN
    contains
      procedure, nopass :: q_limit_test_code
    end type

!-------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------
    public :: q_limit_test_code

contains

!> @brief Interface
!> @param[in]        nlayers             number of layers
!> @param[in,out]    q_wth               vapour specific humidity after time step
!> @param[in,out]    qcl_wth             liquid specific humidity after time step
!> @param[in,out]    qcf_wth             frozen specific humidity after time step
!> @param[in,out]    area_frac_wth       area cloud fraction
!> @param[in,out]    bulk_frac_wth       bulk cloud fraction
!> @param[in,out]    liquid_frac_wth     liquid cloud fraction
!> @param[in,out]    frozen_frac_wth     frozen cloud fraction
!> @param[in]        ndf_wth             the number of degrees of freedom for a field in wtheta space.
!> @param[in]        undf_wth            the number of unique degrees of freedom for a field in wtheta
!> @param[in]        map_wth             dofmap for the cell at the base of the column
!!


subroutine q_limit_test_code( nlayers,                           &
                              q_wth,qcl_wth,qcf_wth ,            &
                              area_frac_wth, bulk_frac_wth,      &
                              frozen_frac_wth, liquid_frac_wth , &
                              ndf_wth, undf_wth, map_wth         &
                             )

      implicit none

      integer(kind=i_def), intent(in) :: nlayers
      integer(kind=i_def), intent(in) :: ndf_wth
      integer(kind=i_def), intent(in) :: undf_wth

      real(kind=r_def), intent(inout),  dimension(undf_wth) :: q_wth
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: qcl_wth
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: qcf_wth
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: area_frac_wth
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: bulk_frac_wth
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: frozen_frac_wth
      real(kind=r_def), intent(inout),  dimension(undf_wth) :: liquid_frac_wth
      integer(kind=i_def), intent(in),  dimension(ndf_wth)  :: map_wth

      integer(kind=i_def) :: k
      real(kind=r_def), parameter :: zero = 0.0_r_def
      real(kind=r_def), parameter :: qmin = 1.0e-08_r_def
      real(kind=r_def), parameter :: tol  = 0.0_r_def

      do k= 0, nlayers

          if ( q_wth(map_wth(1) + k) < qmin )  q_wth(map_wth(1) + k) = qmin

          if ( qcl_wth(map_wth(1) + k) <= zero ) then
            qcl_wth(map_wth(1) + k) = zero
            liquid_frac_wth(map_wth(1) + k) = zero
            if ( qcf_wth(map_wth(1)+k) <= tol ) then
              area_frac_wth(map_wth(1) + k) = zero
              bulk_frac_wth(map_wth(1) + k) = zero
            end if
          end if

          if ( qcf_wth(map_wth(1) + k) <= zero )  then
            qcf_wth(map_wth(1) + k) = zero
            frozen_frac_wth(map_wth(1) + k) = zero
            if ( qcl_wth(map_wth(1) + k) <= tol ) then
              area_frac_wth(map_wth(1) + k) = zero
              bulk_frac_wth(map_wth(1) + k) = zero
            end if
          end if

          if ( ( liquid_frac_wth(map_wth(1) + k) == zero ) .and. &
               ( frozen_frac_wth(map_wth(1) + k) == zero ) ) then
            area_frac_wth(map_wth(1) + k) = zero
            bulk_frac_wth(map_wth(1) + k) = zero
          end if

      end do

end subroutine q_limit_test_code

end module q_limit_test_kernel_mod
