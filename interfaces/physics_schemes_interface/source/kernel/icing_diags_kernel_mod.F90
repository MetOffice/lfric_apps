!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to cloud diagnostics.
!>
module icing_diags_kernel_mod

  use argument_mod,       only : arg_type,                                     &
                                 GH_FIELD, GH_REAL,                            &
                                 GH_READ, GH_READWRITE, CELL_COLUMN
  use constants_mod,      only : r_def, r_double, i_def, i_um, r_um, l_def
  use empty_data_mod,     only : empty_real_data
  use fs_continuity_mod,  only : Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: icing_diags_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! sfip
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! theta
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! exner
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! cld_fraction
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! condensate
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! w_wind
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA)  & ! height_wth
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: icing_diags_code
  end type

  public :: icing_diags_code

contains

  !> @brief Interface to the cloud diagnostics.
  !> @details Calculation of various cloud diagnostics.
  !>
  !> @param[in]     nlayers      Number of layers
  !> @param[in,out] sfip         Simplified Forecast Icing Potential
  !> @param[in]     theta        Potential temperature
  !> @param[in]     exner        Exner function
  !> @param[in]     cld_fraction Cloud fraction
  !> @param[in]     condensate   Condensate amount
  !> @param[in]     w_wind       Vertical component of wind (positive up)
  !> @param[in]     height_wth   Height above sea level in wth
  !> @param[in]     ndf_wth      Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth     Number unique of degrees of freedom for potential temperature space
  !> @param[in]     map_wth      Dofmap for the cell at the base of the column for potential temperature space

  subroutine icing_diags_code( nlayers,          &
                               sfip_wth,         &
                               theta_wth,        &
                               exner_wth,        &
                               cld_fraction_wth, &
                               condensate_wth,   &
                               w_in_wth,         &
                               height_wth,       &
                               ndf_wth,          &
                               undf_wth,         &
                               map_wth           )

    use science_conversions_mod, only: zero_degrees_celsius

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    real(kind=r_def), pointer, intent(inout)            :: sfip_wth(:)

    real(kind=r_def), intent(in), dimension(undf_wth)   :: theta_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: exner_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: cld_fraction_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: condensate_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: w_in_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: height_wth

    integer(kind=i_def) :: k
    real(kind=r_def)    :: temp, lwc, w_wind, cld_fraction, sfip

    ! Membership functions for aviation icing based on:
    real(kind=r_def) :: member_t,  & ! temperature
                        member_w,  & ! vertical velocity
                        member_clw   ! water content

    ! Thresholds for membership functions
    real(kind=r_def), parameter :: t1 = zero_degrees_celsius-28.0_r_def
    real(kind=r_def), parameter :: t2 = zero_degrees_celsius-12.0_r_def
    real(kind=r_def), parameter :: t3 = zero_degrees_celsius -1.0_r_def
    real(kind=r_def), parameter :: t4 = zero_degrees_celsius +1.0_r_def
    real(kind=r_def), parameter :: w1 =   -0.1_r_def
    real(kind=r_def), parameter :: w2 =    0.0_r_def
    real(kind=r_def), parameter :: w3 =    0.05_r_def
    real(kind=r_def), parameter :: q1 =    0.4_r_def
    real(kind=r_def), parameter :: tune_cw=0.2_r_def  ! Tuning factor for vert vel membership fn
    real(kind=r_def), parameter :: tune_cc=0.45_r_def ! Tuning factor for cloud water membership fn
    real(kind=r_def), parameter :: tune_cr=0.35_r_def ! Tuning factor for RH membership fn

    if (.not. associated(sfip_wth, empty_real_data) ) then

      do k = 1, nlayers
        ! Get hold of what we need from LFric space
        temp         = theta_wth(map_wth(1) + k) * &
                       exner_wth(map_wth(1) + k)
        lwc          = condensate_wth(map_wth(1) + k)
        w_wind       = w_in_wth(map_wth(1) + k)
        cld_fraction = cld_fraction_wth(map_wth(1) + k)
        !
        ! Simplified Forecast Icing Potential (SFIP) based on
        ! Morcrette et al (2019) https://doi.org/10.1175/WAF-D-18-0177.1
        ! which is heavily based on
        ! Bernstein et al (2005) https://doi.org/10.1175/JAM2246.1
        ! and Belo-Pereira (2015) https://doi.org/10.1002/met.1505
        !
        ! Temperature membership function
        if (temp < t1) then
          member_t = 0.0_r_def
        else if (temp >= t1 .and. temp < t2) then
          member_t = (temp-t1)/(t2-t1)
        else if (temp >= t2 .and. temp < t3) then
          member_t = 1.0_r_def
        else if (temp >= t3 .and. temp < t4) then
          member_t = 1.0_r_def-((temp-t3)/(t4-t3))
        else if (temp >= t4) then
          member_t = 0.0_r_def
        end if

        ! Cloud liquid water membership function.
        ! Can include contribution from ice and rain
        ! (depending on what gets passed in).
        if ( lwc * 1000.0_r_def < q1 ) then
          member_clw = (1.0_r_def/q1) * lwc * 1000.0_r_def
        else if ( lwc * 1000.0_r_def >=q1 ) then
          member_clw = 1.0_r_def
        end if

        ! Vertical velocity membership function
        ! N.B. there is a typo in equation 4 of
        ! https://doi.org/10.1175/WAF-D-18-0177.1
        ! it should be w3=0.05m/s, like in their Fig 2c.
        if (w_wind < w1) then
          member_w = -0.4_r_def
        else if (w_wind >= w1 .and. w_wind < w2) then
          member_w = -0.4_r_def+(0.4_r_def/(w2-w1))*(w_wind-w1)
        else if (w_wind >= w2 .and. w_wind < w3) then
          member_w = (1.0_r_def/(w3-w2))*(w_wind-w2)
        else if (w_wind >= w3) then
          member_w = 1.0_r_def
        end if

        ! NB as discussed in 2019 paper, we replace the RH membership function
        ! in their equation 1 with just the bulk cloud fraction.
        sfip = ( (tune_cr * cld_fraction ) +                                   &
                 (tune_cc * member_clw   ) +                                   &
                 (tune_cw * member_w     )                                     &
                ) * member_t

        if ( height_wth(map_wth(1) + k) > 2.0e4_r_def ) then
          ! As currently coded, we can have positive SFIP in stratosphere
          ! due to vertical motions (member_w positive) where member_t
          ! warms again above cold-point. So reset to zero above 20 km.
          sfip = 0.0_r_def
        end if

        ! Ensure bounded between 0 and 1 (member_w can create negative numbers)
        sfip = min(max(sfip,0.0_r_def),1.0_r_def)

        ! Put back into wth space
        sfip_wth( map_wth(1)+k ) = sfip

      end do

      sfip_wth( map_wth(1)+0 ) = sfip_wth( map_wth(1)+1 )

    end if

  end subroutine icing_diags_code

end module icing_diags_kernel_mod
