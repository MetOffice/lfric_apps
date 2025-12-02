! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_scatter_conv_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Convert rate of change of water tracers between specific humidity and
!   mixing ratio units or vice versa if required, and uncompress the
!   water tracer arrays.
!
!   For use with the deep, shallow and congestus convection schemes.
!
!   Note, that water tracer specific humidity = tracer ratio * qv
!   and that water tracer mixing ratio = tracer ratio * mv
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'WTRAC_SCATTER_CONV_MOD'

contains

subroutine wtrac_scatter_conv(np_field, npnts, np_c, nlev, n_wtrac,            &
                              iconv_scheme, idx, mt, wtrac_e,                  &
                              dqbydt_wtrac, dqclbydt_wtrac, dqcfbydt_wtrac)

use gen_phys_inputs_mod, only: l_mr_physics
use wtrac_conv_mod,      only: conv_e_wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in)  :: np_field       ! Number of points in full field
integer, intent(in)  :: npnts          ! Number of points in segment
integer, intent(in)  :: np_c           ! Number of convective points
integer, intent(in)  :: nlev           ! Number of model layers
integer, intent(in)  :: n_wtrac        ! Number of water tracers
integer, intent(in)  :: iconv_scheme   ! Type of convective scheme
                                       !  (=1 for G-R scheme)
integer, intent(in)  :: idx(np_c)      ! Index for convection points on
                                       ! full grid

real(kind=real_umphys), intent(in)     :: mt(npnts,nlev)
                                       ! Total water mixing ratio (kg/kg)
                                       ! (rhowet/rhodry = 1 + mt)

type(conv_e_wtrac_type), intent(in) :: wtrac_e(n_wtrac)
                                       ! Structure containing
                                       ! compressed water tracer arrays

real(kind=real_umphys), intent(in out) :: dqbydt_wtrac(np_field,nlev,n_wtrac)
                              ! Rate of change of water tracer vapour (kg/kg)
real(kind=real_umphys), intent(in out) :: dqclbydt_wtrac(np_field,nlev,n_wtrac)
                              ! Rate of change of water tracer liquid (kg/kg)
real(kind=real_umphys), intent(in out) :: dqcfbydt_wtrac(np_field,nlev,n_wtrac)
                              ! Rate of change of water tracer ice (kg/kg)

! Local variables
integer :: i, k, i_wt                ! loop counters
real(kind=real_umphys)    :: denom   ! 1/denominator

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_SCATTER_CONV'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


if (iconv_scheme == 1) then  ! G-R scheme (which requires q)

  if (l_mr_physics) then  !  Conversion required mr -> q

    do i_wt = 1,n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,np_c
          denom = 1.0 + mt(idx(i),k)
          dqbydt_wtrac(idx(i),k,i_wt)   = denom * wtrac_e(i_wt)%dqbydt(i,k)
          dqclbydt_wtrac(idx(i),k,i_wt) = denom * wtrac_e(i_wt)%dqclbydt(i,k)
          dqcfbydt_wtrac(idx(i),k,i_wt) = denom * wtrac_e(i_wt)%dqcfbydt(i,k)
        end do
      end do
    end do

  else         ! input is specific humidity therefore no problems

    do i_wt = 1, n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,np_c
          dqbydt_wtrac(idx(i),k,i_wt)   = wtrac_e(i_wt)%dqbydt(i,k)
          dqclbydt_wtrac(idx(i),k,i_wt) = wtrac_e(i_wt)%dqclbydt(i,k)
          dqcfbydt_wtrac(idx(i),k,i_wt) = wtrac_e(i_wt)%dqcfbydt(i,k)
        end do
      end do
    end do
  end if      ! Test on l_mr_physics

else        ! Future schemes to be coded in mixing ratio

  if (l_mr_physics) then ! ok
    do i_wt = 1, n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,np_c
          dqbydt_wtrac(idx(i),k,i_wt)   = wtrac_e(i_wt)%dqbydt(i,k)
          dqclbydt_wtrac(idx(i),k,i_wt) = wtrac_e(i_wt)%dqclbydt(i,k)
          dqcfbydt_wtrac(idx(i),k,i_wt) = wtrac_e(i_wt)%dqcfbydt(i,k)
        end do
      end do
    end do

  else  ! Output needs to be specific humidity

    do i_wt = 1,n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,np_c
          denom =  1.0 / ( 1.0 + mt(idx(i),k) )
          dqbydt_wtrac(idx(i),k,i_wt)   = denom * wtrac_e(i_wt)%dqbydt(i,k)
          dqclbydt_wtrac(idx(i),k,i_wt) = denom * wtrac_e(i_wt)%dqclbydt(i,k)
          dqcfbydt_wtrac(idx(i),k,i_wt) = denom * wtrac_e(i_wt)%dqcfbydt(i,k)
        end do
      end do
    end do
  end if

end if   ! test on convection scheme

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_scatter_conv

! ---------------------------------------------------------------------

subroutine wtrac_scatter_conv_mid(np_field, npnts, np_c, nlev, n_wtrac,        &
                                  iconv_scheme, idx, mt, wtrac_e,              &
                                  dqbydt_wtrac, dqclbydt_wtrac, dqcfbydt_wtrac)

!
! Description:
!   Convert rate of change of water tracers between specific humidity and
!   mixing ratio units or vice versa if required, and uncompress the
!   water tracer arrays.
!
!   For use with the mid-level convection scheme only. (dqbydt due to
!   mid-level convection must be added to the existing dqbydt values as
!   this type of convection can occur at the same grid point as other types
!   of convection.)
!
!   Note, that water tracer specific humidity = tracer ratio * qv
!   and that water tracer mixing ratio = tracer ratio * mv
!

use gen_phys_inputs_mod, only: l_mr_physics
use wtrac_conv_mod,      only: conv_e_wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim


implicit none

! Subroutine arguments
integer, intent(in)  :: np_field       ! Number of points in full field
integer, intent(in)  :: npnts          ! Number of points in segment
integer, intent(in)  :: np_c           ! Number of convective points
integer, intent(in)  :: nlev           ! Number of model layers
integer, intent(in)  :: n_wtrac        ! Number of water tracers
integer, intent(in)  :: iconv_scheme   ! Type of convective scheme
                                       !  (=1 for G-R scheme)
integer, intent(in)  :: idx(np_c)      ! Index for convection points on
                                       ! full grid

real(kind=real_umphys), intent(in)     :: mt(npnts,nlev)
                                       ! Total water mixing ratio (kg/kg)
                                       ! (rhowet/rhodry = 1 + mt)

type(conv_e_wtrac_type), intent(in) :: wtrac_e(n_wtrac)
                                       ! Structure containing
                                       ! compressed water tracer arrays

real(kind=real_umphys), intent(in out) :: dqbydt_wtrac(np_field,nlev,n_wtrac)
                              ! Rate of change of water tracer vapour (kg/kg)
real(kind=real_umphys), intent(in out) :: dqclbydt_wtrac(np_field,nlev,n_wtrac)
                              ! Rate of change of water tracer liquid (kg/kg)
real(kind=real_umphys), intent(in out) :: dqcfbydt_wtrac(np_field,nlev,n_wtrac)
                              ! Rate of change of water tracer ice (kg/kg)

! Local variables
integer :: i, k, i_wt              ! loop counters
real(kind=real_umphys) :: denom    ! 1/denominator

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_SCATTER_CONV_MID'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


if (iconv_scheme == 1) then  ! G-R scheme (which requires q)

  if (l_mr_physics) then  !  Conversion required mr -> q

    do i_wt = 1,n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,np_c
          denom = 1.0 + mt(idx(i),k)
          dqbydt_wtrac(idx(i),k,i_wt)   =                                      &
                          dqbydt_wtrac(idx(i),k,i_wt)                          &
                          + denom * wtrac_e(i_wt)%dqbydt(i,k)
          dqclbydt_wtrac(idx(i),k,i_wt) =                                      &
                          dqclbydt_wtrac(idx(i),k,i_wt)                        &
                          + denom * wtrac_e(i_wt)%dqclbydt(i,k)
          dqcfbydt_wtrac(idx(i),k,i_wt) =                                      &
                          dqcfbydt_wtrac(idx(i),k,i_wt)                        &
                          + denom * wtrac_e(i_wt)%dqcfbydt(i,k)
        end do
      end do
    end do

  else         ! input is specific humidity therefore no problems

    do i_wt = 1, n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,np_c
          dqbydt_wtrac(idx(i),k,i_wt)   =                                      &
            dqbydt_wtrac(idx(i),k,i_wt)   + wtrac_e(i_wt)%dqbydt(i,k)
          dqclbydt_wtrac(idx(i),k,i_wt) =                                      &
            dqclbydt_wtrac(idx(i),k,i_wt) + wtrac_e(i_wt)%dqclbydt(i,k)
          dqcfbydt_wtrac(idx(i),k,i_wt) =                                      &
            dqcfbydt_wtrac(idx(i),k,i_wt) + wtrac_e(i_wt)%dqcfbydt(i,k)
        end do
      end do
    end do
  end if      ! Test on l_mr_physics

else        ! Future schemes to be coded in mixing ratio

  if (l_mr_physics) then ! ok
    do i_wt = 1, n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,np_c
          dqbydt_wtrac(idx(i),k,i_wt)   =                                      &
            dqbydt_wtrac(idx(i),k,i_wt)   + wtrac_e(i_wt)%dqbydt(i,k)
          dqclbydt_wtrac(idx(i),k,i_wt) =                                      &
            dqclbydt_wtrac(idx(i),k,i_wt) + wtrac_e(i_wt)%dqclbydt(i,k)
          dqcfbydt_wtrac(idx(i),k,i_wt) =                                      &
            dqcfbydt_wtrac(idx(i),k,i_wt) + wtrac_e(i_wt)%dqcfbydt(i,k)
        end do
      end do
    end do

  else  ! Output needs to be specific humidity

    do i_wt = 1,n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,np_c
          denom = 1.0 / ( 1.0 + mt(idx(i),k) )
          dqbydt_wtrac(idx(i),k,i_wt)   =                                      &
                          dqbydt_wtrac(idx(i),k,i_wt)                          &
                          + denom * wtrac_e(i_wt)%dqbydt(i,k)
          dqclbydt_wtrac(idx(i),k,i_wt) =                                      &
                          dqclbydt_wtrac(idx(i),k,i_wt)                        &
                          + denom * wtrac_e(i_wt)%dqclbydt(i,k)
          dqcfbydt_wtrac(idx(i),k,i_wt) =                                      &
                          dqcfbydt_wtrac(idx(i),k,i_wt)                        &
                          + denom * wtrac_e(i_wt)%dqcfbydt(i,k)
        end do
      end do
    end do

  end if
end if   ! test on convection scheme

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_scatter_conv_mid

end module wtrac_scatter_conv_mod
