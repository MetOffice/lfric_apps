! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_gather_conv_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Set water tracers compressed arrays for the various types of convection
!   and convert between specific humidity and mixing ratio or vice versa
!   if required.
!
! Method:
!   Note, that water tracer specific humidity = tracer ratio * qv
!   and that water tracer mixing ratio = tracer ratio * mv
!
!   Therefore, water tracer specific humidity =
!                                water tracer mixing ratio/(1+mt)
!   and water tracer mixing ratio =
!                                water tracer specific humidity*(1+mt)
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName='WTRAC_GATHER_CONV_MOD'

contains

! Subroutine Interface:
subroutine wtrac_gather_conv(np_field, npnts, np_c, nlev, n_wtrac,             &
                             iconv_scheme, idx, mt,                            &
                             q_wtrac, qcl_wtrac, qcf_wtrac,                    &
                             wtrac_e)

use gen_phys_inputs_mod, only: l_mr_physics
use wtrac_conv_mod,      only: conv_e_wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in)  :: np_field      ! Number of points in a full field
integer, intent(in)  :: npnts         ! Number of points in segment
integer, intent(in)  :: np_c          ! Number of convective points
integer, intent(in)  :: nlev          ! Number of model layers
integer, intent(in)  :: n_wtrac       ! Number of water tracers
integer, intent(in)  :: iconv_scheme  ! Type of convective scheme
                                      !  (=1 for G-R scheme)
integer, intent(in)  :: idx(np_c)     ! Index for convection points on
                                      !   full grid

real(kind=real_umphys), intent(in)     :: mt(npnts,nlev)
                                      ! Total water mixing ratio (kg/kg)
                                      ! (rhowet/rhodry = 1 + mt)
real(kind=real_umphys), intent(in)     :: q_wtrac(np_field,nlev,n_wtrac)
                                      ! Water tracer vapour (kg/kg)
real(kind=real_umphys), intent(in)     :: qcl_wtrac(np_field,nlev,n_wtrac)
                                      ! Water tracer liq condensate (kg/kg)
real(kind=real_umphys), intent(in)     :: qcf_wtrac(np_field,nlev,n_wtrac)
                                      ! Water tracer ice condensate (kg/kg)

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                              ! Structure containing the compressed water
                              ! tracer arrays

! Local variables
integer :: j, k, i_wt               ! loop counters
real(kind=real_umphys) :: denom     ! 1/denominator

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_GATHER_CONV'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


if (iconv_scheme == 1) then  ! G-R scheme (which requires q)

  if (l_mr_physics) then  !  Conversion required mr -> q

    do i_wt = 1,n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,np_c
          denom = 1.0 / (1.0+mt(idx(j),k))
          wtrac_e(i_wt)%q(j,k)   = q_wtrac(idx(j),k,i_wt)   * denom
          wtrac_e(i_wt)%qcl(j,k) = qcl_wtrac(idx(j),k,i_wt) * denom
          wtrac_e(i_wt)%qcf(j,k) = qcf_wtrac(idx(j),k,i_wt) * denom
        end do
      end do
    end do

  else         ! input is specific humidity therefore no problems
    do i_wt = 1, n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,np_c
          wtrac_e(i_wt)%q(j,k)   = q_wtrac(idx(j),k,i_wt)
          wtrac_e(i_wt)%qcl(j,k) = qcl_wtrac(idx(j),k,i_wt)
          wtrac_e(i_wt)%qcf(j,k) = qcf_wtrac(idx(j),k,i_wt)
        end do
      end do
    end do
  end if      ! Test on l_mr_physics

else        ! Future schemes to be coded in mixing ratio

  if (l_mr_physics) then   ! Input as required

    do i_wt = 1, n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,np_c
          wtrac_e(i_wt)%q(j,k)   = q_wtrac(idx(j),k,i_wt)
          wtrac_e(i_wt)%qcl(j,k) = qcl_wtrac(idx(j),k,i_wt)
          wtrac_e(i_wt)%qcf(j,k) = qcf_wtrac(idx(j),k,i_wt)
        end do
      end do
    end do

  else      !  Conversion from q to m required

    do i_wt = 1, n_wtrac
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,np_c
          wtrac_e(i_wt)%q(j,k)   = q_wtrac(idx(j),k,i_wt  )*(1.0+mt(idx(j),k))
          wtrac_e(i_wt)%qcl(j,k) = qcl_wtrac(idx(j),k,i_wt)*(1.0+mt(idx(j),k))
          wtrac_e(i_wt)%qcf(j,k) = qcf_wtrac(idx(j),k,i_wt)*(1.0+mt(idx(j),k))
        end do
      end do
    end do

  end if    ! Test on l_mr_physics

end if   ! test on convection scheme

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_gather_conv

end module wtrac_gather_conv_mod
