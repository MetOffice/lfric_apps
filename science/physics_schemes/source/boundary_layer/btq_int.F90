! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To interpolate buoyancy parameters BT and BQ from full
!  levels to half levels

! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 95
!   THIS CODE is WRITTEN to UMDP 3 PROGRAMMING STANDARDS.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module btq_int_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'BTQ_INT_MOD'
contains

subroutine btq_int (                                                           &
! in levels
 bl_levels,                                                                    &
! in fields
 z_tq,z_uv,bq,bt,bq_cld,bt_cld,a_qs,a_dqsdt,                                   &
! out fields
 bqm,btm,bqm_cld,btm_cld,a_qsm,a_dqsdtm                                        &
  )

use atm_fields_bounds_mod, only: pdims,tdims
use bl_option_mod, only: one
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! ARGUMENTS WITH intent in. IE: INPUT VARIABLES.
integer, intent(in) ::                                                         &
 bl_levels              ! in No. of atmospheric levels for which

real(kind=r_bl), intent(in) ::                                                 &
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1),        &
 bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                            ! in A buoyancy parameter for clear air
                            !    on p,T,q-levels (full levels).
 bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                            ! in A buoyancy parameter for clear air
                            !    on p,T,q-levels (full levels).
 bq_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                            ! in A buoyancy parameter for cloudy air
                            !    on p,T,q-levels (full levels).
 bt_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                            ! in A buoyancy parameter for cloudy air
                            !    on p,T,q-levels (full levels).
 a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                            ! in Saturated lapse rate factor
                            !    on p,T,q-levels (full levels).
 a_dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels)
                            ! in Saturated lapse rate factor
                            !    on p,T,q-levels (full levels).

! out fields
real(kind=r_bl), intent(out) ::                                                &
 bqm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),           &
                            ! out A buoyancy parameter for clear air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 btm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),           &
                            ! out A buoyancy parameter for clear air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 bqm_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels),                                                           &
                            ! out A buoyancy parameter for cloudy air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 btm_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels),                                                           &
                            ! out A buoyancy parameter for cloudy air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 a_qsm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),         &
                            ! out Saturated lapse rate factor
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 a_dqsdtm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                 &
          bl_levels)
                            ! out Saturated lapse rate factor
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
!-----------------------------------------------------------------------
!    Local and other symbolic constants :-

!  Define local storage.
!  (b) Scalars.
real(kind=r_bl) ::                                                             &
 wk,                                                                           &
             ! Temporary in weighting factor.
 wkm1,                                                                         &
             ! Temporary in weighting factor.
 weight1,weight2,weight3

integer ::                                                                     &
 i,j,                                                                          &
             ! Loop counter (horizontal field index).
 k       ! Loop counter (vertical level index).

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BTQ_INT'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! 1.  Loop round levels.
!-----------------------------------------------------------------------
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP             private(weight1,weight2,weight3,wkm1,wk,i,j,k)               &
!$OMP             SHARED(btm,bt,bqm,bq,btm_cld,bt_cld,bqm_cld,bq_cld,          &
!$OMP                    a_qsm,a_qs,a_dqsdtm,a_dqsdt,                          &
!$OMP                    bl_levels,pdims,z_tq,z_uv)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

      !---------------------------------------------------------------
      ! 1.1 Calculate buoyancy parameters at half levels,
      !     i.e. at level K-1/2, if current level is level K.
      !---------------------------------------------------------------
      weight1 = one / ( z_tq(i,j,k) -                                          &
                       z_tq(i,j,k-1))
      weight2 = z_tq(i,j,k) -                                                  &
                z_uv(i,j,k)
      weight3 = z_uv(i,j,k) -                                                  &
                z_tq(i,j,k-1)
      wkm1 = weight3 * weight1
      wk = weight2 * weight1

      btm(i,j,k-1) = wkm1*bt(i,j,k) + wk*bt(i,j,k-1)
      bqm(i,j,k-1) = wkm1*bq(i,j,k) + wk*bq(i,j,k-1)
      btm_cld(i,j,k-1) = wkm1*bt_cld(i,j,k) + wk*bt_cld(i,j,k-1)
      bqm_cld(i,j,k-1) = wkm1*bq_cld(i,j,k) + wk*bq_cld(i,j,k-1)
      a_qsm(i,j,k-1) = wkm1*a_qs(i,j,k) + wk*a_qs(i,j,k-1)
      a_dqsdtm(i,j,k-1) = wkm1*a_dqsdt(i,j,k) + wk*a_dqsdt(i,j,k-1)

    end do !
  end do !
end do ! bl_levels
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine btq_int
end module btq_int_mod
