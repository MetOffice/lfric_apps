! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! PC2 inhomogeneous forcing of ice in boundary layer

module pc2_bl_inhom_ice_mod

implicit none
!---------------------------------------------------------------------------
! Description: Called from ni_imp_ctl to do inhomogeneous forcing of ice
!              cloud that has been mixed by the boundary layer

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!---------------------------------------------------------------------------

character(len=*), parameter, private :: ModuleName='PC2_BL_INHOM_ICE_MOD'

contains

subroutine pc2_bl_inhom_ice( qcf_latest, qcf_earliest, cff_earliest,           &
                             cff_latest, cf_latest)


use atm_fields_bounds_mod, only: tdims
use pc2_constants_mod,   only: cloud_rounding_tol, ls_bl0, max_in_cloud_qcf
use science_fixes_mod,   only: l_fix_incloud_qcf
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

implicit none

real, intent(in) ::                                                            &
     qcf_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
                tdims%k_end),                                                  &
                ! in ice cloud water content current value
     qcf_earliest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  tdims%k_end),                                                &
                  ! in ice cloud water content before BL mixing
     cff_earliest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  tdims%k_end)
                  ! in ice cloud fraction before BL mixing
real, intent(in out) ::                                                        &
     cf_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
               tdims%k_end),                                                   &
               ! in out bulk cloud fraction current value to update
     cff_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
                tdims%k_end)
               ! in out ice cloud fraction current value to update

integer :: i, j, k

real :: q4, denom, deltacff, incloudice, temp3

! Max in-cloud qcf used in the denominator in the formula to
! limit qcf/cff to <= max_in_cloud_qcf.
! Should be equal to the parameter max_in_cloud_qcf, but this copy
! can be used to reproduce old buggy behaviour where it was
! accidentally typo'd as 2.0e3 instead of 2.0e-3
real :: max_in_cloud_qcf_d

character(len=*), parameter :: RoutineName='PC2_BL_INHOM_ICE'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle
!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------

if ( l_fix_incloud_qcf ) then
  ! Max allowed value of qcf/cff used in denominator should be equal
  ! to the parameter; use this under the bug-fix.
  max_in_cloud_qcf_d = max_in_cloud_qcf
else
  ! Old bug version; accidentally set to 2.0e3 instead of 2.0e-3.
  max_in_cloud_qcf_d = 2.0e3
end if

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED( tdims, qcf_latest, qcf_earliest,                                 &
!$OMP         cff_earliest, cff_latest, cf_latest, max_in_cloud_qcf_d )        &
!$OMP private( i, j, k, q4, denom, deltacff, incloudice, temp3 )
!$OMP do SCHEDULE(DYNAMIC)
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      ! Calculate Q4. Only perform the calculation if the Q4 is non-zero.
      ! Since ice is mixed by tracer mixing regardless of whether cumulus
      ! is present we still need to provide an ice cloud fraction
      ! increment below cumulus base.

      ! Calculate change in ice cloud fraction differently depending on
      ! whether QCF has increased or decreased.

      q4 = qcf_latest(i,j,k) - qcf_earliest(i,j,k)

      if (q4 > cloud_rounding_tol) then

        ! Source.

        ! Calculate the change in total cloud fraction.
        ! Use a weighted (by ice cloud fraction) average of in-cloud ice
        ! content and a fixed value to specify the plume ice content.
        ! The denominator in the deltaCf calculation is then
        ! ((qcf_earliest/cff_earliest)*cffearliest +
        ! ls_bl0*(1-cff_earliest) - qcf_earliest
        ! and the qcf_earliest terms cancel.

        denom = ls_bl0 * ( 1.0 - cff_earliest(i,j,k) )

        if ( abs(denom) > 1.0e-10 ) then

          denom = q4 / denom

          deltacff = (1.0-cff_latest(i,j,k)) * denom
          cff_latest(i,j,k) = cff_latest(i,j,k) + deltacff

          ! calculate total cf based on minimum overlap
          cf_latest(i,j,k)  = cf_latest(i,j,k)  + deltacff

        else
          cf_latest(i,j,k)  = 1.0
          cff_latest(i,j,k) = 1.0
        end if

      else if (q4 < -cloud_rounding_tol) then

        ! Sink.

        ! Given a reduction in qcf, remove some CFF in order to
        ! maintain the same in-cloud IWC.

        incloudice = qcf_latest(i,j,k) /                                       &
             max(cff_latest(i,j,k),1.0e-6)

        cf_latest(i,j,k) = cf_latest(i,j,k) +                                  &
             q4/max(incloudice, ls_bl0)

        cff_latest(i,j,k) = cff_latest(i,j,k) +                                &
             q4/max(incloudice, ls_bl0)

        if (cff_latest(i,j,k) > 0.0 .and.                                      &
             !Prevent very high in-cloud values by increasing CFF.
             (qcf_latest(i,j,k) / cff_latest(i,j,k) > max_in_cloud_qcf)) then
          temp3=cff_latest(i,j,k)
          cff_latest(i,j,k)=qcf_latest(i,j,k)/max_in_cloud_qcf_d
          temp3=cff_latest(i,j,k)-temp3
          ! Update total cloud fraction.
          cf_latest(i,j,k)=cf_latest(i,j,k)+temp3
        end if

      end if
    end do
  end do
end do
!$OMP end do
!$OMP end PARALLEL

!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
return
end subroutine pc2_bl_inhom_ice

end module pc2_bl_inhom_ice_mod
