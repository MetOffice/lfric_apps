! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Mixes the subcloud temperature and humidity increments
!
module mix_ipert_6a_mod

use um_types, only: real_umphys

implicit none

! Description:
! Mixes the convective increments from the initial parcel
! perturbation in shallow/deep convection throughout the
! boundary layer.
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.

character(len=*), parameter, private :: ModuleName='MIX_IPERT_6A_MOD'

contains

subroutine mix_ipert_6a(npnts, nlev, nbl, n_wtrac, ntml,                       &
                        p_layer_boundaries, exner_layer_centres,               &
                        dthbydt, dqbydt, wtrac_e, flx_init,                    &
                        thpert, qpert, wtrac_p)

use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use wtrac_conv_mod,        only: l_wtrac_conv, conv_e_wtrac_type,              &
                                 conv_p_wtrac_type

implicit none

!----------------------------------------------------------------------
! Variables with intent in
!----------------------------------------------------------------------

integer, intent(in) :: npnts        ! Number of points

integer, intent(in) :: nlev         ! Number of model levels

integer, intent(in) :: nbl          ! in number of model layers
                                    ! potentially in the boundary layer

integer, intent(in) :: n_wtrac      ! Number of water tracers

integer, intent(in) :: ntml(npnts)  ! number of model layers for which
                                    ! increments are to be well-mixed.

real(kind=real_umphys), intent(in)    :: flx_init(npnts)
                                    ! the initial parcel mass flux

real(kind=real_umphys), intent(in)    :: thpert(npnts)
                                    ! the initial parcel temperature
                                    ! perturbation

real(kind=real_umphys), intent(in)    :: qpert(npnts)
                                    ! in the initial parcel humidity
                                    ! perturbation

real(kind=real_umphys), intent(in)    :: p_layer_boundaries(npnts,0:nlev)
                                    ! pressure at layer boundaries

real(kind=real_umphys), intent(in)    :: exner_layer_centres(npnts,0:nlev)
                                    ! exner pressure at layer centres

type(conv_p_wtrac_type), intent(in) :: wtrac_p(n_wtrac)
                                    ! Structure containing parcel water
                                    ! tracer fields

!----------------------------------------------------------------------
! variables with intent inout
!----------------------------------------------------------------------

real(kind=real_umphys), intent(in out) :: dthbydt(npnts,nlev)
                                    ! increment to potential
                                    ! temperature due to convection

real(kind=real_umphys), intent(in out) :: dqbydt(npnts,nlev)
                                    ! increment to mixing ratio
                                    ! due to convection

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                    ! Structure containing environment
                                    ! water tracer fields
!----------------------------------------------------------------------
! Variables that are locally defined
!----------------------------------------------------------------------

integer :: i,k,i_wt                 ! loop counters

real(kind=real_umphys)    :: delpsum(npnts)
                                    ! summation of model layer thicknesses
                                    ! with height. (pa)

real(kind=real_umphys)    :: delpk(npnts,nbl)
                                    ! difference in pressure across a
                                    ! layer (pa)

real(kind=real_umphys)    :: delpexsum(npnts)
                                    ! summation of model layer thicknesses
                                    ! multiplied by the exner pressure (pa).

real(kind=real_umphys)    :: dthbydt_exdp(npnts)
                                    ! increment to potential temperature
                                    ! due to initial perturbation at ntml
                                    ! multiplied by the layer thickness and
                                    ! exner pressure

real(kind=real_umphys)    :: dqbydt_dp(npnts)
                                    ! increment to mixing ratio
                                    ! due to initial perturbation at ntml
                                    ! multiplied by the layer thickness

real(kind=real_umphys), allocatable :: dqbydt_dp_wtrac(:,:)
                                    ! increment to water tracer mixing ratio
                                    ! due to initial perturbation at ntml
                                    ! multiplied by the layer thickness

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MIX_IPERT_6A'

!----------------------------------------------------------------------
! Calculate the layer pressure thickness and sum up.
!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i = 1, npnts
  delpk(i,1)   = p_layer_boundaries(i,0) - p_layer_boundaries(i,1)
  delpsum(i)   = delpk(i,1)
  delpexsum(i) = delpk(i,1) * exner_layer_centres(i,1)
end do

do k = 2, nbl
  do i = 1, npnts
    if (k  <=  ntml(i)) then
      delpk(i,k)   = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)
      delpsum(i)   = delpsum(i)   + delpk(i,k)
      delpexsum(i) = delpexsum(i) + delpk(i,k) * exner_layer_centres(i,k)
    end if
  end do
end do

!----------------------------------------------------------------------
! Calculate the potential temperature and mixing ratio increments due
! to the initial perturbation multiplied be the appropriate thickness
! nb. the delpk(ntml) in the numerator and denominator cancel
!----------------------------------------------------------------------

do i = 1, npnts
  dthbydt_exdp(i) = -flx_init(i) * thpert(i) * exner_layer_centres(i,ntml(i))
  dqbydt_dp(i)    = -flx_init(i) * qpert(i)
end do

! Water tracers
if (l_wtrac_conv) then
  allocate(dqbydt_dp_wtrac(npnts,n_wtrac))
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      dqbydt_dp_wtrac(i,i_wt) = -flx_init(i) * wtrac_p(i_wt)%qpert(i)
    end do
  end do
end if

!----------------------------------------------------------------------
! Mix the increments due to initial perturbation throughout the
! sub-cloud layer.
!----------------------------------------------------------------------

do k = 1, nbl
  do i = 1, npnts
    if (k  <=  ntml(i)) then
      dthbydt(i,k) = dthbydt(i,k) + dthbydt_exdp(i) / delpexsum(i)
      dqbydt(i,k)  = dqbydt(i,k)  + dqbydt_dp(i)    / delpsum(i)
    end if
  end do
end do

! Water tracers
if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do k = 1, nbl
      do i = 1, npnts
        if (k  <=  ntml(i)) then
          wtrac_e(i_wt)%dqbydt(i,k) = wtrac_e(i_wt)%dqbydt(i,k)                &
                                    + dqbydt_dp_wtrac(i,i_wt) / delpsum(i)
        end if
      end do
    end do
  end do
  deallocate(dqbydt_dp_wtrac)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mix_ipert_6a

end module mix_ipert_6a_mod
