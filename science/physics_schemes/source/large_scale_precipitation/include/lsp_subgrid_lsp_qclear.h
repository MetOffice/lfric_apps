! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Conversion between moments of PSD

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! ######################################################################
! Purpose:
! Returns the average relative humidity in the cloud-free (clear-sky)
! portion of gridboxes.
!
! Method:
!   Parametrizes the width of the vapour distribution in the part
!   of the gridbox which does not have liquid water present.
!
! Description of Code:
!   Fortran95
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
! ######################################################################

! Subroutine Arguments
!
!-----------------------------------------------------------------------
! arguments with intent in. ie: input variables.
!-----------------------------------------------------------------------

integer, intent(in) :: npnts
! Points (dimensions) being processed

real (kind=prec),   intent(in)  :: q(npnts)
! water vapour mixing ratio (kg/kg)

real (kind=prec),   intent(in)  :: qsmr(npnts)
! Saturation mixing ratio, w.r.t. ice where T < 0C, (kg/kg)

real (kind=prec),   intent(in)  :: qsmr_wat(npnts)
! Saturation mixing ratio, w.r.t. water, (kg/kg)

real (kind=prec),   intent(in)  :: qcf(npnts)
! Frozen cloud condensate (kg/kg)

real (kind=prec),   intent(in)  :: cloud_liq_frac(npnts)
! Fraction of gridbox occupied by liquid or mixed phase cloud

real (kind=prec),   intent(in)  :: cloud_frac(npnts) !
! Fraction of gridbox occupied by cloud or any type

real (kind=prec),   intent(in)  :: rhcrit(npnts)
! Critical relative humidity for large-scale cloud formation (fraction)

!-----------------------------------------------------------------------
! arguments with intent out. ie: output variables.
!-----------------------------------------------------------------------

! relative humidity of the cloud-free portion of the gridboxes (fraction)
real (kind=prec), intent(out) :: q_clear(npnts)

!-----------------------------------------------------------------------

!  Local variables

real (kind=prec) :: area_ice        ! Ice-only cloud fraction
real (kind=prec) :: width      ! Width of water vapour PDF
real (kind=prec) :: qa         ! q in the portion of gridbox not containing
                               ! liquid cloud
integer :: k       ! loop integer

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! loop over points
do k = 1, npnts

  area_ice = max(cloud_frac(k) - cloud_liq_frac(k), 0.0_prec)

  if (cloud_liq_frac(k) < 1.0_prec) then

    ! Derived from equation 34 in UMDP 26.
    qa = (q(k) - cloud_liq_frac(k)*qsmr_wat(k)) /                              &
         (1.0_prec - cloud_liq_frac(k))

    width = 2.0_prec *(1.0_prec-rhcrit(k))*qsmr_wat(k)                         &
            *max((1.0_prec-0.5_prec*qcf(k)/                                    &
                  (real(ice_width,kind=prec) * qsmr_wat(k))),                  &
                  0.001_prec)

    ! The full width cannot be greater than 2q because otherwise
    ! part of the gridbox would have negative q. Also ensure that
    ! the full width is not zero (possible if rhcpt is 1).
    ! 0.001 is to avoid divide by zero problems
    width = min(width , max(2.0_prec*q(k),0.001_prec*qsmr(k)))

    ! Adjust for ice-only clouds, if subgrid partition
    ! of qv is selected. Else, use the same value of qv
    ! everywhere outside the liquid-cloud

    if (area_ice  >  0.0_prec) then
      if ( l_subgrid_qv ) then
        q_clear(k) = qa - 0.5_prec*width * area_ice
      else
        q_clear(k) = qa
      end if
    else
      q_clear(k) = qa
    end if ! area_ice gt 0

  else ! cf_liq lt 1

    ! -----------------------------------------------
    ! If cloud-free or overcast then set q_clear = q
    ! -----------------------------------------------
    q_clear(k) = q(k)
  end if ! cf_liq lt 1

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
