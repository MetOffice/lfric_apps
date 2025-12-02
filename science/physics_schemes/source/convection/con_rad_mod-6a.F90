! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Convective cloud microphysics routine
!
module con_rad_6a_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Calculates convective cloud top, base and amount
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


character(len=*), parameter, private :: ModuleName='CON_RAD_6A_MOD'

contains

! Subroutine Interface:

subroutine con_rad_6a (k, npnts, start_lev,                                    &
                       ccwk, ccwkp1, tmp_ccwkp1, flxkp1, delpkp1,              &
                       l_q_interact, bterm,                                    &
                       tcw, cclwp, lcca, lcbase, lctop,                        &
                       iccb, icct, cca, idx, ni)

use science_fixes_mod,      only: l_fix_ccb_cct
use ereport_mod,            only: ereport
use errormessagelength_mod, only: errormessagelength
use umPrintMgr,             only: PrintStatus, PrStatus_Diag

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use planet_constants_mod, only: g

implicit none



!----------------------------------------------------------------------
! Arguments with intent(in)
!----------------------------------------------------------------------
integer,intent(in) :: k               ! Present model layer number
integer,intent(in) :: npnts           ! Vector length
integer,intent(in) :: start_lev(npnts)! Level at which convection initiated

real(kind=real_umphys),intent(in) :: ccwk(npnts)
                                    ! Total condensate in level k (kg/kg)
real(kind=real_umphys),intent(in) :: ccwkp1(npnts)
                                    ! Total condensate in level k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: tmp_ccwkp1(npnts)
                                    ! Total condensate in level k+1 (kg/kg)
                                    ! before precipitation
real(kind=real_umphys),intent(in) :: flxkp1(npnts)
                                    ! Parcel mass flux in layer k+1 (Pa/s)
real(kind=real_umphys),intent(in) :: delpkp1(npnts)
                                    ! pressure difference across layer k+1 (Pa)

logical,intent(in) :: l_q_interact  ! True if PC2 is switched on
logical,intent(in) :: bterm(npnts)  ! Mask for parcels which terminate
                                    ! in layer k+1

! ----------------------------------------------------------------------
! Arguments with intent(in/out)
! ----------------------------------------------------------------------
real(kind=real_umphys),intent(in out) :: tcw(npnts)
                                     ! Total condensed water (kg/m**2/s)
                                    ! in summed to layer k
                                    ! out summed to layer k+1
real(kind=real_umphys),intent(in out) :: cclwp(npnts)
                                     ! Condensed water path (kg/m**2)
                                    ! in summed to layer k
                                    ! out summed to layer k+1
real(kind=real_umphys),intent(in out) :: lcca(npnts)
                                     ! Lowest conv. cloud amount (%)

integer,intent(in out) :: lcbase(npnts)! Lowest conv. cloud base level
integer,intent(in out) :: lctop(npnts) ! Lowest conv. cloud top level

integer,intent(in out) :: iccb(npnts)  ! convective cloud base_level
integer,intent(in out) :: icct(npnts)  ! convective cloud top level
! Need to be inout as need to be remembered between subsequent calls to convec2

real(kind=real_umphys),intent(in out) :: cca(npnts)
                                       ! convective cloud amount (%)


! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------
integer :: i,m          ! Loop counter

integer :: errorc  ! Error code for ereport
character (len=errormessagelength) :: cmessage  ! String to store error message

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CON_RAD_6A'

! ----------------------------------------------------------------------
!  Calculate cloud base and lowest cloud base
!  when cloud base set zero total condensed water
! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  Calculate cloud base and cloud top (and lowest base and top)
! ----------------------------------------------------------------------

if ( l_fix_ccb_cct ) then
  ! Corrected version of the ccb / cct calculation; consistently
  ! set the convective base and top to just be the start-level and
  ! termination level of the ascent, ignoring cloud.
  ! This is consistent with how iccb and icct are used in the 6A
  ! convection scheme, to determine the depth of the ascent, not just
  ! the cloudy part.
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    if ( k == start_lev(i) ) then
      iccb(i) = k
      cclwp(i)  = 0.0
      tcw(i)   = 0.0
      if ( lcbase(i) == 0 )  lcbase(i) = iccb(i)
    end if
    if (bterm(i)) then
      icct(i)  = k+1
      if ( lctop(i) == 0 )  lctop(i) = icct(i)
    end if
  end do

else
  ! Original, confusing code; sets the base to be the first
  ! cloudy level in the ascent, except if using PC2 (l_q_interact),
  ! in which case the base is forced to be at the ascent start level
  ! even if the parcel contains no cloud at k or k+1.
  ! Then the conv cloud top is set to just be the top of the ascent,
  ! regardless of whether the parcel contained cloud, or if using PC2.
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    if ( ccwk(i)  <=  0.0 .and. tmp_ccwkp1(i)  >   0.0 ) then
      !       Assuming initial parcel condensate is zero
      iccb(i)   = k+1
      cclwp(i)  = 0.0

      if ( lcbase(i)  ==  0 ) then
        !         Lowest cloud base
        lcbase(i) = k+1
      end if  ! If_lcbase_1

      ! Note that, if not l_q_interact, ccwk must always be zero
      ! in the 1st level of ascent, so if any condensation occurs,
      ! the above check will set the cloud-base.  However, if using PC2
      ! (l_q_interact), the initial parcel can contain condensate (ccwk>0),
      ! so need an additional check...

    else if ( l_q_interact .and. k  ==  start_lev(i) ) then

      ! ...If the initial parcel contains condensate, then cloud-base level
      ! is the parcel initial level, not k+1.

      ! Non-zero initial parcel condensate (initialized to environment)
      iccb(i)  = k
      cclwp(i) = 0.0
      tcw(i)   = 0.0

      if ( lcbase(i)  ==  0 ) then
        !         Lowest cloud base
        lcbase(i) = k
      end if  ! If_lcbase_2

    end if

    if (bterm(i)) then
      icct(i)  = k+1
    end if

    if (bterm(i) .and.  lctop(i) == 0 ) then
      lctop(i) = k+1
    end if

  end do  ! i = 1, npnts

end if  ! l_fix_ccb_cct

! If full diagnostic output requested, or not using the fix...
if ( PrintStatus==PrStatus_Diag .or. (.not. l_fix_ccb_cct) ) then
  ! ...print a warning if cloud-base and cloud-top are not set consistently
  errorc = -111
  do m=1, ni
    i = idx(m)
    if ( bterm(i) ) then
      if (      ( iccb(i)>0   .and. icct(i)<iccb(i)    )                       &
           .or. ( icct(i)>0   .and. iccb(i)==0         )                       &
           .or. ( lcbase(i)>0 .and. lctop(i)<lcbase(i) )                       &
           .or. ( lctop(i)>0  .and. lcbase(i)==0       )                       &
           .or. ( iccb(i)>0   .and. lcbase(i)==0       )                       &
           .or. ( iccb(i)==0  .and. lcbase(i)>0        )  ) then
        ! - If the base has been set, the top must also be set (and vice versa)
        ! - The top mustn't be below the base.
        ! - If the current base/top has been set, the lowest base/top must
        !   also be set (and vice versa).
        write(cmessage,"(A,4(A,I4))")                                          &
          "Convective cloud-base and cloud-top levels are inconsistent!",      &
          " iccb = ",   iccb(i),    " icct = ",  icct(i),                      &
          " lcbase = ", lcbase(i),  " lctop = ", lctop(i)
        call ereport( RoutineName, errorc, cmessage )
      end if
    end if
  end do
end if

!DIR$ IVDEP
do m=1, ni
  i = idx(m)

  if ( flxkp1(i)  >   0.0) then

    !---------------------------------------------------------------------
    ! Sum total condensed water per second - assumes that the initial
    ! convective layer is unsaturated.
    ! Uses the CCW before precipitation
    !---------------------------------------------------------------------
    tcw(i)   = tcw(i)   + flxkp1(i) * tmp_ccwkp1(i) / g

    !---------------------------------------------------------------------
    ! Sum conv condensed water path - assumes that the initial
    ! convective layer is unsaturated
    ! Uses the CCW after precipitation
    !---------------------------------------------------------------------
    cclwp(i) = cclwp(i) + ccwkp1(i) * delpkp1(i) / g

  end if

  !---------------------------------------------------------------------
  ! calculate convective cloud amount if convection terminates in
  ! layer k and total condensed water path over a time step
  !
  ! UM documentation paper p27
  ! section (9), equation (37)
  !---------------------------------------------------------------------
  if ( bterm(i) .and. tcw(i) >  0.0 ) then

    if ( tcw(i)  <   2.002e-6 ) tcw(i) = 2.002e-6

    cca(i) = 0.7873 + 0.06 * log(tcw(i))
    if (cca(i)  >   1.0) cca(i) = 1.0

    if (lcca(i) <= 0.0) then
      lcca(i) = 0.7873 + 0.06 * log(tcw(i))
      if (lcca(i)  >   1.0) lcca(i) = 1.0
    end if

    tcw(i) = 0.0

  end if
end do ! I loop over NPNTS

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine con_rad_6a
end module con_rad_6a_mod
