! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_FIT_FSAT-------------------------------------------
!
!   Purpose: To speed up the large scale hydrology code (LTOP=true)
!            dramatically. This is done by fitting exponential
!            functions to the incomplete gamma function for each grid
!            box and the complete range of possible "water table"
!            (top_crit) cases - see documentation.
!            Estimates the fitted parameters for Fsat=function(ZW)
!            and  Fwet=function(ZW) for each land grid point.
!            (Calculating the incomplete gamma function for each grid
!            box at each time step was very time consuming).
!                                                             !
! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25
!
! Code Description:
!   Language: Originally FORTRAN 77 + common extensions.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

module calc_fit_fsat_mod

implicit none

character(len=*), parameter, private :: ModuleName = 'CALC_FIT_FSAT_MOD'

contains

subroutine calc_fit_fsat(soil_pts,soil_index,npnts                             &
  ,fexp,ti_mean,ti_sig,gamtot,zdepth                                           &
  ,a_fsat,c_fsat,a_fwet,c_fwet)


use calc_fsat_mod,         only: calc_fsat

use umPrintMgr, only:                                                          &
    umPrint,                                                                   &
    umMessage,                                                                 &
    PrintStatus,                                                               &
    PrStatus_Diag

use Ereport_Mod, only:                                                         &
  Ereport

! Replaces c_topog.h
use jules_hydrology_mod, only: zw_max, nfita

use errormessagelength_mod, only: errormessagelength

use yomhook,   only: & ! DrHook
    lhook, dr_hook

use parkind1,  only:                                                           &
    jprb, jpim

implicit none

! Subroutine arguments
!   Scalar arguments with intent(in) :
integer ::                                                                     &
 npnts                                                                         &
                  ! in No. of land points.
,soil_pts         ! in No. of land soil points.

real ::                                                                        &
 zdepth           ! in Standard Soil model DEPTH.

!   Array arguments with intent(in) :
integer ::                                                                     &
 soil_index(npnts)! in Array of soil points.

real ::                                                                        &
 ti_mean(npnts)                                                                &
                  ! in Gridbox mean topographic index.
,ti_sig(npnts)                                                                 &
                  ! in Std. deviation in topographic index.
,fexp(npnts)                                                                   &
                  ! in Exp. decay in deep layer.
,gamtot(npnts)    ! in Total gamma function.

!   Array arguments with intent(out) :
real ::                                                                        &
 a_fsat(npnts)                                                                 &
                  ! out Fitting parameter for Fsat.
,c_fsat(npnts)                                                                 &
                  ! out Fitting parameter for Fsat.
,a_fwet(npnts)                                                                 &
                  ! out Fitting parameter for Fwet.
,c_fwet(npnts)  ! out Fitting parameter for Fwet.

! Local scalars:
integer :: nzw ! Number of ZW values used in fit.
parameter(nzw=200)  ! Maximum value for a significant improvement
!                         ! in the fit.

integer ::                                                                     &
 i,j,iz                                                                        &
             ! Loop counters.
,ifita
             ! Loop counters.


real :: dzw  ! WORK ZW increment ; defined by ZW_MAX and NZW.

real ::                                                                        &
 rms                                                                           &
             ! WORK RMS errors for given fsat fit values.
,rmsw                                                                          &
             ! WORK RMS errors for given fwet fit values.
,rmsold                                                                        &
             ! WORK RMS errors for given fsat fit values.
             !      for best fit so far.
,rmswold                                                                       &
             ! WORK RMS errors for given fwet fit values
             !      for best fit so far.
,cfitmin                                                                       &
             ! WORK Minimum possible value for Cfit.
,cfitmax                                                                       &
             ! WORK Maximum possible value for Cfit.
,cfit                                                                          &
             ! WORK CFit value for given loop.
,thr_err
             ! WORK Error threshold value


parameter(cfitmin=0.0)
parameter(thr_err=5.0e-3)

! Local arrays:
real ::                                                                        &
 fsat_calc(nzw)                                                                &
                        ! WORK Surface saturation fraction.        OMP Private
,fsat_fit(nzw)                                                                 &
                        ! WORK Fitted surface saturation fraction. OMP Private
,fwet_calc(nzw)                                                                &
                        ! WORK Wetland fraction.                   OMP Private
,fwet_fit(nzw)                                                                 &
                        ! WORK Fitted wetland fraction.            OMP Private
,top_crit(nzw)                                                                 &
                        ! WORK log(QBASE_MAX/QBASE) -see document. OMP Private
,dumzw(nzw)
                        ! WORK Dummy water table depth (m).        OMP Shared

real ::                                                                        &
 dumfsat(1)                                                                    &
                        ! WORK Dummy surface saturation fraction.
,dumfwetl(1)
                        ! WORK Dummy wetland fraction.

real ::                                                                        &
 top_crit1z(nzw)                                                               &
                        ! WORK As above but for an individual zw.
,top_min                                                                       &
                        ! WORK value for when zw=zw_max.
,wutot(1)                                                                      &
                        ! WORK Dummy (set to 1.0).
,gamtot1(1)             ! Temp variable to allow interface with calc_fsat

integer :: soil_index1(1) = 1 !For the call to calc_fsat to pass a single 1 in
                              !array format
integer :: indexi(soil_pts)

logical :: L_Gamtot=.false.

real :: temp1
integer :: jtot
integer :: errorstatus
integer            ::     error_count
integer, parameter :: max_error_count = 50
character (len=errormessagelength) :: cmessage
character (len=*), parameter :: RoutineName ='CALC_FIT_FSAT'
! Dr Hook
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_count = 0

cfitmax=0.15*nfita

! Define the water table depths to be used in the fitting process:
dzw=1.0/real(nzw)*zw_max
do iz=1,nzw
  dumzw(iz)=real(iz-1)*dzw
end do

jtot = 0
do j=1,soil_pts
  i=soil_index(j)
  if (ti_mean(i) >  0.0 .and. ti_sig(i) >  0.0 ) then
    jtot = jtot + 1
    indexi(jtot) = i
  end if
end do

! Calculate TOP_CRIT for the water tables depths:
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(soil_pts,zw_max,zdepth,l_gamtot,nfita,cfitmax,error_count,        &
!$OMP        soil_index1,                                                      &
!$OMP        soil_index,ti_mean,ti_sig,fexp,gamtot,                            &
!$OMP        a_fsat,a_fwet,c_fsat,c_fwet,                                      &
!$OMP        dumzw, indexi, jtot, temp1, top_min,                              &
!$OMP        PrintStatus                                                       &
!$OMP        )                                                                 &
!$OMP private(j,i,iz,ifita,                                                    &
!$OMP         cfit,rmsold,rmswold,rms,rmsw,                                    &
!$OMP         top_crit1z,wutot,gamtot1,dumfsat,dumfwetl,                       &
!$OMP         fsat_calc,fwet_calc,top_crit,fsat_fit,fwet_fit,                  &
!$OMP         ErrorStatus,CMessage                                             &
!$OMP         )

!Privately defined for each thread
wutot(1)    = 1.0
dumfsat(1)  = 0.0
dumfwetl(1) = 0.0

!$OMP do SCHEDULE(STATIC)
do j=1,jtot
  i = indexi(j)
  top_min = 1.0/fexp(i)*exp(-fexp(i)*(zw_max-zdepth))
  temp1 = zdepth+1.0/fexp(i)-top_min

  do iz=1,nzw

    if (dumzw(iz) <= zdepth) then
      top_crit1z(iz) = -log(1-(dumzw(iz)/temp1))
    else
      top_crit1z(iz) =  log(temp1/                                             &
                       (1.0/fexp(i) * exp(-fexp(i) * (dumzw(iz)-zdepth))       &
                       -top_min))
    end if
  end do

  do iz=1,nzw

    ! Calculate FSAT and FWET for one ZW at one soil pnt:

    !This call to calc_fsat is somewhat abused, so we need to use arguments
    !that are arrays of size 1
    gamtot1(1) = gamtot(i)  !INOUT, so remember to copy back afterwards

    call calc_fsat(1,1,l_gamtot,soil_index1,ti_mean(i),ti_sig(i),              &
                   top_crit1z(iz),wutot,gamtot1,dumfsat,dumfwetl)

    gamtot(i) = gamtot1(1)

    fsat_calc(iz)=dumfsat(1)
    fwet_calc(iz)=dumfwetl(1)
    top_crit(iz)=top_crit1z(iz)
    if (iz == 1) then  ! Values at zw=0m
      a_fsat(i)=fsat_calc(iz)
      a_fwet(i)=fwet_calc(iz)
    end if
  end do

  rmsold=1.0e10
  rmswold=1.0e10

  do ifita=1,nfita
    cfit=cfitmax*(ifita)/real(nfita)

    ! This isnt really root mean squared - just a measure of fitness.
    rms=0.0
    rmsw=0.0
    !TOP_CRIT=TI_MAX when zw=zw_max
    do iz=1,nzw
      fsat_fit(iz)=a_fsat(i)*exp(-cfit*top_crit(iz))
      fwet_fit(iz)=a_fwet(i)*exp(-cfit*top_crit(iz))
      rms=rms+(fsat_calc(iz)-fsat_fit(iz))**2
      rmsw=rmsw+(fwet_calc(iz)-fwet_fit(iz))**2
    end do            !ZW
    rms=rms/real(nzw)
    rmsw=rmsw/real(nzw)

    if (rms < rmsold) then
      rmsold=rms
      c_fsat(i)=cfit
    end if
    if (rmsw <  rmswold) then
      rmswold=rmsw
      c_fwet(i)=cfit
    end if
  end do

  if (rmsold >= thr_err**2) then
    if (c_fsat(i) <= cfitmin .or. c_fsat(i) >= cfitmax) then

      do iz=1,nzw
        fsat_fit(iz)=a_fsat(i)*exp(-c_fsat(i)*top_crit(iz))
      end do               !ZW

      write(umMessage,'(a,i4,3f20.10)')                                        &
        'error CFIT FSAT',i,c_fsat(i),cfitmin,cfitmax
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,3f20.10)')'fsat_calc='                               &
        ,fsat_calc(1),fsat_calc(3)                                             &
        ,fsat_calc(5)
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,3f20.10)')'fsat_fit='                                &
        ,fsat_fit(1),fsat_fit(3),                                              &
        fsat_fit(5)
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,f20.10)')'RMS=',sqrt(rmsold)
      call umPrint(umMessage,src='calc_fit_fsat')
      ErrorStatus = 35
      write(CMessage, '(A)') 'Error in CFIT FSAT in LSH model setup'

      call Ereport ( RoutineName, ErrorStatus, CMessage)
    end if
  end if

  if (rmswold >= thr_err**2) then
    if (c_fwet(i) <= cfitmin .or. c_fwet(i) >= cfitmax) then

      do iz=1,nzw
        fsat_fit(iz)=a_fsat(i)*exp(-c_fsat(i)*top_crit(iz))
        fwet_fit(iz)=a_fwet(i)*exp(-c_fwet(i)*top_crit(iz))
      end do               !ZW

      write(umMessage,'(a,i4,3f20.10)')                                        &
        'error CFIT FWET',i,c_fwet(i),cfitmin,cfitmax
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,3f20.10)')'fwet_calc='                               &
       ,fwet_calc(1),fwet_calc(3)                                              &
       ,fwet_calc(5)
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,3f20.10)')'fwet_fit='                                &
       ,fwet_fit(1),fwet_fit(3)                                                &
       ,fwet_fit(5)
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,f20.10)')'RMSW=',sqrt(rmswold)
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,3f20.10)')'(fsat_calc=)'                             &
       ,fsat_calc(1),fsat_calc(3)                                              &
       ,fsat_calc(5)
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,3f20.10)')'(fsat_fit=)'                              &
       ,fsat_fit(1),fsat_fit(3)                                                &
       ,fsat_fit(5)
      call umPrint(umMessage,src='calc_fit_fsat')
      write(umMessage,'(a,f20.10)')'(RMS=)',sqrt(rmsold)
      call umPrint(umMessage,src='calc_fit_fsat')
      ErrorStatus = 40
      write(CMessage, '(A)') 'Error in CFIT FWET in LSH model setup'

      call Ereport ( RoutineName, ErrorStatus, CMessage)
    end if
  end if

  if (PrintStatus >= PrStatus_Diag) then
!$OMP CRITICAL (error_printing)
    if ( (rmsold >= thr_err**2 .or. rmswold >= thr_err**2)                     &
    .and. (error_count <= max_error_count)  ) then
      error_count = error_count + 1
      if (error_count < max_error_count) then
        write(umMessage,'(a,2f20.10)')                                         &
          'Warning LSH RMS Error in fit:'                                      &
          ,sqrt(rmsold),sqrt(rmswold)
        call umPrint(umMessage,src='calc_fit_fsat')
      else if (error_count == max_error_count) then
        write(umMessage,'(a)')                                                 &
          'LSH RMS Error printing DISCONTINUED.'
        call umPrint(umMessage,src='calc_fit_fsat')
      end if
    end if
!$OMP end CRITICAL (error_printing)
  end if

end do                     ! NPNTS
!$OMP end do
!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine calc_fit_fsat
end module calc_fit_fsat_mod
