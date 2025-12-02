! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module ls_cld_c_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='LS_CLD_C_MOD'

contains
!  Large-scale Cloud Scheme.
! Subroutine Interface:
! ======================================================================

!  Large-scale Cloud Scheme Compression routine (Cloud points only).
! Subroutine Interface:
subroutine ls_cld_c(                                                           &
 p_f,rhcrit,qsl_f,qn_f,q_f,t_f,                                                &
 qcl_f,cf_f,grid_qc_f,bs_f,                                                    &
 indx,points,rhc_row_length,rhc_rows,                                          &
 bl_levels,k, l_mixing_ratio)

use water_constants_mod,  only: lc
use planet_constants_mod, only: lcrcp, r, repsilon
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim
use atm_fields_bounds_mod,only: tdims
use cloud_inputs_mod,     only: i_eacf, all_clouds
use qsat_mod,             only: qsat_wat, qsat_wat_mix

implicit none

! Purpose: Calculates liquid cloud water amounts and cloud amounts,
!          temperature and specific humidity from cloud-conserved and
!          other model variables. This is done for one model level.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.29

!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
integer ::                                                                     &
                      !, intent(in)
  rhc_row_length,rhc_rows,                                                     &
!       No. of gridpoints being processed.
    bl_levels,                                                                 &
!       No. of boundary layer levels
    k,                                                                         &
!       Level no.
   points,                                                                     &
!       No. of gridpoints with non-zero cloud
   indx(tdims%j_len*tdims%i_len,2)
!       index for  points with non-zero cloud from lowest model level.

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 rhcrit(rhc_row_length,rhc_rows),                                              &
!       Critical relative humidity.  See the paragraph incorporating
!       eqs P292.11 to P292.14.
   p_f(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       pressure (Pa).
   qsl_f(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       saturated humidity at temperature TL, and pressure P_F
   qn_f(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       Normalised super/subsaturation ( = QC/BS).

logical ::                                                                     &
                      !, intent(in)
 l_mixing_ratio   !  Use mixing ratio formulation

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 q_f(           tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       On input : Vapour + liquid water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
   t_f(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       On input : Liquid water temperature (TL) (K).
!       On output: Temperature at processed levels (K).

real(kind=real_umphys) ::                                                      &
                      !, intent(out)
 qcl_f(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       Cloud liquid water content at processed levels (kg per kg air).
   cf_f(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Liquid cloud fraction at processed levels.
   grid_qc_f(     tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Super/subsaturation on processed levels. Input initially RMDI.
   bs_f(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       Value of bs at processed levels. Input initialized to RMDI.

!  Local parameters and other physical constants------------------------
real(kind=real_umphys) :: alphl      ! For liquid AlphaL calculation.
real(kind=real_umphys) :: wtn        ! Weighting for ALPHAL iteration
integer ::                                                                     &
 its                              ! Total number of iterations
parameter (its=5,wtn=0.75)

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
real(kind=real_umphys) ::                                                      &
 al,                                                                           &
                       ! LOCAL AL (see equation P292.6).
 alphal,                                                                       &
                       ! LOCAL ALPHAL (see equation P292.5).
 qn_adj,                                                                       &
 rhcritx          ! scalar copy of RHCRIT(I,J)
integer ::                                                                     &
 multrhc          ! Zero if (rhc_row_length*rhc_rows) le 1, else 1

!  (b) Others.
integer ::   i,ii,ij,n   ! Loop counters:I,II-horizontal field index.
!                                       : N - iteration counter.


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LS_CLD_C'

!  Local dynamic arrays-------------------------------------------------
!    7 blocks of real workspace are required.
real(kind=real_umphys) ::                                                      &
   qs,                                                                         &
!       Saturated spec humidity for temp T.
   qcn(points),                                                                &
!       Cloud water normalised with BS.
   t(points),                                                                  &
!       temperature.
   q(points),                                                                  &
!       specific humidity.
   bs(points),                                                                 &
!       Sigmas*sqrt(6): sigmas the parametric standard deviation of
!       local cloud water content fluctuations.
   alphal_nm1(points)
!       ALPHAL at previous iteration.

!- End of Header


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------
! Operate on INDEXed points with non-zero cloud fraction.
! ----------------------------------------------------------------------
if ( (rhc_row_length * rhc_rows)  >   1) then
  multrhc = 1
else
  multrhc = 0
end if
alphl=repsilon*lc/r

!        RHCRITX = RHCRIT(1,1)
! Points_do1:
!CDIR NODEP
do i = 1, points
  ii = indx(i,1)
  ij = indx(i,2)
  if ( multrhc== 1) then
    rhcritx = rhcrit(ii,ij)
  else
    rhcritx = rhcrit(1,1)
  end if
  qcn(i)= qn_f(ii,ij)
  ! ----------------------------------------------------------------------
  ! 1. Calculate ALPHAL (eq P292.5) and AL (P292.6).
  !    CAUTION: T_F acts as TL (input value) until update in final section
  !    CAUTION: Q_F acts as QW (input value) until update in final section
  ! ----------------------------------------------------------------------

  alphal = alphl * qsl_f(ii,ij) / (t_f(ii,ij) * t_f(ii,ij)) !P292.5
  al = 1.0 / (1.0 + (lcrcp * alphal))                    ! P292.6
  alphal_nm1(i) = alphal

  ! Rhcrit_if1:
  if (rhcritx  <   1.0) then
    ! ----------------------------------------------------------------------
    ! 2. Calculate BS (ie. sigma*sqrt(6), where sigma is
    !    as in P292.14) and normalised cloud water QCN=qc/BS, using eqs
    !    P292.15 & 16 if RHcrit < 1.
    ! N.B. QN (input) is initially in QCN
    ! N.B. QN does not depend on AL and so CF and QCN can be calculated
    !      outside the iteration (which is performed in LS_CLD_C).
    !      QN is > -1 for all points processed so CF > 0.
    ! ----------------------------------------------------------------------

    bs(i) = (1.0-rhcritx) * al * qsl_f(ii,ij)  ! P292.14
    if (qcn(i)  <=  0.0) then
      cf_f(ii,ij) = 0.5 * (1.0 + qcn(i)) * (1.0 + qcn(i))
      qcn(i)= (1.0 + qcn(i)) * (1.0 + qcn(i)) * (1.0 + qcn(i)) / 6.0
    else if (qcn(i)  <   1.0) then
      cf_f(ii,ij) = 1.0 - 0.5 * (1.0 - qcn(i)) * (1.0 - qcn(i))
      qcn(i)=qcn(i) + (1.0-qcn(i)) * (1.0-qcn(i)) * (1.0-qcn(i))/6.0
    else ! QN  >=  1
      cf_f(ii,ij) = 1.0
    end if ! QCN_if

    ! ----------------------------------------------------------------------
    ! 3.b If necessary, modify cloud fraction using empirically adjusted
    !     cloud fraction parametrization, but keep liquid content the same.
    ! ----------------------------------------------------------------------
    if (i_eacf >= all_clouds) then
      ! really this should be i_eacf=all_clouds or i_eacf=not_mixph,
      ! but writing it like that loses bit-comparison in highly optimised jobs

      ! Adjust QN according to EACF parametrization
      if (k <= bl_levels) then
        qn_adj=(qn_f(ii,ij)+0.184)/(1.0-0.184)
      else
        qn_adj=(qn_f(ii,ij)+0.0955)/(1.0-0.0955)
      end if

      !         Calculate cloud fraction using adjusted QN
      if (qn_adj  <=  0.0) then
        cf_f(ii,ij) = 0.5 * (1.0 + qn_adj) * (1.0 + qn_adj)
      else if (qn_adj  <   1.0) then
        cf_f(ii,ij) = 1.0 - 0.5 * (1.0 - qn_adj) * (1.0 - qn_adj)
      else ! QN_ADJ  >=  1
        cf_f(ii,ij) = 1.0
      end if ! QN_ADJ_if

    end if  ! i_eacf

  else ! i.e. if RHcrit = 1
    ! ----------------------------------------------------------------------
    ! 3.a If RHcrit = 1., all points processed have QN > 0 and CF = 1.
    ! ----------------------------------------------------------------------
    bs(i) = al
    cf_f(ii,ij) = 1.0
  end if ! Rhcrit_if1

  ! ----------------------------------------------------------------------
  ! 3.1 Calculate 1st approx. to qc (store in QCL)
  ! ----------------------------------------------------------------------

  qcl_f(ii,ij) = qcn(i) * bs(i)

  ! ----------------------------------------------------------------------
  ! 3.2 Calculate 1st approx. specific humidity (total minus cloud water)
  ! ----------------------------------------------------------------------

  q(i) = q_f(ii,ij) - qcl_f(ii,ij)

  ! ----------------------------------------------------------------------
  ! 3.3 Calculate 1st approx. to temperature, adjusting for latent heating
  ! ----------------------------------------------------------------------

  t(i) = t_f(ii,ij) + lcrcp*qcl_f(ii,ij)
end do ! Points_do1

! ----------------------------------------------------------------------
! 4. Iteration to find better cloud water values.
! ----------------------------------------------------------------------
! Its_if:
if (its  >=  2) then
  ! Its_do:
  do n = 2, its

    ! Points_do2:
    rhcritx = rhcrit(1,1)
    do i = 1, points
      ii = indx(i,1)
      ij = indx(i,2)
      if ( multrhc== 1) then
        rhcritx = rhcrit(ii,ij)
      else
        rhcritx = rhcrit(1,1)
      end if
      ! T_if:
      if (t(i)  >   t_f(ii,ij)) then
        !           NB. T > TL implies cloud fraction > 0.

        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qs,t(i),p_f(ii,ij))
        else
          call qsat_wat(qs,t(i),p_f(ii,ij))
        end if

        alphal = (qs - qsl_f(ii,ij)) / (t(i) - t_f(ii,ij))
        alphal = wtn * alphal + (1.0 - wtn) * alphal_nm1(i)
        alphal_nm1(i) = alphal
        al = 1.0 / (1.0 + (lcrcp * alphal))
        ! Rhcrit_if2:
        if (rhcritx  <   1.0) then
          bs(i) = (1.0-rhcritx) * al * qsl_f(ii,ij)
          !                                                             P292.14
        else
          bs(i) = al
        end if  ! Rhcrit_if2

        ! ----------------------------------------------------------------------
        ! 4.1 Calculate Nth approx. to qc (store in QCL).
        ! ----------------------------------------------------------------------

        qcl_f(ii,ij) = qcn(i) * bs(i)

        ! ----------------------------------------------------------------------
        ! 4.2 Calculate Nth approx. spec. humidity (total minus cloud water).
        ! ----------------------------------------------------------------------

        q(i) = q_f(ii,ij) - qcl_f(ii,ij)

        ! ----------------------------------------------------------------------
        ! 4.3 Calculate Nth approx. to temperature, adjusting for latent heating
        ! ----------------------------------------------------------------------

        t(i) = t_f(ii,ij) + lcrcp * qcl_f(ii,ij)

      end if ! T_if
    end do ! Points_do2
  end do ! Its_do
end if ! Its_if

! ----------------------------------------------------------------------
! 5. Finally scatter back cloud point results to full field arrays.
!    CAUTION: T_F updated from TL (input) to T (output)
!    CAUTION: Q_F updated from QW (input) to Q (output)
! ----------------------------------------------------------------------

!DIR$ IVDEP
! Points_do3:
do i = 1, points
  ii = indx(i,1)
  ij = indx(i,2)
  q_f(ii,ij) = q(i)
  t_f(ii,ij) = t(i)
  grid_qc_f(ii,ij) = bs(i) * qn_f(ii,ij)
  bs_f(ii,ij) = bs(i)
end do ! Points_do3


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_cld_c
! ======================================================================
end module ls_cld_c_mod
