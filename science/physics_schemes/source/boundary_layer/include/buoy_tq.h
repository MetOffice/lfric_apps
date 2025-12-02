! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

! PURPOSE: To calculate buoyancy parameters on p,T,q-levels

! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 95
!   THIS CODE is WRITTEN to UMDP 3 PROGRAMMING STANDARDS.

! ARGUMENTS WITH intent in. IE: INPUT VARIABLES.
integer, intent(in) ::                                                         &
 bl_levels              ! in No. of atmospheric levels for which

real(kind=prec), intent(in) ::                                                 &
 p(tdims%i_start:tdims%i_end,                                                  &
   tdims%j_start:tdims%j_end,0:bl_levels+1),                                   &
                                    ! in Pressure at pressure points.
 t(tdims%i_start:tdims%i_end,                                                  &
   tdims%j_start:tdims%j_end,bl_levels),                                       &
                                    ! in Temperature (K). At P points
 q(tdims_l%i_start:tdims_l%i_end,                                              &
   tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),                   &
                              ! in Sp humidity (kg water per kg air).
 qcl(tdims_l%i_start:tdims_l%i_end,                                            &
     tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),                 &
                              ! in Cloud liq water (kg per kg air).
 qcf(tdims_l%i_start:tdims_l%i_end,                                            &
     tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),                 &
                              ! in Cloud liq water (kg per kg air).
 cf_bulk(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end, bl_levels)
                              ! in Cloud fraction (decimal).

! ARGUMENTS WITH intent out. IE: OUTPUT VARIABLES.
real(kind=prec), intent(out) ::                                                &
 bq(tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end,bl_levels),                                      &
                              ! out A buoyancy parameter for clear air
 bt(tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end,bl_levels),                                      &
                              ! out A buoyancy parameter for clear air
 bq_cld(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end,bl_levels),                                  &
                            ! out A buoyancy parameter for cloudy air
 bt_cld(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end,bl_levels),                                  &
                            ! out A buoyancy parameter for cloudy air
 bq_gb(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end,bl_levels),                                   &
                              ! out A grid-box mean buoyancy parameter
 bt_gb(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end,bl_levels),                                   &
                              ! out A grid-box mean buoyancy parameter
 a_qs(tdims%i_start:tdims%i_end,                                               &
       tdims%j_start:tdims%j_end,bl_levels),                                   &
                            ! out Saturated lapse rate factor
 a_dqsdt(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end,bl_levels),                                 &
                            ! out Saturated lapse rate factor
 dqsdt(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end,bl_levels)
                            ! out Derivative of q_SAT w.r.t. T

! LOCAL VARIABLES.
real(kind=prec) ::                                                             &
 qs(tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end), & ! WORK Saturated mixing ratio.
 tmp1(tdims%i_start:tdims%i_end), & ! TEMP array to contain lc or ls

 tmp2(tdims%i_start:tdims%i_end) ! TEMP array to contain lcrcp or lsrcp

integer ::                                                                     &
  i,j,                                                                         &
  k

real(kind=prec) ::                                                             &
  bc

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! 1.  Loop round levels.
!-----------------------------------------------------------------------

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP private(i, j, k, bc, qs, tmp1, tmp2)                                     &
!$OMP SHARED(bl_levels, p, t, q, qcf, qcl, cf_bulk, bt, bq, bt_cld, bq_cld,    &
!$OMP        bt_gb, bq_gb, a_qs, a_dqsdt, dqsdt, tdims, l_mr_physics, r,       &
!$OMP        repsilon, c_virtual, etar, lcrcp, ls, lsrcp, l_noice_in_turb)

do k = 1, bl_levels

  !-----------------------------------------------------------------------
  ! 1.1 Calculate saturated specific humidity at pressure and
  !     temperature of current level.
  !-----------------------------------------------------------------------

  !No halos in these variables
  if (l_noice_in_turb) then
    ! use qsat_wat
    if ( l_mr_physics ) then
      call qsat_wat_mix(qs,t(:,:,k),p(:,:,k),tdims%i_end,tdims%j_end)
    else
      call qsat_wat(qs,t(:,:,k),p(:,:,k),tdims%i_end,tdims%j_end)
    end if
  else
    if ( l_mr_physics ) then
      call qsat_mix(qs,t(:,:,k),p(:,:,k),tdims%i_end,tdims%j_end)
    else
      call qsat(qs,t(:,:,k),p(:,:,k),tdims%i_end,tdims%j_end)
    end if
  end if ! test on l_noice_in_turb

  !   Using the temp arrays and splitting the i index helps vectorisation,
  !   Before the index split, no vectorisation was taking place as each
  !   conditional was computationally heavy with only minor differences
  !   between each side (differences being the lc, ls, lcrcp & lsrcp variables)

  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if (t(i,j,k) > tm .or. l_noice_in_turb) then
        tmp1(i) = lc
        tmp2(i) = lcrcp
        !            ...  (Clausius-Clapeyron) for T above freezing
      else
        tmp1(i) = ls
        tmp2(i) = lsrcp
        !            ...  (Clausius-Clapeyron) for T below freezing
      end if
    end do ! p_points,i

      !---------------------------------------------------------------
      ! 1.2 Calculate buoyancy parameters BT and BQ, required for the
      !     calculation of stability.
      !---------------------------------------------------------------

    do i = tdims%i_start, tdims%i_end
      bt(i,j,k) = 1.0_prec/t(i,j,k)
      bq(i,j,k) =                                                              &
        c_virtual/(1.0_prec+c_virtual*q(i,j,k)-qcl(i,j,k)-qcf(i,j,k))


      dqsdt(i,j,k) = (repsilon * tmp1(i) * qs(i,j))                            &
                   / ( r * t(i,j,k) * t(i,j,k) )


      a_qs(i,j,k) = 1.0_prec / (1.0_prec + tmp2(i)*dqsdt(i,j,k))

      a_dqsdt(i,j,k) = a_qs(i,j,k) * dqsdt(i,j,k)

      bc = tmp2(i)*bt(i,j,k) - etar*bq(i,j,k)

      !--------------------------------------------------------------
      ! 1.3 Calculate in-cloud buoyancy parameters.
      !--------------------------------------------------------------

      bt_cld(i,j,k) = bt(i,j,k) - a_dqsdt(i,j,k) * bc
      bq_cld(i,j,k) = bq(i,j,k) + a_qs(i,j,k) * bc

      !--------------------------------------------------------------
      ! 1.4 Calculate grid-box mean buoyancy parameters.
      !--------------------------------------------------------------

      bt_gb(i,j,k) = bt(i,j,k) +                                               &
                     cf_bulk(i,j,k)*( bt_cld(i,j,k) - bt(i,j,k) )
      bq_gb(i,j,k) = bq(i,j,k) +                                               &
                     cf_bulk(i,j,k)*( bq_cld(i,j,k) - bq(i,j,k) )

    end do ! p_points,i
  end do ! p_points,j
end do ! bl_levels

!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
