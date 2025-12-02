! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Capture of raindrops by ice
! Subroutine Interface:
module lsp_capture_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_CAPTURE_MOD'

contains

subroutine lsp_capture(                                                        &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qcf, qrain, qgraup, t,                                                       &
                                          ! Water contents and temp
  cficei,                                                                      &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rainfrac, rain_liq, rain_mix,                                                &
                                          ! Rain fraction information
  rain_ice, rain_clear,                                                        &

  rho, rhor, m0, tcg, tcgi,                                                    &
                                          ! Parametrization information
  corr, dhir, ice_nofall, rain_nofall,                                         &
  lfrcp , ice_type,                                                            &
                                          ! Microphysical information
  l_psd,                                                                       &
                                          ! Code options
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  l_use_agg_vt,                                                                &
                                          ! Fallspeed branch choice
  dqprec_mix, dqprec_ice,                                                      &
                                          ! Total precip incs within ice-cloud
  precfrac_k,                                                                  &
                                          ! Prognostic precip fraction
  wtrac_mp_cpr_old                                                             &
                                          ! Store used for water tracers
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: zerodegc, cx, constp, bi,                                &
                      zero, one, two, small_number

! Microphysics modules
use mphys_constants_mod,  only: ice_type_offset
use mphys_inputs_mod,     only: l_mcr_qrain, l_diff_icevt,                     &
                                graupel_option, gr_srcols,                     &
                                l_mcr_precfrac, l_subgrid_graupel_frac,        &
                                i_update_precfrac, i_homog_areas, i_sg_correl

! Dr Hook modules
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Large scale precipitation modules
use lsp_moments_mod,      only: lsp_moments

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type

implicit none

! Purpose:
!   Update cloud prognostics as a result of capture of raindrops
!   by ice particles

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out a specified
!   distribution of raindrops.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.


! Subroutine Arguments

integer, intent(in) ::                                                         &
  points,                                                                      &
                        ! Number of points to calculate
  ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                        !              3 - graupel)

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  cficei(points),                                                              &
                        ! 1/Fraction of gridbox with ice cloud
    rho(points),                                                               &
                          ! Air density / kg m-3
    rhor(points),                                                              &
                          ! 1/Air density / m3 kg-1
    m0,                                                                        &
                          ! Seed ice water content / kg kg-1
    tcg(points),                                                               &
                          ! T dependent function in ice size distribution
    tcgi(points),                                                              &
                          ! 1/tcg (no units)
    corr(points),                                                              &
                          ! Fall velocity correction factor (no units)
    dhir(points),                                                              &
                          ! Depth of layer / timestep  / m s-1
    ice_nofall(points),                                                        &
                          ! Fraction of the ice-mass that is not falling out
    rain_nofall(points),                                                       &
                          ! Fraction of the rain-mass that is not falling out
    lfrcp,                                                                     &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
    one_over_tsi          ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  qcf(points),                                                                 &
                        ! Ice water content    / kg kg-1
  qrain(points),                                                               &
                        ! Rain mixing ratio / kg kg-1
  qgraup(points),                                                              &
                        ! Graupel mixing ratio / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  rainfrac(points),                                                            &
                        ! Rain fraction
  rain_liq(points),                                                            &
                        ! Overlap fraction of rain and liquid
  rain_mix(points),                                                            &
                        ! Overlap frac of rain with mixed phase cloud
  rain_ice(points),                                                            &
                        ! Overlap fraction of rain and ice
  rain_clear(points),                                                          &
                        ! Overlap fraction of rain and clear-sky
  ptransfer(points)
                        ! Mass rimed in this timestep / kg kg-1

logical, intent(in) ::                                                         &
  l_psd,                                                                       &
                        ! Use generic ice particle size distribution
  l_use_agg_vt(points)
                        ! Determines which vt-D parameters to use

! Stuff for updating the prognostic precip fraction;
! increments to rain+graupel in the mixed-phase and ice-only cloud regions
real(kind=real_lsprec), intent(in out) :: dqprec_mix(points)
real(kind=real_lsprec), intent(in out) :: dqprec_ice(points)
! Prognostic precip fraction
real(kind=real_lsprec), intent(in out) :: precfrac_k(points)

type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old
                        ! Water tracers store

! Local Variables

integer ::                                                                     &
  i,                                                                           &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

real (kind=real_lsprec) ::                                                     &
  dpr(points),                                                                 &
                        ! Transfer of mixing ratio  / kg kg-1
  dqcf(points),                                                                &
                        ! Increment to qcf / kg kg-1
  vr1,                                                                         &
                        ! Average fall speed of raindrops  / m s-1
  vi1,                                                                         &
                        ! Average fall speed of ice particles  / m s-1
  fv1(points),                                                                 &
                        ! Average fall speed difference between
                        ! ice particles and raindrops  / m s-1
  lamr1(points),                                                               &
                        ! Reciprocal of slope parameter in raindrop
                        ! size distribution  / m
  lami1(points),                                                               &
                        ! Reciprocal of slope parameter in ice
                        ! particle size distribution  / m
  lam4,                                                                        &
                        ! 4.0 * lamr1
  m_bi_rat,                                                                    &
                        ! Ratio of moments m_bi_1 and m_bi
  lamfac1,                                                                     &
                        ! Combination of lamr1 and lamr2
  m_0(points), m_1(points), m_2(points),                                       &
                        ! Moments of the ice particle size distributn.
  m_bi(points), m_bi_1(points),                                                &
                        ! bi and bi+1 moment of the generic ice psd.
  m_bi_di(points),                                                             &
                        ! bi+di moment of the generic ice size distn
  m_bic_dic(points)
                        ! bic+dic moment of the generic ice size distn

integer, parameter :: it_cry = 0 ! defines ice_type to be crystals

real (kind=real_lsprec), parameter    :: gr_thr = 1.0e-4_real_lsprec
                        ! threshold for converting rain water to
                        ! graupel

! Fraction of increment occuring in the ice-only vs mixed-phase cloud
real(kind=real_lsprec) :: frac

! Amount of qcf that is not falling out
real (kind=real_lsprec) :: qcf_nofall(points)

! Local compression variable
integer ::                                                                     &
  npts,                                                                        &
                        ! Number of compressed points to compute
  c,                                                                           &
                        ! Compressed point index
  ix(points)
                        ! Original index

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_CAPTURE'


    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif


cry_offset = ice_type * ice_type_offset

! Set qcf amount that is not falling out
do i = 1, points
  qcf_nofall(i) = qcf(i) * ice_nofall(i)
end do

    !-----------------------------------------------
    ! Calculate moments of size distribution if appropriate
    !-----------------------------------------------
if (l_psd) then
      ! Calculate the 0th, 1st, 2nd and bi+di (cx(82)) moments of the
      ! ice particle size distribution
  call lsp_moments(points,rho,t,qcf_nofall,cficei,zero,m_0)
  call lsp_moments(points,rho,t,qcf_nofall,cficei,one,m_1)
  call lsp_moments(points,rho,t,qcf_nofall,cficei,two,m_2)
  call lsp_moments(points,rho,t,qcf_nofall,cficei,bi,m_bi)
  call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(83),m_bi_1)
  if (.not. l_diff_icevt) then
          ! Only one set of vt-D parameters
    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(82),m_bi_di)
  else
    ! The vt-D relation which gives the
    ! least mass weighted mean fallspeed will be used so
    ! calculate both required moments
    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(182),m_bic_dic)
                        ! psd moment with crystal parameters
    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(82),m_bi_di)
                        ! psd moment with aggregate parameters
  end if
end if

! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  if (qcf_nofall(i)  >   m0 .and. qrain(i)  >   zero                           &
      .and. (rain_mix(i)+rain_ice(i))  >   zero                                &
      .and. t(i)  <   zerodegc) then

    npts = npts + 1
    ix(npts) = i

  end if ! qcf_nofall(i) >  m0 etc.

end do

! In all the following loops, the 'i' variable is used for the
! original index while the 'c' variable is used for the compressed
! index.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)


      !-----------------------------------------------
      ! Calculate reciprocal of lambda for rain
      !-----------------------------------------------
  if (l_mcr_qrain) then
        ! rain is a mixing ratio (kg kg-1)
    lamr1(i) = (rho(i) * constp(50) * qrain(i)*rain_nofall(i) / rainfrac(i))   &
               **(cx(52))

  else
        ! rain is a flux (kg m-2 s-1)
    lamr1(i) = (qrain(i) * rho(i) * dhir(i)                                    &
               /( rainfrac(i)*constp(42)*corr(i)) )**(cx(42))

  end if  ! l_mcr_qrain

      !-----------------------------------------------
      ! Calculate reciprocal of lambda for ice crystals
      !-----------------------------------------------
  if (.not. l_psd) then
    lami1(i) = (rho(i) * qcf_nofall(i) * cficei(i)                             &
           *constp(5+cry_offset)*tcgi(i))**(-cx(7+cry_offset))
  end if


end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

      !-----------------------------------------------
      ! Calculate rain mass-weighted fallspeed
      !-----------------------------------------------
  !          If (l_mcr_qrain) then
  !            ! rain is a mixing ratio (kg kg-1)
  !            vr1    = (constp(41)/6.0) * corr(i) *                       &
  !     &             ( rho(i)*constp(50)*qrain(i)/rainfrac(i) )**cx(51)

  if (l_mcr_qrain) then

    vr1 = ((lamr1(i)**cx(45))/(constp(53)*rainfrac(i))) *                      &
     ( ( constp(54) / ((lamr1(i)+cx(56))**cx(59))) +                           &
       ( constp(55) / ((lamr1(i)+cx(57))**cx(60))))
  else

    vr1 = ((lamr1(i)**cx(45))/constp(53))   *                                  &
     ( ( constp(54) / ((lamr1(i)+cx(56))**cx(59))) +                           &
        ( constp(55) / ((lamr1(i)+cx(57))**cx(60))))

  end if !l_mcr_qrain

      !-----------------------------------------------
      ! Calculate ice mass-weighted fallspeed
      !-----------------------------------------------
  if (l_psd) then
        ! Use the generic PSD
    if (.not. l_diff_icevt) then
       ! Only one set of vt-D parameters
       ! constp(82) = ci*ai
      vi1 = constp(82) * corr(i) * m_bi_di(i)                                  &
          / (rho(i) * qcf_nofall(i) * cficei(i))
    else
      if (l_use_agg_vt(i)) then
          ! Use aggregate parameters
        vi1 = constp(82) * corr(i) * m_bi_di(i)                                &
            / (rho(i) * qcf_nofall(i) * cficei(i))
      else
          ! Use crystal parameters
        vi1 = constp(182) * corr(i) * m_bic_dic(i)                             &
            / (rho(i) * qcf_nofall(i) * cficei(i))
      end if
    end if

  else
    vi1 = constp(4+cry_offset) * corr(i) * (rho(i) * qcf_nofall(i) * cficei(i) &
           * constp(5+cry_offset) * tcgi(i))**cx(3+cry_offset)
  end if  ! l_psd

      !-----------------------------------------------
      ! Estimate the mean absolute differences in velocities
      !-----------------------------------------------
  fv1(i) = max(abs(vr1-vi1),(vr1+vi1)/8.0_real_lsprec)

end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

      !------------------------------------------------
      ! Calculate transfer
      !------------------------------------------------
  if (l_psd) then

        ! Use the generic ice particle size distribution
        ! constp(86)=pi**2/24 x1r rho_water
        ! constp(87)=gamma(6+x4r)
        ! constp(88)=2 gamma(5+x4r)
        ! constp(89)=gamma(4+x4r)

    dpr(i) = constp(86) * fv1(i) *  timestep * rhor(i) *                       &
            (rain_mix(i)+rain_ice(i)) * lamr1(i)**(-cx(46)) *                  &
            (lamr1(i)**cx(43) * constp(87) * m_0(i) +                          &
             lamr1(i)**cx(44) * constp(88) * m_1(i) +                          &
             lamr1(i)**cx(45) * constp(89) * m_2(i) )

  else

    lamfac1 = constp(10+cry_offset) * constp(43) *                             &
             (lamr1(i)**cx(43) * lami1(i)**cx(8+cry_offset)) +                 &
             constp(11+cry_offset) * constp(44) *                              &
             (lamr1(i)**cx(44) * lami1(i)**cx(9+cry_offset)) +                 &
             constp(12+cry_offset) * constp(45) *                              &
             (lamr1(i)**cx(45) * lami1(i)**cx(10+cry_offset))

    dpr(i) = tcg(i) * constp(13+cry_offset) *                                  &
           lami1(i)**(-cx(11+cry_offset)) *                                    &
           lamr1(i)**(-cx(46)) * fv1(i) * lamfac1 *                            &
           timestep * rhor(i) * (rain_mix(i)+rain_ice(i))

  end if  ! l_psd

end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

        ! Limit transfer to the mass of rain that is available

  if ( rain_liq(i) > zero .or. rain_clear(i) > zero) then

    dpr(i) = min(dpr(i),qrain(i) *                                             &
                      (rain_mix(i)+rain_ice(i))/rainfrac(i))

  else

        ! All rain is through ice-only or mixed-phase cloud, so
        ! (rain_mix + rain_ice) / rainfrac must be 1.

    dpr(i) = min(dpr(i),qrain(i))

  end if

  ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi


end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

      !------------------------------------------------
      ! Adjust ice and rain contents
      !------------------------------------------------

  if ( graupel_option == gr_srcols ) then

    ! Move contents to graupel rather than snow/ice

    if (l_psd) then

      ! With the generic ice PSD active, we examine the
      ! mass-weighted mean diameters of the ice and rain.
      ! If the ice is bigger, we form snow. If the rain is
      ! bigger we form graupel.

      ! Calculate m_bi+1 / m_bi (mass-weighted diameter of ice)

      m_bi_rat = m_bi_1(i) / m_bi(i)

      ! Calculate 4 /lambda  (mass-weighted diameter of rain)
      ! As lamr1 = 1/lambda can just multiply this by 4:

      lam4 = 4.0_real_lsprec * lamr1(i)

      if ( m_bi_rat >= lam4 ) then

        ! Snow m.m.d. bigger than rain
        ! Put transfer quantities as snow

        qcf(i) = qcf(i) + dpr(i)
        dqcf(i) = dpr(i)

      else

        ! Rain m.m.d. bigger than snow
        ! Put transfer quantities as graupel

        qgraup(i) = qgraup(i) + dpr(i)
        dqcf(i) = zero

      end if

    else ! l_psd

      ! Without the generic ice PSD, we have two rain
      ! categories. In this case, we will follow a
      ! similar method to the Met Office LEM, transferring
      ! crystals colliding with small rain amounts to
      ! crystals and anything else to graupel.

      if ( ice_type == it_cry .and. qrain(i)*rain_nofall(i) < gr_thr ) then

        ! For crystals capturing low amounts of liquid
        ! water we shall move the result of the capture
        ! process to snow.

        ! If qr >= gr_thr (heavy rain) the quantities are
        ! assumed to freeze rapidly and produce graupel.

        qcf(i) = qcf(i) + dpr(i)
        dqcf(i) = dpr(i)

      else

        ! The capture process is one of
        ! i  ) Crystals capturing large rain amounts or
        ! ii ) aggregates capturing any rain amount
        ! so we shall send the result of this capture process
        ! to graupel.

        qgraup(i) = qgraup(i) + dpr(i)
        dqcf(i) = zero

      end if ! ice_type eq 0

    end if ! l_psd

  else    ! l_mcr_qgraup / l_sr2graup

    ! Graupel is not active, so move the amounts to
    ! the appropriate ice category.

    qcf(i) = qcf(i) + dpr(i)
    dqcf(i) = dpr(i)

  end if  ! l_mcr_qgraup / l_sr2graup


end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

  ! Now remove the appropriate transfer amount
  ! from the rain category and adjust the
  ! temperature

  qrain(i) = qrain(i) - dpr(i)
  t(i)     = t(i)     + dpr(i) * lfrcp

end do  ! Points

if (l_wtrac) then     ! Store phase change amount for water tracer use
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts
    i = ix(c)
    wtrac_mp_cpr_old%qchange(i) = dpr(i)
  end do
end if

if ( l_mcr_precfrac ) then
  ! If using prognostic precipitation fraction...
  ! Need to store contribution of capture of rain by ice to the
  ! total precip (rain+graupel) increment in the ice-cloud partition.

  if ( l_subgrid_graupel_frac ) then
    ! Graupel is included in the prognostic precip fraction;
    ! so only conversion of rain to qcf alters total precip (qrain+qgraup);
    ! conversion of rain to graupel makes no difference.

    if ( i_update_precfrac == i_homog_areas ) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do c = 1, npts
        i = ix(c)
        ! Increment applies in both the rain_ice and rain_mix area partitions;
        ! divvy it up between the two in proportion to area
        frac = rain_ice(i) / max( rain_ice(i) + rain_mix(i), small_number )
        dqprec_ice(i) = dqprec_ice(i) -      frac  * dqcf(i)
        dqprec_mix(i) = dqprec_mix(i) - (one-frac) * dqcf(i)
      end do
    else if ( i_update_precfrac == i_sg_correl ) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do c = 1, npts
        i = ix(c)
        if ( dqcf(i) > zero ) then
          ! Use fraction of total rain area in which capture is occuring
          ! to blend between homogeneous and anti-correlated increment:
          frac = one - (rain_ice(i) + rain_mix(i)) / rainfrac(i)
          ! In the limit that capture occurs in the whole of rainfrac (frac=0),
          ! assume homogeneous so precfrac stays the same.
          ! In the limit that capture occurs in a small fraction of rainfrac
          ! (frac close to 1), assume precfrac decreases in proportion to
          ! precip mass.
          precfrac_k(i) = precfrac_k(i)                                        &
                        * ( one - frac*dqcf(i)/( qrain(i) + dqcf(i)            &
                                               + max( qgraup(i), zero ) ) )
          ! qrain guaranteed to be > 0 as is a condition for being in the list
          ! of points ix.  But need to limit graupel for safety.
        end if
      end do
    end if

  else
    ! Graupel not included in the prognostic precip fraction;
    ! all capture of rain reduces the amount of "precip".

    if ( i_update_precfrac == i_homog_areas ) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do c = 1, npts
        i = ix(c)
        ! Increment applies in both the rain_ice and rain_mix area partitions;
        ! divvy it up between the two in proportion to area
        frac = rain_ice(i) / max( rain_ice(i) + rain_mix(i), small_number )
        dqprec_ice(i) = dqprec_ice(i) -      frac  * dpr(i)
        dqprec_mix(i) = dqprec_mix(i) - (one-frac) * dpr(i)
      end do
    else if ( i_update_precfrac == i_sg_correl ) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do c = 1, npts
        i = ix(c)
        ! Use fraction of total rain area in which capture is occuring
        ! to blend between homogeneous and anti-correlated increment
        frac = one - (rain_ice(i) + rain_mix(i)) / rainfrac(i)
        ! In the limit that capture occurs in the whole of rainfrac (frac=0),
        ! assume homogeneous so precfrac stays the same.
        ! In the limit that capture occurs in a small fraction of rainfrac
        ! (frac close to 1), assume precfrac decreases in proportion to
        ! precip mass.
        precfrac_k(i) = precfrac_k(i)                                          &
                      * ( one - frac*dpr(i)/( qrain(i) + dpr(i) ) )
      end do
    end if

  end if  ! ( l_subgrid_graupel_frac )

end if  ! ( l_mcr_precfrac )


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_capture
end module lsp_capture_mod
