! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Gravity Wave (Ultra-Simple Spectral Parametrization) Scheme Core.
! Subroutine Interface:

module gw_ussp_core_mod

use um_types, only: real_usprec

implicit none

character(len=*), parameter, private :: ModuleName='GW_USSP_CORE_MOD'

contains

subroutine gw_ussp_core(points, launchlev, k_start, k_end, ussp_launch_factor, &
  wavelstar, cgw_scale_factor, two_omega, sin_theta_latitude, rho_th,          &
  nbv, udotk, totalppn, fptot, L_add_cgw)

! purpose: This subroutine calculates the vertically propagated flux of
!          pseudomomentum due to gravity waves as parametrised by the
!          Ultra Simple Spectral gravity wave Parametrization
!         (originally Warner and McIntyre, since reengineered for use in
!          the UM and extended).

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: gravity_wave_drag

! code description:
!   language: fortran 90
!   this code is written to umdp3 programming standards.
!   documentation: Unified Model Documentation Paper 34 (Non-Orog GW)

use gw_ussp_prec_mod, only: pi
use gw_ussp_params_mod, only: idir
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

! ----------------------------------------------------------------------+-------
! Subroutines defined from modules
! ----------------------------------------------------------------------+-------

! ----------------------------------------------------------------------+-------

implicit none

! ----------------------------------------------------------------------+-------
! Subroutine arguments of gw_ussp_core
! ----------------------------------------------------------------------+-------
! Number of horizontal points in field
integer, intent(in) :: points

! Level for gw launch - in principle this could be launchlev(i,j)
integer, intent(in) :: launchlev

! Starting level
integer, intent(in) :: k_start

! Top model level
integer, intent(in) :: k_end

! Factor enhancement for invariant global wave launch amplitude
real(kind=real_usprec), intent(in)    :: ussp_launch_factor

! omega is Angular speed of Earth's rotation, 2 * omega
real(kind=real_usprec), intent(in)    :: two_omega

! Scale conversion factor from convective rain to GW flux
real(kind=real_usprec), intent(in)    :: cgw_scale_factor

! Characteristic (spectrum peak) wavelength (m)
real(kind=real_usprec), intent(in)    :: wavelstar

! Latitudes on grid
real(kind=real_usprec), intent(in)    :: sin_theta_latitude(points)

! Rho (kg m^-3)
real(kind=real_usprec), intent(in)    :: rho_th(points,k_start:k_end)

! Buoyancy [Brunt Vaisala] Frequency (rad s^-1)
real(kind=real_usprec), intent(in)    :: nbv(points,k_start:k_end)

! Component of wind in phi_jdir direction (m s^-1)
real(kind=real_usprec), intent(in)    :: udotk(points,k_start:k_end,           &
                                                                  idir)

! Total Precipitation for CGW source
real(kind=real_usprec), intent(in)    :: totalppn(points)

! Pseudomomentum flux integrated over azimuthal sector (kg m^-1 s^-2)
real(kind=real_usprec), intent(out)   :: fptot(points,k_start:k_end,           &
                                                                   idir)

! Flag to indicate CGW flux calculation should be carried out
logical, intent(in) :: L_add_cgw

! ----------------------------------------------------------------------+-------
! Local parameters
! ----------------------------------------------------------------------+-------

! Maximum number of iterations of Newton Raphson do (While) loop
integer, parameter :: maxwhile       = 9

real(kind=real_usprec), parameter :: zero = 0.0_real_usprec
real(kind=real_usprec), parameter :: one = 1.0_real_usprec

! Parameter beta in the launch spectrum total energy equation
real(kind=real_usprec), parameter :: beta_e0    = 1.0227987125e-1_real_usprec

! Parameter p in B_0(p) for launch spectrum intrinsic frequency
! NOTE: This parameter determines the intrinsic frequency spectrum
!       shape and hence the integral form in 4.1, which is strictly
!       valid only for p > 1. !!if contemplating changes BE WARNED!!
real(kind=real_usprec), parameter :: psat       = 5.0_real_usprec /            &
                                                  3.0_real_usprec

! Psat - 1
real(kind=real_usprec), parameter :: psatm1     = psat - 1.0_real_usprec

! 2 - Psat
real(kind=real_usprec), parameter :: twompsat   = 2.0_real_usprec - psat

! Reciprocal of max wavelength at launch (/m)
real(kind=real_usprec), parameter :: lminl      = 1.0_real_usprec/20000

! Power s of vertical wavenumber spectrum A_0(s,t) at low m
real(kind=real_usprec), parameter :: ss         = 1.0_real_usprec

! s + 1
real(kind=real_usprec), parameter :: ssp1       = ss + 1.0_real_usprec

! Power t=t_sat of vertical wavenumber spectrum at large m due to
! saturation by gravity wave breaking (and shape of chopping fn)
real(kind=real_usprec), parameter :: tt         = 3.0_real_usprec

! t - 1, t - 2, 1 / (t-2), (t-3) / (t-2), 2 - t
real(kind=real_usprec), parameter :: ttm1       = tt - 1.0_real_usprec
real(kind=real_usprec), parameter :: ttm2       = tt - 2.0_real_usprec
real(kind=real_usprec), parameter :: rttm2      = 1.0_real_usprec / ttm2
real(kind=real_usprec), parameter :: ttrat      = (tt - 3.0_real_usprec) * rttm2
real(kind=real_usprec), parameter :: twomtt     = 2.0_real_usprec - tt

! s + t, 1 / (s+t)
real(kind=real_usprec), parameter :: ssptt      = ss + tt
real(kind=real_usprec), parameter :: rssptt     = 1.0_real_usprec / ssptt

! Weight for (n+1)th guess at mNlX in iteration solution
real(kind=real_usprec), parameter :: mweight    = 0.8_real_usprec

! Strength coefficient constant for Launch spectrum (CCL / A0)
real(kind=real_usprec), parameter :: ccl0       = 3.41910625e-9_real_usprec

! Scale conversion factor from convective rain to GW flux (mPa)
real(kind=real_usprec), parameter :: cor2flux   = 7.20507e-04_real_usprec

! ----------------------------------------------------------------------+-------
! Security parameters
! ----------------------------------------------------------------------+-------

! Minimum allowed non-zero value of precipitation (0.1 mm / day)
real(kind=real_usprec), parameter :: rppnmin   = 1.0_real_usprec /             &
                                                 1.1574e-06_real_usprec

! Minimum allowed non-zero value of Curvature Coefficient A
real(kind=real_usprec), parameter ::  asecp     =  1.0e-20_real_usprec
real(kind=real_usprec), parameter ::  asecn     = -(asecp)

! ----------------------------------------------------------------------+-------
! Local Constants (a) Physical
! ----------------------------------------------------------------------+-------

! Normalised minimum vertical wavenumber at launch
real(kind=real_usprec)            ::  mnlmin

! Wavenumber at peak in spectrum (/m)
real(kind=real_usprec)            ::  mstar

! Reciprocal of mstar (m) and mstar^2 (m^2)
real(kind=real_usprec)            ::  rmstar
real(kind=real_usprec)            ::  rmstarsq

! Equatorial planetary vorticity gradient parameter B_eq (/m /s )
real(kind=real_usprec)            :: beta_eq_rmstar

! Normalisation factor for the launch spectrum vertical wavenumber A0/(s+1)(t-1)
real(kind=real_usprec)            :: a0_r_sp1tm1

! Normalisation factor for the launch spectrum vertical wavenumber -A0/(t-1)
real(kind=real_usprec)            :: a0_r_1mt

! Integral segment [(s+1) * mNLmin**(1-t)]
real(kind=real_usprec)            :: tail_chop2b

! Integral segment [(t-1) * mNLmin**(s+1)]
real(kind=real_usprec)            :: head_chop2a

! Isotropic GW Launch Flux = mStar^(-2) * [CCL0 * ussp_launch_factor]
real(kind=real_usprec)            :: glob_launch_flux

! Scale conversion factor from convective rain to GW flux
! [cor2flux * cgw_scale_factor]
real(kind=real_usprec)            :: cgw_convert_factor

! ----------------------------------------------------------------------+-------
! Local variables (scalars) used in GW_USSP
! Some effectively expanded to workspace (using vector registers)
! ----------------------------------------------------------------------+-------

! Points index
integer         :: i

! Level index
integer         :: k

! Index values for chop case loops
integer         :: nni

! Azimuthal direction index
integer         :: jdir

! Counter for while loop
integer         :: jwhile


! Number of spectra with low m intersect
integer         :: nchop2

! I location of chop type points
integer         :: indexi(points)

! Number of spiral A solution points
! integer         :: nspira

! ----------------------------------------------------------------------+-------

! Azimuthal sector for launch spectrum integral Delta Phi / 2
real(kind=real_usprec)            :: ddphir2

! Inertial frequency at current latitude (rad s^-1)
real(kind=real_usprec)            :: f_f

! High wavenumber intersection point
real(kind=real_usprec)            :: mnly

! Wavenumber where Doppler transformed spectrum is reduced to zero
real(kind=real_usprec)            :: mkill

! omega_min(launch) / N (k)
real(kind=real_usprec)            :: ominrnbv

! Constant component of saturation spectrum
real(kind=real_usprec)            :: ccs0rmstarsq

! Minimum range value of function f
real(kind=real_usprec)            :: fminus

! Intermediate  value of function f
real(kind=real_usprec)            :: fterm

! Maximum range value of function f
real(kind=real_usprec)            :: fplus

! Minimum range value of function g
real(kind=real_usprec)            :: gminus

! Intermediate  value of function g
real(kind=real_usprec)            :: gterm

! Maximum range value of function g
real(kind=real_usprec)            :: gplus

! Indicates spectra with low m intersect
logical         :: L_chop2
! ----------------------------------------------------------------------+-------
! Local variables (dynamic arrays) used in GW_USSP
! ----------------------------------------------------------------------+-------

! Delta U=udotk(launch)-udotk(k)
real(kind=real_usprec)            :: ddu_a

! Record of total flux
real(kind=real_usprec)            :: fpfac(points)

! Coefficient A in intersect point equation
real(kind=real_usprec)            :: acoeff

! Term (A / B) in intersect point equation
real(kind=real_usprec)            :: curvature(points)

! Coefficient B in intersect point equation
real(kind=real_usprec)            :: attenuation(points)

! Chop function B*[1 + (A/B)]^(t-2)
real(kind=real_usprec)            :: intercept1(points)

! Chop function B*[1 + (A/B)*mNlmin]^(t-2)
real(kind=real_usprec)            :: mintercept(points)

! Starting value of vertical wavenumber for crossing point search
real(kind=real_usprec)            :: mguess_a(points)

! Intersect mNlX estimates
real(kind=real_usprec)            :: mnlx(points,0:maxwhile)

! Compressed attenuation array
real(kind=real_usprec)            :: atte_c(points)

! Compressed curvature array
real(kind=real_usprec)            :: curv_c(points)

! Weighting of n term in iter
real(kind=real_usprec)            :: wgtn(points)

! Either f_f or the equatorial minimum frequency, whichever is less (rad s^-1)
real(kind=real_usprec)            :: omin(points)

! [CGW Flux]_klaunch
real(kind=real_usprec)            :: cgwll(points)

! [Rho . Cl]_klaunch
real(kind=real_usprec)            :: rhocl(points)

! [Rho(z) . Csat(z) / m*^2]_k
real(kind=real_usprec)            :: fsatk_scale(points, launchlev+1:k_end)

! Indicator for direction of spiral solution
logical         :: L_ftheng(points)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='GW_USSP_CORE'

!  End of Header

! ==Main Block==--------------------------------------------------------+-------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------
! Local Constants (a) Physical
! ----------------------------
mstar    = 2.0_real_usprec * pi / wavelstar
rmstar   = one / mstar
rmstarsq = rmstar * rmstar
beta_eq_rmstar = 2.3e-11_real_usprec * rmstar

mnlmin   = wavelstar * lminl
fminus = mnlmin**ssptt
tail_chop2b = ssp1 / (mnlmin**ttm1)
head_chop2a = ttm1 * (mnlmin**ssp1)
a0_r_sp1tm1 = one / ( ss + tt - head_chop2a )
a0_r_1mt    = (-(ssp1) ) * a0_r_sp1tm1
ddphir2 = pi / idir
ccs0rmstarsq = rmstarsq * beta_e0 * sin(ddphir2) * psatm1 /(pi*twompsat)

! ----------------------------------------------------------------------+-------
! Core calculations on Physics Grid now within separate subroutine
! ----------------------------------------------------------------------+-------

! L_add_cgw = F : calculate standard USSP isotropic GW launch flux
! L_add_cgw = T : calculate variable CGW launch flux

if (L_add_cgw) then
  glob_launch_flux = zero
  cgw_convert_factor = cor2flux * cgw_scale_factor
else
  glob_launch_flux = rmstarsq * ccl0 * ussp_launch_factor
  cgw_convert_factor = zero
end if

! ----------------------------------------------------------------------+-------
! 3.0 Initialize gravity wave spectrum variables
! ----------------------------------------------------------------------+-------


! ----------------------------------------------------------------------+-------
! 3.1 Compute minimum intrinsic frequency Omin from inertial frequency
!     squared f_f. See UMDP 34, eqn. (A.7).
!-----------------------------------------------------------------------+-------

do i=1, points
  f_f = (two_omega * sin_theta_latitude(i))**2
  omin(i) = nbv(i,launchlev) * beta_eq_rmstar
  omin(i) = max( omin(i), f_f )
  omin(i) = sqrt(omin(i))

  ! ----------------------------------------------------------------------+-----
  !   Convection source of flux at Launch
  ! ----------------------------------------------------------------------+-----
  cgwll(i) = totalppn(i) * rppnmin
  cgwll(i) = max( cgwll(i), one )
  cgwll(i) = sqrt(cgwll(i)) * cgw_convert_factor
end do

! ----------------------------------------------------------------------+-------
! 3.2 Compute [rho(z) . C(z)]_k / m*^2 scaling factor for quasi-saturated
!     spectrum (as per ... UMDP 34: 1.15,1.16)
! ----------------------------------------------------------------------+-------

! if (abs(psatm1) >= 0.1) then
! For current setting of parameter psat this test is always true

do k=launchlev+1,k_end
  do i=1, points
    ominrnbv = omin(i) / nbv(i,k)
    fsatk_scale(i,k) = rho_th(i,k) * ccs0rmstarsq *                            &
     (nbv(i,k))**2 * (ominrnbv**psatm1) *                                      &
     (one - (ominrnbv**twompsat)) / (one - (ominrnbv**psatm1))
  end do
end do

! else
! Require a different functional form for normalisation factor B0
!   BBS = 1.0 / alog(nbv(i,k) / omin(i,j))
! end if


! ----------------------------------------------------------------------+-------
! 4.0 Calculations carried out for each azimuthal direction
! ----------------------------------------------------------------------+-------


do i=1, points
  ! ----------------------------------------------------------------------+-----
  !     Initialize gravity wave spectrum variables for the launch level
  ! ----------------------------------------------------------------------+-----

  ! ----------------------------------------------------------------------+-----
  !     Globally invariant value of Total vertical Flux of horizontal
  !     pseudomomentum at Launch level, analytic integral under curve
  !     = rho_l * C_l / m*^2 ... UMDP 34: 1.14
  ! ----------------------------------------------------------------------+-----
  rhocl(i) = rho_th(i,launchlev) * (glob_launch_flux + cgwll(i))
end do


do jdir=1,idir
  do k=k_start,launchlev-1
    do i=1, points
      ! ----------------------------------------------------------------------+-
      !       Set total flux of horizontal pseudomomentum at bottom
      ! ----------------------------------------------------------------------+-
      fptot(i,k,jdir) = zero
    end do
  end do

  do k=launchlev,k_end
    do i=1, points
      ! ----------------------------------------------------------------------+-
      !     Note: the total flux at levels above is initialised to the sum of
      !     all sources propagated up through the column.
      ! ----------------------------------------------------------------------+-
      if (k == launchlev) then
        fptot(i,k,jdir) = rhocl(i)
      else
        fptot(i,k,jdir) = zero
      end if
      fptot(i,k,jdir) = fptot(i,k-1,jdir) + fptot(i,k,jdir)
    end do
  end do
end do

! ----------------------------------------------------------------------+-------
!     Loop over directions and levels and calculate horizontal component
!     of the vertical flux of pseudomomentum for each azimuthal
!     direction and for each altitude
! ----------------------------------------------------------------------+-------
IDir_do2a: do jdir=1, idir
  Levels_do92: do k=launchlev+1,k_end
    L_chop2 = .false.
    Points_do92: do i=1, points
      ! ----------------------------------------------------------------------+-
      !       Initialise MGUESS (start point for iterative searches if needed)
      ! ----------------------------------------------------------------------+-
      mguess_a(i) = zero
      Fptot_if1: if (fptot(i,k-1,jdir) >  zero) then
        ! ----------------------------------------------------------------------
        ! 4.2       Calculate variables that define the Chop Type Cases.
        ! ----------------------------------------------------------------------
        ddu_a = udotk(i,launchlev,jdir) -                                      &
                udotk(i,k,jdir)
        ! ----------------------------------------------------------------------
        !        UMDP 34: 1.23 coefficient B
        !        Using ratio of flux scalings rather than densities is more
        !        robust because total flux imported from other schemes may
        !        have different relationship to launch density than 1.14
        !        Note: the total flux is initialised to the sum of all sources
        !        propagated up through the column.
        ! ----------------------------------------------------------------------
        attenuation(i) =                                                       &
           ( fsatk_scale(i,k) / fptot(i,k,jdir) ) *                            &
           ( nbv(i,launchlev) / nbv(i,k) )**ttm1
        ! ----------------------------------------------------------------------
        !           UMDP 34: 1.22 coefficient A = (A/B) * B
        ! ----------------------------------------------------------------------
        curvature(i) = ddu_a * mstar /  nbv(i,launchlev)

        acoeff = curvature(i) * attenuation(i)

        mintercept(i) = attenuation(i) *                                       &
                ( (mnlmin * curvature(i)) + one )**ttm2

        Curv_if1: if (curvature(i) <  asecn) then
          ! --------------------------------------------------------------------
          !           Negative Doppler Shift : factor will hit zero (kill point)
          ! --------------------------------------------------------------------
          mkill = one / abs(curvature(i))

          Mkill_if1: if (mkill <= mnlmin ) then
            ! ------------------------------------------------------------------
            !             Chop Type IV : No flux propagates
            ! ------------------------------------------------------------------
            fptot(i,k,jdir) = zero
          else
            if (mkill >  one) then
              intercept1(i) = attenuation(i) *                                 &
                             ( one + curvature(i) )**ttm2
            else
              ! ----------------------------------------------------------------
              !          Doppler factor minimum (kill point) situated below
              !          mstar in the low-m part of the launch spectrum
              ! ----------------------------------------------------------------
              intercept1(i) = zero
            end if

            Lowend_if1: if (intercept1(i) >= one) then
              ! ----------------------------------------------------------------
              !    Chop Type I: Intersection in high wavenumber part only
              ! ----------------------------------------------------------------
              mnly = ( attenuation(i)**ttrat -                                 &
                       attenuation(i)) / acoeff

              fptot(i,k,jdir) = fptot(i,k,jdir) *                              &
           (one - (a0_r_1mt * curvature(i) * mnly**twomtt))
            else
              if (mintercept(i) <= fminus) then
                ! --------------------------------------------------------------
                !      Chop Type IIb: Low wavenumber intersect only below min
                ! --------------------------------------------------------------
                fptot(i,k,jdir) = fptot(i,k,jdir) *                            &
                 a0_r_sp1tm1 * tail_chop2b * mintercept(i)                     &
                       * ( (mnlmin * curvature(i)) + one )
              else
                ! --------------------------------------------------------------
                !      Chop Type IIa: Low wavenumber intersect only
                ! --------------------------------------------------------------
                L_chop2 = .true.
                mguess_a(i) = min(mkill, one)
                fpfac(i) = fptot(i,k,jdir) * a0_r_sp1tm1
                fptot(i,k,jdir) = zero
              end if
            end if  Lowend_if1

          end if  Mkill_if1

        else if (curvature(i) >  asecp) then
          ! --------------------------------------------------------------------
          !           Positive Doppler Shift : non-zero factor (no kill point)
          ! --------------------------------------------------------------------
          intercept1(i) = attenuation(i) *                                     &
                             ( one + curvature(i) )**ttm2

          Chop3_if1: if (intercept1(i) <  one) then
            ! ------------------------------------------------------------------
            !             Chop Type III: Intersection in both wavenumber parts
            ! ------------------------------------------------------------------
            fpfac(i) = fptot(i,k,jdir) * a0_r_sp1tm1

            ! ------------------------------------------------------------------
            !               First find intersect in high wavenumber part
            !               UMDP 34: 1.25
            ! ------------------------------------------------------------------
            mnly = ( attenuation(i)**ttrat -                                   &
                     attenuation(i)) / acoeff

            fptot(i,k,jdir) = fptot(i,k,jdir) * a0_r_1mt *                     &
                            curvature(i) * mnly**twomtt

            ! ------------------------------------------------------------------
            !               Then find intersect in low wavenumber part to reckon
            !               its flux contribution for addition when available
            ! ------------------------------------------------------------------
            if (mintercept(i) <= fminus) then
              ! ----------------------------------------------------------------
              !             Chop Type IIIb: Low wavenumber intersect below min
              ! ----------------------------------------------------------------
              fptot(i,k,jdir) = fptot(i,k,jdir) +                              &
                                ( fpfac(i) * tail_chop2b *                     &
                             mintercept(i) *                                   &
                  ( (mnlmin * curvature(i)) + one ) )
            else
              ! ----------------------------------------------------------------
              !               Chop Type IIIa: Low wavenumber intersect
              ! ----------------------------------------------------------------
              L_chop2 = .true.
              mguess_a(i) = one
            end if

            ! ------------------------------------------------------------------
            !             else Chop Type 0: No intersection (spectrum unaltered)
            ! ------------------------------------------------------------------
          end if  Chop3_if1
        else
          ! --------------------------------------------------------------------
          !           Negligible Doppler shift
          ! --------------------------------------------------------------------
          !             Strictly this is analytic solution mNLX.  UMDP 34: 1.27
          mnly = attenuation(i)**rssptt
          if (mnly <= mnlmin) then
            ! ------------------------------------------------------------------
            !            Chop Type IIb: Low wavenumber intersect only below min
            ! ------------------------------------------------------------------
            fptot(i,k,jdir) = fptot(i,k,jdir) * a0_r_sp1tm1 *                  &
                             tail_chop2b * attenuation(i)
          else
            if (mnly <  one)  fptot(i,k,jdir) =                                &
            ! ------------------------------------------------------------------
            !            Chop Type IIc: Low wavenumber intersect only (analytic)
            ! ------------------------------------------------------------------
                          fptot(i,k,jdir) * a0_r_sp1tm1 *                      &
                                   ( (ssptt * (mnly**ssp1)) - head_chop2a )
            ! ------------------------------------------------------------------
            !            else Chop Type 0: No intersection (spectrum unaltered)
            ! ------------------------------------------------------------------
          end if
        end if  Curv_if1

      end if  Fptot_if1
    end do  Points_do92

    Lchop2_if1: if (L_chop2) then
      ! ----------------------------------------------------------------------+-
      !   Process low wavenumber contribution: evaluate intersect mNX
      ! ----------------------------------------------------------------------+-
      nchop2 = 0

      do i=1, points
        if (mguess_a(i) >  zero) then
          nchop2 = nchop2 + 1
          indexi(nchop2)  = i
        end if
      end do

      Nchop2_do1: do i=1,nchop2
        ! ----------------------------------------------------------------------
        !     Chop Type IIa : / Full solution required for mNlX
        ! or  Chop Type IIIa: \
        ! ----------------------------------------------------------------------
        nni  = int(indexi(i))

        fplus  = mguess_a(nni)**ssptt
        gplus  = intercept1(nni)
        !       fminus = mnlmin**ssptt    Defined as a constant
        gminus = mintercept(nni)
        atte_c(i) = attenuation(nni)
        curv_c(i) = curvature(nni)

        fterm = ( ((fminus / atte_c(i))**rttm2) - one ) / curv_c(i)
        gterm = gminus**rssptt
        L_ftheng(i) = .false.

        Curv_if2: if (curvature(nni) >  asecp) then
          ! --------------------------------------------------------------------
          !       Positive Doppler Shift
          ! --------------------------------------------------------------------
          wgtn(i) = zero
        else
          ! --------------------------------------------------------------------
          !       Negative Doppler Shift
          ! --------------------------------------------------------------------
          wgtn(i) = one - mweight

          if (fplus <= gminus  .and.  gplus >  fminus) then
            fterm = (((fplus / atte_c(i))**rttm2) - one)/ curv_c(i)
            gterm = gplus**rssptt
            L_ftheng(i) = (gterm <  fterm)

          else if (fplus >  gminus  .and.  gplus <= fminus) then
            L_ftheng(i) = (gterm >= fterm)

          else if (fplus <= gminus  .and.  gplus <= fminus) then
            L_ftheng(i) = .true.

            !         else Use default settings
          end if
        end if  Curv_if2

        if (L_ftheng(i)) then
          !         nspira = nspira + 1
          mnlx(i,0) = fterm
        else
          mnlx(i,0) = gterm
        end if
      end do  Nchop2_do1

      do jwhile=0,maxwhile-1
#if defined (CRAYFTN_VERSION)
        ! This directive specifies that it is safe to execute all memory
        ! references and arithmetic operations within all conditional branches
        ! of a loop.

!DIR$ SAFE_CONDITIONAL
#endif
        do i=1,nchop2

          if (L_ftheng(i)) then
            ! ------------------------------------------------------------------
            !         Obtain m_n+1 from g_n+1  = f_n (m_n)
            ! ------------------------------------------------------------------
            mnlx(i,jwhile+1) = (                                               &
                 (((mnlx(i,jwhile)**ssptt) / atte_c(i))**rttm2) - one )        &
                 / curv_c(i)
          else
            ! ------------------------------------------------------------------
            !         Obtain m_n+1 from f_n+1  = g_n (m_n)
            ! ------------------------------------------------------------------
            mnlx(i,jwhile+1) = ( (atte_c(i) *                                  &
                ((one + (curv_c(i) * mnlx(i,jwhile)))**ttm2))**rssptt )
          end if

          mnlx(i,jwhile+1) = ((one - wgtn(i)) * mnlx(i,jwhile+1)) +            &
                                    (wgtn(i)  * mnlx(i,jwhile))

        end do
      end do

!DIR$ IVDEP
#if defined (CRAYFTN_VERSION)
!DIR$ VECTOR ALWAYS
#endif
      do i=1,nchop2
        nni  = int(indexi(i))

        fptot(nni,k,jdir) = fptot(nni,k,jdir) + (                              &
        fpfac(nni) * ( ((mnlx(i,maxwhile)**ssp1) * ( ssptt +                   &
        (ssp1 * mnlx(i,maxwhile) * curv_c(i)) ) ) - head_chop2a ) )

      end do
    end if  Lchop2_if1

    do i=1, points
      !-----------------------------------------------------------------------+-
      !       Now correct pseudomomentum flux in the evolved spectrum if the
      !       new value is non-physical (pseudomomentum flux cannot increase
      !       with altitude)
      !-----------------------------------------------------------------------+-
      fptot(i,k,jdir) = min(fptot(i,k,jdir),fptot(i,k-1,jdir))

    end do

  end do  Levels_do92
end do  IDir_do2a

! ----------------------------------------------------------------------+-------
! 4.5   Apply Opaque Upper Boundary to set fluxes to zero at top
! ----------------------------------------------------------------------+-------
do jdir=1,idir
  do i=1, points
    fptot(i,k_end,jdir) =  zero
  end do
end do

! ----------------------------------------------------------------------+-------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine gw_ussp_core

end module gw_ussp_core_mod
