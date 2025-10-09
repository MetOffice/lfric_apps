! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!  data module for switches/options concerned with the cloud scheme.

  ! Module holding constants for the PC2 cloud scheme
  !   A description of what each switch or number refers to is provided
  !   with the namelist
  !
  !   Any routine wishing to use these options may do so with the 'USE'
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !

module pc2_constants_mod

use um_types, only: real_umphys

implicit none

! Which cloud scheme?
integer, parameter :: i_cld_off = 0
! No cloud scheme, just use whatever happens in the microphysics
integer, parameter :: i_cld_smith = 1
! The Smith (1990) diagnostics parametrization
integer, parameter :: i_cld_pc2 = 2
! The Wilson et al (2008) prognostic parametrization
integer, parameter :: i_cld_bimodal = 3
! The Van Weverberg et al (2020) diagnostic parametrization

!=================================================================
! PC2 Cloud Scheme Terms
!=================================================================

! Constants migrated from c_lspmic.h
! --------------------------------------------------------------------
real(kind=real_umphys), parameter :: wind_shear_factor = 1.5e-4
! Parameter that governs the rate of spread of ice cloud fraction
! due to windshear

real(kind=real_umphys), parameter :: cloud_rounding_tol = 1.0e-12
! Tolerance to use when comparing cloud amounts, to allow for
! possible effects of rounding error.

! Constants migrated from ni_imp_ctl.F90
! --------------------------------------------------------------------
real(kind=real_umphys), parameter :: ls_bl0 = 1.0e-4
! Specified in-plume value of ice content when there is not enough
! information to calculate this parameter another way (kg kg-1)

! Erosion options
! --------------------------------------------------------------------
integer, parameter :: pc2eros_exp_rh = 1
! Original method described in Wilson et al. (2008). Uses the
! notion of the relative rate of narrowing of the moisture PDF.
! This rate is related pragmatically to RH, using an ad hoc
! exponential. Change to qcl and CFL are both calculated as a
! result of narrowing the PDF.

! Method 2 has been retired.

integer, parameter :: pc2eros_hybrid_sidesonly = 3
! As pc2_erosion_hybrid_method_all_faces but surface area of
! exposed cloud is calculated considering the lateral sides only.

! Erosion numerical method options
! --------------------------------------------------------------------
integer, parameter :: i_pc2_erosion_explicit = 1  ! Explicit forwards-in-time
integer, parameter :: i_pc2_erosion_implicit = 2  ! Implicit wrt cloud-fraction
integer, parameter :: i_pc2_erosion_analytic = 3  ! Analytical solution

! Options for PC2 initiation
! --------------------------------------------------------------------
integer, parameter :: pc2init_smith = 1
! Original method described in Wilson et al. (2008). Uses a
! triangular moisture PDF with a width determined by the critical
! relative humidity to initiate cloud from clear sky or remove cloud
! from overcast conditions.
integer, parameter :: pc2init_bimodal = 2
! Use the bimodal diagnostic cloud scheme to initiation cloud from
! clear sky or remove cloud from overcast conditions. This scheme
! reconstructs a bimodal moisture PDF within entrainment zones near
! sharp inversions

! Options for PC2 initiation logic
! --------------------------------------------------------------------
integer, parameter :: pc2init_logic_original = 1
! Original logic; initiation is only allowed to occur where the existing
! cloud-fraction is very close to 0 or 1, and is not allowed below the
! diagnosed cloud-base in well-mixed layers.  Includes various arbitrary
! thresholds on cloud-fraction, temperature, whether cumulus is diagnosed,
! whether RH is increasing in an Eulerian sense etc.
integer, parameter :: pc2init_logic_simplified = 2
! Slightly simplified version of the above; removes the dependence
! on temperature and the cumulus flag.
integer, parameter :: pc2init_logic_smooth = 3
! Much simpler logic that should yield smoother behaviour;
! The diagnostic cloud scheme calculations are performed at any
! grid-point where they are expected to yield nonzero liquid water
! content.  The diagnosed qcl is then taken as a minimum limit applied
! to PC2's prognostic qcl.

! Options for updating ice cloud fraction due to ice cloud fraction
! falling in from layer above.
! --------------------------------------------------------------------
integer, parameter :: original_but_wrong = 1
! Original code used in Wilson et al (2008a,b). It incorrectly
! calculates the fall velocity of ice into the current layer from
! the layer above, using the fall velocity in the current layer
! (rather than the velocity in layer above). This allows ice to
! fall into layer when there is no ice above! This option also uses
! a hard-wired, globally constant wind_shear_factor, rather than
! the shear derived from the model winds.
integer, parameter :: ignore_shear = 0
! Correctly determines the fall velocity into the layer using the
! fall velocity of ice in the layer above.
integer, parameter :: real_shear = 2
! Correctly determines the fall velocity into the layer using the
! fall velocity of ice in the layer above. Adjusts the overhang by
! using the true wind-shear calculated from the wind and also
! translates the lateral displacement of the falling ice cloud due
! to shear into a cloud fraction increment by considering the size
! of the grid-box.
! --------------------------------------------------------------------
! Options for representing forced cumulus:
!   0=off, 1=on=capping convective boundary layers only
integer, parameter :: cbl_and_cu = 2  ! also used at cumulus cloud base
integer, parameter :: forced_cu_cca = 3 ! As cbl_and_cu, but include
                                        ! forced cumulus cloud in the
                                        ! diagnosed convective cloud, instead
                                        ! of the prognostic cloud fields.

! Options for critical relative humidity parametrization
! --------------------------------------------------------------------
integer, parameter :: rhcpt_off = 0
! No parametrization of rhcrit, it is specified as a list of values
! in the namelist
integer, parameter :: rhcpt_tke_based = 2
! New method based on parametrizing the sub-grid variance of T and q
! from the gradients, eddy diffusivities and TKE diagnosed in the
! boundary layer scheme

! Options for the cloud area parametrization
! --------------------------------------------------------------------
integer, parameter :: acf_off = 0
! No cloud area parametrization, the area cloud fraction equals the
! bulk cloud fraction
integer, parameter :: acf_cusack = 1
! Cusack method based on temperature and moisture gradients, modified
! for use with PC2 and described in Boutle & Morcrette (2010, ASL)
integer, parameter :: acf_brooks = 2
! Brooks method based on an empirical relationship between cloud area
! and volume derived from radar and lidar

! Options for PC2 homogeneous forcing calculation of the PDF amplitude
! at the saturation boundary g(-qc)
! --------------------------------------------------------------------
integer, parameter :: i_pc2_homog_g_cf = 1
! Blend solutions from liquid cloud water and saturation defecit
! weighted by a function of cloud-fraction.
integer, parameter :: i_pc2_homog_g_width = 2
! Blend solutions from liquid cloud water and saturation defecit
! weighted by their respective PDF-widths.

! Constants migrated from c_pc2pdf.h

      ! Number of iterations in the initiation
integer, parameter :: init_iterations = 10

! Tolerance of critical relative humidity for initiation of cloud
real(kind=real_umphys), parameter :: rhcrit_tol = 0.01

real(kind=real_umphys), parameter :: bm_tiny = 1.0e-10
! Small minimum value used at various places in the bimodal cloud scheme

real(kind=real_umphys), parameter :: bm_negative_init = -999.0
! Negative default value used at various places in the bimodal cloud scheme

! Power that is used to weight the two values of G when they are
! merged together
real(kind=real_umphys), parameter :: pdf_merge_power = 0.5

! Power that describes the way G varies with s near the boundary
! of the cloud probability density function. For a "top-hat" it
! is equal to 0, for a "triangular" distribution it is equal to 1.
real(kind=real_umphys), parameter :: pdf_power = 0.0

! Parameters that govern the turbulent decrease in width of the
! PDF.  (dbs/dt)/bs = (DBSDTBS_TURB_0 + dQc/dt DBSDTBS_TURB_1)
!                 * exp( - dbsdtbs_exp Q_c / (a_L qsat(T_L)))
! dbsdtbs_turb_0 is set in the namelist run_cloud PC2 options
real(kind=real_umphys), parameter :: dbsdtbs_turb_1 = 0.0
real(kind=real_umphys), parameter :: dbsdtbs_conv   = 0.0
real(kind=real_umphys), parameter :: dbsdtbs_exp    = 10.05

! Constant migrated from pc2_checks
real(kind=real_umphys), parameter :: one_over_qcf0 = 1.0e4
! One_over_qcf0 is reciprocal of a reasonable in-cloud ice content.
real(kind=real_umphys), parameter :: min_in_cloud_qcf = 1.0e-6
! Minium in-cloud ice water content ensured by reducing CFF
! Ensure value below is 1.0/value above
real(kind=real_umphys), parameter :: one_over_min_in_cloud_qcf = 1.0e6
! One over MIN_IN_CLOUD_QCF.
real(kind=real_umphys), parameter :: max_in_cloud_qcf = 2.0e-3
! Max allowed in-cloud ice water content after inhomogeneous BL forcing.
real(kind=real_umphys), parameter :: condensate_limit = 1.0e-10
! Minimum value of condensate
real(kind=real_umphys), parameter :: wcgrow           = 5.0e-4

end module pc2_constants_mod
