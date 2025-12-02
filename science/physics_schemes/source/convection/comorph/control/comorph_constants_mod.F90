! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module comorph_constants_mod

implicit none

save

!################################################################
! Constants for the CoMorph convection scheme
!################################################################

! Flag indicating whether the variables contained in this
! module have been set
logical :: l_init_constants = .false.


!----------------------------------------------------------------
! Hardwired integer indicators
!----------------------------------------------------------------

! Indicators for each hydrometeor species
integer, parameter :: i_cl = 1           ! Liquid cloud
integer, parameter :: i_rain = 2         ! Rain
integer, parameter :: i_cf = 3           ! Ice cloud
integer, parameter :: i_snow = 4         ! Snow
integer, parameter :: i_graup = 5        ! Graupel
integer, parameter :: i_frogs = 10       ! Frogs swept up by gust
integer, parameter :: i_blueberries = 11 ! Blueberries

! Indicators for different sub-grid fractions
integer, parameter :: i_sg_homog = 0     ! Homogeneous over whole grid-box
integer, parameter :: i_sg_frac_liq = 1  ! Only in liquid cloud fraction
integer, parameter :: i_sg_frac_ice = 2  ! Only in ice cloud fraction
integer, parameter :: i_sg_frac_prec = 3 ! Only in the precipitation fraction

! Precision, range and kind for 64 bit real
integer, parameter :: prec64  = 15
integer, parameter :: range64 = 307
integer, parameter :: real64_kind = selected_real_kind(prec64,range64)
! Precision, range and kind for 32 bit real
integer, parameter :: prec32  = 6
integer, parameter :: range32 = 37
integer, parameter :: real32_kind = selected_real_kind(prec32,range32)

! Precision switches: (either 32 or 64 bit):
!   For the convection scheme's internal variables
integer, parameter :: real_cvprec = real32_kind
!   For variables passed in or out to the host model
integer, parameter :: real_hmprec = real64_kind

! Note all real constants below are declared with
! the convection scheme's internal precision, and need
! to be converted back to host-model precision if used
! directly with the full fields passed into CoMorph

! Length of character strings storing diagnostic names etc
integer, parameter :: name_length = 100



!----------------------------------------------------------------
! Magic numbers stored in native precision
!----------------------------------------------------------------

! pi = circle circumference / circle diameter
real(kind=real_cvprec), parameter :: pi = 3.14159_real_cvprec

real(kind=real_cvprec), parameter :: zero  = 0.0_real_cvprec
real(kind=real_cvprec), parameter :: one   = 1.0_real_cvprec
real(kind=real_cvprec), parameter :: two   = 2.0_real_cvprec
real(kind=real_cvprec), parameter :: three = 3.0_real_cvprec
real(kind=real_cvprec), parameter :: four  = 4.0_real_cvprec
real(kind=real_cvprec), parameter :: six   = 6.0_real_cvprec

real(kind=real_cvprec), parameter :: half = 0.5_real_cvprec
real(kind=real_cvprec), parameter :: quarter = 0.25_real_cvprec
real(kind=real_cvprec), parameter :: third = one/three
real(kind=real_cvprec), parameter :: four_thirds = four/three
real(kind=real_cvprec), parameter :: three_over_eight = three/8.0_real_cvprec

! Smallest possible floating point number
real(kind=real_cvprec), parameter :: min_float = tiny(zero)
! Square-root of the above (a safe bet as a min threshold for safe division)
real(kind=real_cvprec), parameter :: sqrt_min_float = sqrt(min_float)

! Smallest representable increment between 1 and 1+delta
real(kind=real_cvprec), parameter :: min_delta = epsilon(zero)
! Square-root of the above
real(kind=real_cvprec), parameter :: sqrt_min_delta = sqrt(min_delta)


!---------------------------------------------------------------
! Things set at run-time by the host-model
!---------------------------------------------------------------
! (where given, initialisation values are defaults used in the
!  stand-alone unit tests).


! ARRAY DIMENSIONS

! Full 2-D array dimensions
integer :: nx_full = 0
integer :: ny_full = 0

! Lowest and highest model-levels where convection is allowed.
integer :: k_bot_conv = 0
integer :: k_top_conv = 0

! Note: convection will operate on the points
! ( 1:nx_full, 1:ny_full, k_bot_conv:k_top_conv ),
! even if the input environment profile arrays have different
! dimensions to this (which may well be correct, since the
! input arrays may have halos and/or boundary-condition points
! where we don't want to call convection).
! The input and output array arguments are all passed in/out
! via derived type structures containing pointers to the actual
! 3D arrays.  This allows the convection code to not have to
! specify their actual dimensions, so that halos etc can be
! added or removed in the host model without having to change
! the convection code.
! However, this means it is up to the calling routine to ensure
! that nx_full,ny_full,k_bot_conv,k_top_conv are set consistent
! with the input arrays!

! Highest model-level where convection is allowed to initiate
integer :: k_top_init = 0
! Some of the fields used in the convective initiation
! calculation may not be available on all model-levels
! (e.g. turbulence fields only available up to the
!  highest boundary-layer level).
! You can set k_top_init lower than k_top_conv, to allow the
! scheme to do convection up to k_top_conv without referencing
! such fields out-of-bounds.
! Note: any fields on rho-levels used in the convective
! initiation calculations need to be available up to and
! including rho-level k_top_init + 1.

! Number of passive tracers to be carried by convection
integer :: n_tracers = 0

! Flag for whether each tracer is positive-only
! (most tracers probably shouldn't be allowed to go negative,
!  but some user-defined tracers might, so make this an option
!  that can be passed in)
logical, allocatable :: tracer_positive(:)


! BASIC MODEL SETTINGS

! Timestep used by the convection scheme
real(kind=real_cvprec) :: comorph_timestep = 600.0_real_cvprec

! Switch for approximating the dry adiabat using R_dry / cp_dry
! (needs to be set consistent with the host model's treatment)
logical :: l_approx_dry_adiabat = .true.

! Flag for using spherical coordinates
! (means the grid area increases with height,
!  for the purposes of conservation)
logical :: l_spherical_coord = .true.

! Flag for whether the convection scheme interacts with
! partial cloud fractions
! (convection increments the cloud-fractions if true).
logical :: l_cv_cloudfrac = .true.

! Flag for raining-out tracers within the convective drafts
logical :: l_tracer_scav = .false.

! Flags for whether to calculate CAPE / mass-flux-weighted CAPE
! (maybe needed either for diagnostics or for some as-yet unwritten
!  convective closure option)
logical :: l_calc_cape = .false.
logical :: l_calc_mfw_cape = .false.


! PHYSICS CONSTANTS; to BE OVERWRITTEN WITH HOST MODEL VALUES

! Acceleration due to gravity
real(kind=real_cvprec) :: gravity = 9.81_real_cvprec

! Melting point temperature of water / K
real(kind=real_cvprec) :: melt_temp = 273.0_real_cvprec

! Gas constants / J kg-1 K-1
real(kind=real_cvprec) :: R_dry = 287.058_real_cvprec ! Dry air
real(kind=real_cvprec) :: R_vap = 461.500_real_cvprec ! Water vapour

! Specific heat capacities / J kg-1 K-1
real(kind=real_cvprec) :: cp_dry = 1004.0_real_cvprec ! Dry air
real(kind=real_cvprec) :: cp_vap = 1860.0_real_cvprec ! Water vapour
real(kind=real_cvprec) :: cp_liq = 4181.0_real_cvprec ! Liquid water
real(kind=real_cvprec) :: cp_ice = 2110.0_real_cvprec ! Ice
! Note: If the heat capacities of vapour, ice and liquid water
! are set different to eachother, the temperature dependence
! of the latent heats of condensation, sublimation and fusion
! will naturally be accounted for.  If constant latent heats
! are desired, you must set cp_ice = cp_liq = cp_vap

! Latent heats at a given reference temperature / J kg-1
!                         Condensation / evaporation
real(kind=real_cvprec) :: L_con_ref = 2.501E6_real_cvprec
!                         Freezing / melting
real(kind=real_cvprec) :: L_fus_ref = 3.340E5_real_cvprec

! Reference temperature at which the above reference latent heats
! are valid (they optionally vary with temperature)
real(kind=real_cvprec) :: ref_temp_l = 273.0

! Densities of liquid water and ice / kg m-3
real(kind=real_cvprec) :: rho_liq = 1000.0_real_cvprec
real(kind=real_cvprec) :: rho_ice = 916.7_real_cvprec



!---------------------------------------------------------------
! CoMorph settings currently hardwired to defaults
!---------------------------------------------------------------


! GENERAL CONTROL OPTIONS...

! Number of updraft types for cloud spectrum model option
integer, parameter :: n_updraft_types = 1

! Number of independent downdraft types
! (distinct from "fall-back downdrafts", which are each tied
!  to a given updraft type).
integer :: n_dndraft_types = 1

! Max number of distinct convection layers that can be output
! separately for diagnostics
integer, parameter :: n_conv_layers_diag = 2

! Switches for whether to include fall-back flows for updrafts
! and downdrafts.
! The fall-back flow from an updraft consists of overshooting
! air detrained from the updraft which has sufficient negative
! buoyancy to fall down into the model-layers below its
! detrainment level.  A separate downdraft calculation
! (called the updraft fall-back) is done
! to return this detrained air to its neutral buoyancy level.
! Similarly, "undershooting" air detrained from downdrafts may
! need to rise back up to its neutral buoyancy level
! in a separate updraft calculation
! (called the downdraft fall-back)
logical, parameter :: l_updraft_fallback = .false.
logical, parameter :: l_dndraft_fallback = .false.
! Where there are multiple types of updraft / downdraft,
! each type is assumed to have its own separate fall-back flow.

! Max fraction of the layer-mass that updrafts, downdrafts,
! and their respective fall-back flows are allowed to entrain
! in one timestep (for numerical stability)
real(kind=real_cvprec), parameter :: max_ent_frac_up    = 0.4_real_cvprec
real(kind=real_cvprec), parameter :: max_ent_frac_up_fb = 0.0_real_cvprec
real(kind=real_cvprec), parameter :: max_ent_frac_dn    = 0.4_real_cvprec
real(kind=real_cvprec), parameter :: max_ent_frac_dn_fb = 0.0_real_cvprec
! Note that since updrafts and downdrafts (and their respective
! fall-back flows) may all act simultaneously, the max frac that
! might get entrained where all are active is the sum of these
! 4 numbers.
! Therefore, their sum must not exceed 1.0, or else a hell of
! nightmare numerics might be unleashed!
!
! Note that the fall-back mass-source doesn't need to be
! included in the quota, because it is simply a portion of the
! detrainment from the primary updraft/downdraft and so doesn't
! draw from the environment.
!
! Also note that these fractions impose limits on the fraction
! of the precip in a layer that can fall into the parcels, and
! on the momentum increments for numerical stability.
! Also, the convective closure rescaling must be restricted so
! that these limits are not exceeded
! (i.e. CFL limit on the closure)

! Max CFL number for entrainment and mass intitiation from layer
real(kind=real_cvprec), parameter :: max_cfl = 0.9_real_cvprec
! Note: this is used instead of the above max_ent_frac
! parameters if i_cfl_check is set to apply the CFL
! constraint afterwards via closure rescaling instead
! of applying the local limit on entrainment from each
! level during the ascent

! Flags for optional water species (whether they are represented
! inside the convection scheme and incremented by it)
logical, parameter :: l_cv_rain = .true.
logical, parameter :: l_cv_cf = .true.
logical, parameter :: l_cv_snow = .false.
logical, parameter :: l_cv_graup = .true.

! Switch for applying the "compensating subsidence" term
! inside the convection scheme
! (alternatively it could be done by the dynamical core
!  instead, as in the Kuell & Bott HYMACS scheme)
logical, parameter :: l_column_mass_rearrange = .true.

! Switch for adding initiating parcel perturbations based on
! the input turbulent fluxes and TKE
logical, parameter :: l_turb_par_gen = .true.

! Switch for vertically homogenising the convection increments
! below the boundary-layer top.
! This resets the resolved-scale source terms below z_bl_top,
! consistent with a linear increase of mass-flux with height,
! entraining air with a constant-with-height perturbation
! relative to the grid-mean profile.
logical, parameter :: l_homog_conv_bl = .true.

! Flag for converting precip that falls out of the parcel into
! cloud-water, in order to let the large-scale microphysics
! decide how much autoconversion to do instead.
logical, parameter :: l_precip_to_cloud = .false.

! Options for type of CFL-limit to apply:
! No CFL limit
integer, parameter :: i_cfl_check_none = 0
! Vertically-local limit applied to initiating and entrained mass
integer, parameter :: i_cfl_check_local = 1
! Rescale entire vertical profile of resolved-scale source terms
integer, parameter :: i_cfl_check_closure = 2
! Hybrid: use vertically-local CFL limit when above the BL-top,
! but apply the closure CFL limit after homogenising the
! convection within the BL
! (to be used in conjunction with l_homog_conv_bl).
integer, parameter :: i_cfl_check_hybrid = 3

! Set to one of the above options:
integer, parameter :: i_cfl_check = i_cfl_check_hybrid

! Options for diagnosed convective cloud:
! No convective cloud
integer, parameter :: i_convcloud_none = 0
! Liquid convective cloud only
integer, parameter :: i_convcloud_liqonly = 1
! Liquid, ice and mixed-phase convective cloud
integer, parameter :: i_convcloud_mph = 2

! Set to one of the above options
integer, parameter :: i_convcloud = i_convcloud_liqonly

! Switches for checking inputs and outputs for bad values
! (e.g. NaN, Inf, or negative values for positive-only fields)...

! Allowed values for the switches:
integer, parameter :: i_check_bad_none = 0  ! No bad-value checking
integer, parameter :: i_check_bad_warn = 1  ! Print warning if bad value found
integer, parameter :: i_check_bad_fatal = 2 ! Fatal error if bad value found

! Switch for checking the full 3-D fields input and output from comorph
integer, parameter :: i_check_bad_values_3d = i_check_bad_none
! Switch for checking compressed fields inside comorph
integer, parameter :: i_check_bad_values_cmpr = i_check_bad_none

! Switch for checking consistency of turbulence fields on input
integer, parameter :: i_check_turb_consistent = i_check_bad_none

! Switch for printing caps on turb perts
integer, parameter :: i_print_turb_perts_exceeded = i_check_bad_none

! Indirect indexing threshold
real(kind=real_cvprec) :: indi_thresh = 0.5_real_cvprec
! This is used in various places where a calculation only
! needs to be done at a sub-set of points.  If it is needed
! at the majority of points, it is most computationally
! efficient to just do it at all points.  But if it is only
! needed at a small minority of points it is faster to
! loop over a list of stored indices for the required points only.
! There are various calculations which use code of the form:
!
! ! If calculation needed at majorirty of points
! if ( real(nc) > indi_thresh * real(n_points) ) then
!
!   ! Do calculation at all points
!   ! (fully vectorised and so faster on many architectures)
!   do ic = 1, n_points
!     ...
!   end do
!
! ! If calculation only needed at a minoirty of points
! else
!
!   ! Do same calculation using indirect indexing list
!   ! (doesn't vectorise as efficiently).
!   do ic2 = 1, nc
!     ic = index_ic(ic2)
!     ...
!   end do
!
! end if
!
! This should allow better performance on average than just
! either always looping over all points or always using the
! indirect indexing list.  The value of indi_thresh should be tuned
! to optimise performance on a given platform.

! Compression threshold
real(kind=real_cvprec) :: cmpr_thresh = 0.5_real_cvprec
! Similar to indi_thresh above, but applies where instead of indirect
! indexing being used to deal with a calculation only needed
! at a subset of points, the required data are optionally
! copied into compression lists to make the work points
! contiguous in memory.  If the fraction of points where
! the calculation is needed exceeds cmpr_thresh, the gather /
! scatter is omitted, and instead the calculation is just done
! at all points.

! NOTE: Setting indi_thresh between 0 and 1 can cause the model
! to lose bit-reproducibility between different processor decompositions,
! because on some platforms the vectorised vs non-vectorised do-loops
! give slightly different answers, even for exactly the same
! floating point operations!


! MICROPHYSICS SETTINGS

! Properties of each condensed water species...

! Type definition to store the properties for each species:
type :: cond_params_type

  ! Abbreviated name of this species, used for labelling
  ! diagnostics and error messages relating to it.
  character(len=name_length) :: cond_name

  ! Flag for whether this species is ice or liquid
  logical :: l_ice

  ! Heat capacity of condensate (set to heat capacity of liquid
  ! water or ice depending on l_ice switch)
  real(kind=real_cvprec) :: cp = zero

  ! Density of hydrometeor / kg m-3
  ! (set equal to density of liquid water or ice, except for
  !  special case of lower-density rimed ice in graupel)
  real(kind=real_cvprec) :: rho = zero

  ! Number concentration of hydrometeor per unit dry-mass / kg-1
  real(kind=real_cvprec) :: n

  ! Flag for temperature-dependent number concentration.
  ! If true, n above is taken to be the number concentration
  ! at the melting point, with a parameterised increase
  ! at colder temperatures.
  logical :: l_tdep_n

  ! Effective area coefficient for hydrometeors; ratio of actual
  ! area to the area you get by assuming a sphere
  ! (so = 1 for cloud-droplets)
  real(kind=real_cvprec) :: area_coef

  ! Radius of a CCN particle on which the hydrometeors grow / m
  ! Sets the minimum allowed radius for moisture and heat
  ! exchange, and is used for freshly activated cloud droplets
  real(kind=real_cvprec) :: r_min

  ! Indicator for which sub-grid fraction this condensed water
  ! species exists within
  integer :: i_sg

  ! Indicator for the other species which this species converts
  ! into when it melts or freezes
  integer :: i_frzmlt

  ! Super-array address of the condensate species corresponding
  ! to i_frzmlt
  integer :: i_cond_frzmlt = 0

  ! Note that all the fields in this type which are set later
  ! depending on the values of other fields must be given default
  ! initialisation values, otherwise the use of the structure
  ! constructor to set the values for each condensed water
  ! species below gives a compile error due to these fields
  ! being omitted from the constructor.

end type cond_params_type

! Set properties for each species...
! (note: these structures would ideally be parameters,
!  but can't be because some of the components in the structures
!  depend on other settings and so have to be set at run-time,
!  in set_dependent_constants).

! Properties of liquid cloud:
type(cond_params_type), target :: params_cl = cond_params_type(                &
  cond_name = "cl",                                                            &
  l_ice     = .false.,                                                         &
  n         = 1.0E8_real_cvprec,     &  ! Liq cloud number conc.
  l_tdep_n  = .false.,               &  ! Use T-dependent number
  area_coef = 1.0_real_cvprec,       &  ! 1.0 for spheres
  r_min     = 1.0e-6_real_cvprec,    &  ! 1 micron CCN size
  i_sg      = i_sg_frac_liq,         &  ! Lives in the liquid cloud fraction
  i_frzmlt  = i_cf                  )   ! Freeze into ice cloud

! Properties of rain
type(cond_params_type), target :: params_rain = cond_params_type(              &
  cond_name = "rain",                                                          &
  l_ice     = .false.,                                                         &
  n         = 1000.0_real_cvprec,    &  ! Rain number conc. ~ 1 per l
  l_tdep_n  = .false.,               &  ! Use T-dependent number
  area_coef = 1.0_real_cvprec,       &  ! 1.0 for spheres
  r_min     = 0.0_real_cvprec,       &  ! No CCN for rain
  i_sg      = i_sg_frac_prec,        &  ! Lives in the precip fraction
  i_frzmlt  = i_graup               )   ! Freeze into graupel

! Properties of ice-cloud
type(cond_params_type), target :: params_cf = cond_params_type(                &
  cond_name = "cf",                                                            &
  l_ice     = .true.,                                                          &
  n         = 300.0_real_cvprec,     &  ! Ice cloud number at 0oC
  l_tdep_n  = .true.,                &  ! Use T-dependent number
  area_coef = 10.0_real_cvprec,      &  ! 10.0 for crystals
  r_min     = 5.0e-6_real_cvprec,    &  ! 5 micron CCN size
  i_sg      = i_sg_frac_ice,         &  ! Lives in the ice cloud fraction
  i_frzmlt  = i_rain                )   ! Melt into rain

! Properties of snow
type(cond_params_type), target :: params_snow = cond_params_type(              &
  cond_name = "snow",                                                          &
  l_ice     = .true.,                                                          &
  n         = 300.0_real_cvprec,     &  ! Snow number conc. ~0.3 per l
  l_tdep_n  = .false.,               &  ! T-depentent number conc. off.
  area_coef = 10.0_real_cvprec,      &  ! 10.0 for aggregates
  r_min     = 0.0_real_cvprec,       &  ! No CCN for snow
  i_sg      = i_sg_frac_ice,         &  ! Lives in the ice cloud fraction
  i_frzmlt  = i_rain                )   ! Melt into rain

! Properties of graupel
type(cond_params_type), target :: params_graup =cond_params_type(              &
  cond_name = "graup",                                                         &
  l_ice     = .true.,                                                          &
  n         = 100.0_real_cvprec,     &  ! Graupel number conc. ~0.1 per l
  l_tdep_n  = .false.,               &  ! T-depentent number conc. off.
  area_coef = 1.0_real_cvprec,       &  ! 1.0 for spheres
  r_min     = 0.0_real_cvprec,       &  ! No CCN for graupel
  i_sg      = i_sg_homog,            &  ! Assumed homogeneous across grid-box
  i_frzmlt  = i_rain                )   ! Melt into rain

! Density of rimed ice (used for graupel)
real(kind=real_cvprec) :: rho_rim = 600.0_real_cvprec

! Temperature-dependent ice number concentration slope
! The number concentration n(T) will be given by:
! n(T) = n0 exp( fac_tdep_n ( T - Tmelt ) )
! ( but limited above Tmelt and below T_homnuc)
real(kind=real_cvprec) :: fac_tdep_n = -one/8.18_real_cvprec
! set to 1/8.18 K-1, consistent with the UM microphysics.

! Ice nucleation
! Homogeneous freezing temperature / K
! All liquid is instantly frozen below this
real(kind=real_cvprec), parameter :: homnuc_temp = 233.0_real_cvprec
! Heterogeneous nucleation temeprature / K
! Gradual freezing starts below this
real(kind=real_cvprec) :: hetnuc_temp = 263.0_real_cvprec
! Heterogeneous freezing rate coefficient / s-1
! Freezing rate = coef_hetnuc * q_cl
real(kind=real_cvprec), parameter :: coef_hetnuc = 1.0e-7_real_cvprec
! Note: heterogeneous nucleation is only meant to act as a "seed"
! to initiate ice growth by other processes, so this coef should
! be set very small.  If the solution has much sensitivity
! to the value of this coef then something has gone wrong!

! Precipitation fall-speeds
! Asymptotic drag coefficient for a sphere at high Reynolds
! number limit
real(kind=real_cvprec) :: drag_coef_cond = 0.5_real_cvprec
! Kinematic viscocity of air / m2 s-1;
! controls drag at low Reynolds number
real(kind=real_cvprec), parameter :: kin_visc = 2.0e-5_real_cvprec

! Molecular diffusivities for exchange of heat and moisture
! with hydrometeors
! Diffusivity of water vapour in air / m2 s-1
real(kind=real_cvprec), parameter :: mol_diff_q = 2.8e-5_real_cvprec
! Thermal diffusivity of air / m2 s-1
real(kind=real_cvprec), parameter :: mol_diff_t = 1.9e-5_real_cvprec
! Coefficent scaling a term added onto the vapour and heat diffusion,
! for additional exchange due to fall-speed ventilation
real(kind=real_cvprec) :: vent_factor = 0.25_real_cvprec

! Accretion
! Coefficient for collision rate due to the spread of different
! fall-speeds (this is an approximation to avoid getting zero collision
! rate when the mean fall-speeds of 2 different hydrometeor
! species are equal; really we should integrate over size
! distributions to achieve this!)
real(kind=real_cvprec), parameter :: coef_wf_spread = 0.5_real_cvprec
! Coefficient for reduction of collection efficiency by
! deflection flow around hydrometeors
real(kind=real_cvprec) :: col_eff_coef = 1.0_real_cvprec

! Switch for autoconversion, of liquid cloud to rain
integer, parameter :: autoc_linear    = 1 ! coef_auto*( q_cl - q_cl_auto )
integer, parameter :: autoc_quadratic = 2 ! coef_auto*q_cl^2
integer :: autoc_opt = autoc_quadratic

! Autoconversion threshold, when using linear form / kg kg-1
real(kind=real_cvprec) :: q_cl_auto = 1.0e-3_real_cvprec

! Autoconversion rate coefficient / s-1
real(kind=real_cvprec) :: coef_auto = 1.0e-2_real_cvprec

! Environment temperature below which not to bother checking
! for melting, for computational efficiency / K
real(kind=real_cvprec), parameter :: melt_check_temp = 233.0_real_cvprec

! Miniscule condensed water mixing ratio set to activate
! new condensation
real(kind=real_cvprec), parameter :: q_activate = min_float

! Max number of iterations, tolerance for convergence,
! and error handling switch for convergence failure,
! in iterative solution of precip fall-speed
integer, parameter :: solve_wf_n_iter = 5
real(kind=real_cvprec), parameter :: solve_wf_tolerance = 0.05_real_cvprec
integer, parameter :: i_solve_wf_check_converge = i_check_bad_none

! Switch for doing run-time checks that conservation of
! heat, total water and momentum is being maintained
! (has same allowed settings as i_check_bad_values)
integer, parameter :: i_check_conservation = i_check_bad_none

! Flag for run-time checks for internal consistency of the
! implicit phase-change solution
! (has same allowed settings as i_check_bad_values)
integer, parameter :: i_check_imp_consistent = i_check_bad_none


! PLUME-MODEL SETTINGS

! Switch for including separate "core" properties in the parcel
logical, parameter :: l_par_core = .true.

! Switch to account for initiation mass-sources and mass-fluxes
! from other convection types/layers in the implicit detrainment
! calculation.
! (should improve the accuracy and stability, but makes the
! different convection types/layers interdependent and so breaks
! any assumption of their independence during the timestep).
logical, parameter :: l_inter_impl_det = .true.

! Implicitness weighting in the detrainment calculation
! (the environment buoyancy used has alpha * a pre-estimate of
!  the compensating subsidence increment added on, to improve
!  numerical stability and reduce intermittent / noisy behaviour)
real(kind=real_cvprec), parameter :: alpha_detrain = 1.0_real_cvprec

! Dimensionless constant scaling the initiation mass-sources
real(kind=real_cvprec) :: par_gen_mass_fac = 0.25_real_cvprec

! Prescribed draft vertical velocity excess, used when
! the vertical momentum equation is disabled / m s-1
real(kind=real_cvprec) :: wind_w_fac = 1.0_real_cvprec

! Switch to make the diagnosed convective cloud fraction depend
! on parcel buoyancy.  Can be used when not using the
! w-equation, but still want to make the convective fraction
! bigger (smaller) when slower (more rapid) ascent is expected.
logical, parameter :: l_cf_conv_buoy = .true.
! Tuning constant for buoyancy-dependent convective fraction:
! Assuming w' = fac * sqrt( buoyancy * radius )
real(kind=real_cvprec) :: wind_w_buoy_fac = 1.0_real_cvprec
! Minimum allowed implied w
real(kind=real_cvprec), parameter :: w_min = 0.1_real_cvprec

! Scaling factor for the convective cloud fraction, to account for
! in-plume w-variation
! (the slower-rising parts of the plume have much greater area)
real(kind=real_cvprec) :: cf_conv_fac = 2.0_real_cvprec

! Constants controlling the initial parcel perturbations:

! a) Those used if not using turbulence-based parcel properties
!    (l_turb_par_gen = .false.)

! Prescribed initial parcel radius
real(kind=real_cvprec), parameter :: par_gen_radius = 500.0_real_cvprec
! Fractional moisture perturbation to apply to the initial parcel
real(kind=real_cvprec), parameter :: par_gen_qpert = 0.05_real_cvprec

! b) Those used if using turbulence-based parcel properties
!    (l_turb_par_gen = .true.)

! Scaling factors for turbulence-based parcel perturbations
real(kind=real_cvprec), parameter :: par_gen_w_fac = 1.0_real_cvprec
real(kind=real_cvprec) :: par_gen_pert_fac = 0.667_real_cvprec
real(kind=real_cvprec) :: par_gen_radius_fac = 8.0_real_cvprec
! Background non-turbulent moisture perturbation
real(kind=real_cvprec) :: par_gen_rhpert = 0.05_real_cvprec

! Minimum parcel initial radius
! (assymptotic value above the BL-top; reduced near the surface)
real(kind=real_cvprec) :: ass_min_radius = 500.0_real_cvprec
! Factor for linear ramp of min radius near the surface
real(kind=real_cvprec), parameter ::  min_radius_fac = 0.25_real_cvprec

! Maximum allowed fractional water-vapour perturbation
real(kind=real_cvprec), parameter :: max_qpert = 1.0_real_cvprec


! Minimum allowed parcel radius, for safety
real(kind=real_cvprec), parameter :: min_radius = 5.0_real_cvprec

! Switch for how parcel radius evolves with height in the plume:
integer, parameter :: par_radius_evol_const = 0       ! Constant with height
integer, parameter :: par_radius_evol_volume = 1      ! (massflux/rho)^(1/3)
integer, parameter :: par_radius_evol_no_decrease = 2 ! as 1 but cannot decline
integer, parameter :: par_radius_evol_no_detrain = 3  ! as 1 but neglect detrain
integer :: par_radius_evol_method = par_radius_evol_no_decrease

! Scaling factor for par_gen core perturbations relative to
! the parcel mean properties (used if l_par_core = .true.)
real(kind=real_cvprec) :: par_gen_core_fac = 3.0_real_cvprec

! Minimum and maximum allowed values of the core/mean ratio of parcel
! buoyancies, which sets the power of the assumed power-law PDF of
! in-parcel buoyancy used for detrainment
real(kind=real_cvprec), parameter :: min_cmr = 2.0_real_cvprec
real(kind=real_cvprec), parameter :: max_cmr = 6.0_real_cvprec

! Maximum allowed convective area fraction
real(kind=real_cvprec), parameter :: max_sigma = 0.5_real_cvprec

! Method to use for setting the in-parcel cloud-fractions
! (if l_cv_cloudfrac is switched on)
! 1) Condensate assumed to be homogeneous across the parcel
!    (cloud fractions always either 0 or 1)
integer, parameter :: i_par_cloudfrac_hom = 1
! 2) Liquid and ice have partial fractions when both present
integer, parameter :: i_par_cloudfrac_mph = 2

integer, parameter :: i_par_cloudfrac = i_par_cloudfrac_mph
! This option makes convection detrain at least some cloud
! which is just liquid or just ice (not fully overlapped)
! when the parcel is mixed-phase.

! Power for overlap between liquid and ice cloud fractions
! inside the parcel
! overlap_power => 1 (no overlap)
! overlap_power => 0 (total overlap)
real(kind=real_cvprec) :: overlap_power = 0.5_real_cvprec
! Note: do not set to exactly zero as this causes a singularity!
! Instead, change i_par_cloudfrac = i_par_cloudfrac_hom
! to get total overlap between liquid and ice in the parcel.

! Max number of iterations, tolerance for convergence,
! and error handling switch for convergence failure,
! in iterative solution of detrained fraction
integer, parameter :: solve_det_n_iter = 10
real(kind=real_cvprec), parameter :: solve_det_tolerance = 0.001_real_cvprec
integer, parameter :: i_solve_det_check_converge = i_check_bad_none

! Pressure drag on thermals:
! Constant scaling the overall drag
real(kind=real_cvprec) :: drag_coef_par = 0.5_real_cvprec
! "Shape factor" scaling the stability-dependent term in the drag
real(kind=real_cvprec), parameter :: Sd = 3.0_real_cvprec

! Entrainment:
! Mixing entrainment rate (m-1) = ent_coef / parcel radius
real(kind=real_cvprec) :: ent_coef = 0.2_real_cvprec

! Factor for parcel core entrainment
real(kind=real_cvprec) :: core_ent_fac = 1.0_real_cvprec
! Switch for whether to include factor of 1/core_mean_ratio in core dilution
logical :: l_core_ent_cmr = .true.
! core_mean_ratio = core buoyancy / mean buoyancy of the parcel,
! and is used to construct the power-law PDF used for detrainment.
!
! Rate of entrainment of environment air into the parcel core is
! parameterised as...
! if l_core_ent_cmr:
!   core_ent = core_ent_fac * mean_ent / core_mean_ratio
! else:
!   core_ent = core_ent_fac * mean_ent


!---------------------------------------------------------------
! Things set at run-time, dependent on any of the above...
!---------------------------------------------------------------
! (these are set in set_dependent_constants)

! The following are set later depending on the flags
! l_cv_rain, l_cv_cf, etc:

! Number of condensed water species fields in use
integer :: n_cond_species = 0
! Number of liquid water species fields in use
integer :: n_cond_species_liq = 0
! Number of ice species fields in use
integer :: n_cond_species_ice = 0

! Addresses of condensed water species in condensate-only array
integer :: i_cond_cl = 0
integer :: i_cond_rain = 0
integer :: i_cond_cf = 0
integer :: i_cond_snow = 0
integer :: i_cond_graup = 0

! Latent heat of deposition / sublimation at the reference temperature
real(kind=real_cvprec) :: L_sub_ref  ! = L_con_ref + L_fus_ref

! Extrapolated latent heats at absolute zero
! (differ from L_con_ref etc if accounting for temperature-dependence).
real(kind=real_cvprec) :: L_con_0
real(kind=real_cvprec) :: L_fus_0
real(kind=real_cvprec) :: L_sub_0

! List of pointers to the properties of each condensed
! water species...
! Define a type which just contains a pointer to one of the
! condensate species properties structures above
type :: cond_params_pt_type
  type(cond_params_type), pointer :: pt => null()
end type cond_params_pt_type
! Declare a list of these; one for each species
type(cond_params_pt_type), allocatable :: cond_params(:)

! Newline character used in error messages etc
! (needs to be set using a call to the new_line intrinsic,
!  but this must be done in an executable run-time statement
!  so can't be hardwired in this module)
character(len=1) :: newline



end module comorph_constants_mod
