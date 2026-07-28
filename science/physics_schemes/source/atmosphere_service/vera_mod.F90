! *********************************COPYRIGHT************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *********************************COPYRIGHT************************************
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_vera
!
! Part of the Vera scheme for diagnosing visual range.
!
! This is the driver for the Vera visibility scheme.
!

module vera_mod

use vera_kind_mod, only: wp => vera_real, wi => vera_integer

use umPrintMgr,    only: umPrint

implicit none

! Description:
!   Implements the Vera scheme on set of inputs of (P, T, q, qcl, [am]),
!   and, if required, using synthetic input noise to estimate the error
!   in the computed visibility.
!
! Method:
!   This module comprises just the the single subroutine:
!
!   vera
!     Inputs are the meteorological fields (P, T, q, qcl, [am]):
!
!     P    pressure                       [Pa]
!     T    temperature                    [K]
!     q    specific humidity              [kg/kg]
!     qcl  cloud liquid water             [kg/kg]
!
!     An optional input is the aerosol mass mixing ratio, am:
!
!     am   MURK aerosol mass mixing ratio [micrograms/kg]
!
!     The output is the computed visibility [m].
!
!     The scattering coefficient computed from the precipitation can
!     also be used as an optional input - this allows Vera to include
!     the contribution of precipitation in the computation of visibility.
!
!     There are optional inputs to control the phantom aerosol
!     distribution and the treatment of synthetic noise to estimate
!     error in the computed visibility.
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_MOD'

private

! Switch controlling the one-off report of the Vera parameters to the log.
! Cleared by the first call to vera, under a critical region.
logical, save :: report_vera_parameters = .true.

! make the Vera driver a public routine
public :: vera

contains

subroutine vera( n_points, p_in, t_in, q_in, qcl_in        ,                   &
                 aerosol_mmr                               ,                   &
                 vera_config                               ,                   &
                 n_noise                                   ,                   &
                 covariance                                ,                   &
                 centiles_use                              ,                   &
                 ranges_use                                ,                   &
                 scattering_ls                             ,                   &
                 scattering_c                              ,                   &
                 fractional_ls                             ,                   &
                 fractional_c                              ,                   &
                 centiles                                  ,                   &
                 centiles_precip                           ,                   &
                 ranges                                    ,                   &
                 ranges_precip                             ,                   &
                 vis_median                                ,                   &
                 vis_median_precip                         ,                   &
                 vis_geometric_mean                        ,                   &
                 vis_geometric_mean_precip                   )

use vera_global_mod,  only: vera_constants          ,                          &
                            vera_aerosol            ,                          &
                            vera_phantom            ,                          &
                            vera_outputs            ,                          &
                            vera_noise_covariance   ,                          &
                            vera_noise_control      ,                          &
                            vera_ns                 ,                          &
                            vera_aerosol_population ,                          &
                            one_i => vera_one_i     ,                          &
                            one_r => vera_one

use vera_noise_mod ,  only: vera_generate_noise   ,                            &
                            vera_sort             ,                            &
                            vera_centile          ,                            &
                            vera_range

use vera_scheme_mod,  only: vera_scheme    ,                                   &
                            vera_murk_cast

use vera_phantom_list_mod, only: vera_phantom_list

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera
!
! NOTE: The atmospheric inputs (P, T, q, qcl, [am]) are all vectors
!       i.e. 1-d arrays and must all have the same number of elements.
!
!       The output visibility is also a 1-d array, matching the number
!       of elements matching the number of input sets.
!
!       The optional inputs to vera:
!
!       aerosol_mmr - aerosol loading in micrograms/kg. If this input
!                     is not a 1-d array with size matching the input P.
!                     then the background aerosol loading is used,
!                     i.e. vera_aerosol%am_background in the Vera global
!                     data module.
!
!       vera_config - the configuration to use, as defined
!                     in the module vera_phantom_list_mod. The default is to
!                     use the value given by vera_phantom%vera_config in
!                     the Vera global data module.
!
!       n_noise     - how many synthetic noise pertubations generate for
!                     each metorological input. If this is set to zero, then
!                     no synthetic noise is generated and the input sets
!                     of (P, T, q, qcl, [am]) are processed unperturbed
!                     by Vera.
!
!                     The default is to use the value given by
!                     vera_noise_control%n_noise in the Vera global data
!                     module, which is set to zero by default.
!
!       covariance  - if correlated synthetic noise is required, then this
!                     input provides a square covariance matrix. A possible
!                     covariance matrix is given by
!                     vera_noise_control%noise_covariance in the Vera global
!                     data module. If this optional input is ommitted, then
!                     uncorrelated noise is generated.
!
!       centiles_use - if synthetic noise is used, then this argument
!                     can be used to specify the centiles to use when
!                     computing centiles from a noisy visibility
!                     population. This input is a 1-d real array. If this
!                     argument is ommitted then the default values in the
!                     Vera global module will be used.
!
!       ranges_use   - if synthetic noise is used, then this argument can be
!                     used to specify the visible ranges to use when
!                     computing probabilities from a noisy visibility
!                     population. This input is a 1-d real array. If this
!                     argument is ommitted then the default values in the
!                     Vera global module will be used.
!
!       scattering_ls - scattering coefficients computed from the large-
!                     scale precipitation. This is used to compute the
!                     total visibility in precipitation using Koschmeider's
!                     Law.
!
!       scattering_c  - scattering coefficients computed from the convective
!                     scale precipitation. This is used to compute the
!                     total visibility in precipitation using Koschmeider's
!                     Law.
!
!       fractional_ls - fractional coverage of large-scale cloud.
!
!       fractional_c  - fractional coverage of convective cloud.
!
!       The optional outputs from vera:
!
!       vis_median   - median visibility computed from a set of noisy
!                      visibilities, using just the contributions from the
!                      aerosol and Rayleigh scattering.
!
!       vis_median_precip - median visibility computed from a set of noisy
!                      visibilities, using the contributions from the
!                      aerosol and Rayleigh scattering together with
!                      scattering due to precipitation. Computing this
!                      requires the optional inputs scattering_ls and
!                      scattering_c.
!
!       vis_geometric_mean - geometric mean visibility computed from a set
!                      of noisy visibilities, using just the contributions
!                       from the aerosol and Rayleigh scattering.
!
!       vis_geometric_mean_precip - geometric mean visibility computed from
!                      a set of noisy visibilities, using the contributions
!                      from the aerosol and Rayleigh scattering together
!                      with scattering due to precipitation. Computing this
!                      requires the optional inputs scattering_ls and
!                      scattering_c.
!
!       centiles - if synthetic noise is used, then this argument can be
!                     used to specify the output array for the centiles
!                     computed from the visibility population. The size of
!                     this ouptut array is given by
!                     (n_centiles, n_inputs) where n_centiles is the
!                     number of centiles to compute, and n_inputs is the
!                     number of input sets of (P, T, q, qcl, [am]). The
!                     centiles to compute are specified by the 1-d array
!                     vera_noise_control%centiles in the Vera global data
!                     module.
!
!       centiles_precip - similar to the optional output centiles, but
!                     including the scattering contribution from
!                     precipitation. Computing this output requires the
!                     optional inputs scattering_ls and scattering_c.
!
!       ranges      - if synthetic noise is used, then this argument can be
!                     used to specify the output array for the probabilities
!                     of visual ranges computed from the noisy visibility
!                     population. The size of this ouptut array is given by
!                     (n_ranges, n_inputs) where n_ranges is the
!                     number of ranges to compute, and n_inputs is the
!                     number of input sets of (P, T, q, qcl, [am]). The
!                     ranges to compute are specified by the 1-d array
!                     vera_noise_control%ranges in the Vera global data
!                     module.
!
!       ranges_precip - similar to the optional output ranges, but
!                     including the scattering contribution from
!                     precipitation. Computing this output requires the
!                     optional inputs scattering_ls and scattering_c.
!
!===========================================================================

! number of meteorlogical inputs sets of (P, T, q, qcl, [am])
integer (wi),           intent (in)    :: n_points

! atmospheric pressure [Pa]
real    (wp),           intent (in)    :: p_in(n_points)

! atmospheric temperature [K]
real    (wp),           intent (in)    :: t_in(n_points)

! atmospheric specific humidity [kg/kg]
real    (wp),           intent (in)    :: q_in(n_points)

! atmospheric liquid water [kg/kg]
real    (wp),           intent (in)    :: qcl_in(n_points)

! configuration to use
integer (wi), optional, intent (in)    :: vera_config

! how many instances of synthetic noise to produce?
! If n_noise = 0, then no synthetic noise is produced. The default is
! defined by vera_noise_control%n_noise in vera_global_mod.F90
integer (wi), optional, intent (in)    :: n_noise

! aerosol mass mixing ratio [micrograms/kg]
real    (wp), optional, intent (in)    :: aerosol_mmr(n_points)

! covariance matrix for the synthetic noise
real    (wp), optional, intent (in)    :: covariance(:, :)

! large scale precipitation scattering coefficients used to
!  compute total visibility in precipitation
real    (wp), optional, intent (in)    :: scattering_ls(n_points)

! convective precipitation scattering coefficients used to
!  compute total visibility in precipitation
real    (wp), optional, intent (in)    :: scattering_c(n_points)

! total large scale fractional cloud coverage used to compute total
! visibility in precipitation
real    (wp), optional, intent (in)    :: fractional_ls(n_points)

! convective fractional cloud coverage used to compute total
! visibility in precipitation
real    (wp), optional, intent (in)    :: fractional_c(n_points)

! centiles to compute from the pdf of noisy visibilities
real    (wp), optional, intent (in) :: centiles_use(:)

! visible ranges to estimate probabilities from using the
! pdf of noisy visibilities
real    (wp), optional, intent (in) :: ranges_use(:)

! centiles computed from set of noisy visibilities [m]
real    (wp), optional, intent (in out):: centiles(:, :)

! centiles computed from set of noisy visibilities [m],
! including a contribution from precipitation
real    (wp), optional, intent (in out):: centiles_precip(:, :)

! probabilities of visibility being less than threshold ranges
real    (wp), optional, intent (in out):: ranges(:, :)

! probabilities of visibility being less than threshold ranges,
! including a contribution from precipitation
real    (wp), optional, intent (in out):: ranges_precip(:, :)

! median visibility [m] computed from a set of noisy visibilities
real    (wp), optional, intent (in out):: vis_median(n_points)

! median visibility [m] including a contribution from precipitation
real    (wp), optional, intent (in out):: vis_median_precip(n_points)

! geometric_mean visibility [m] computed from a set of noisy visibilities
real    (wp), optional, intent (in out):: vis_geometric_mean(n_points)

! geometric_mean visibility [m] including a contribution from precipitation
real    (wp), optional,intent (in out):: vis_geometric_mean_precip(n_points)

!===========================================================================
! local variables for vera
!===========================================================================

! atmospheric pressure [Pa]
real    (wp)                    :: p(n_points)

! atmospheric temperature [K]
real    (wp)                    :: t(n_points)

! atmospheric specific humidity [kg/kg]
real    (wp)                    :: q(n_points)

! atmospheric liquid water [kg/kg]
real    (wp)                    :: qcl(n_points)

! aerosol mass mixing ratio [micrograms/kg]
real    (wp)                    :: am(n_points)

! computed visibility [m], using just aerosol and Rayleigh scattering,
! i.e. no contribution from precipitation
real    (wp)                    :: vis_no_precip(n_points)

! computed visibility [m], including the contribution from precipitation
real    (wp)                    :: vis_precip(n_points)

! large scale precipitation scattering coefficients used to
!  compute total visibility in precipitation
real    (wp)                    :: precip_scattering_ls(n_points)

! convective precipitation scattering coefficients used to
!  compute total visibility in precipitation
real    (wp)                    :: precip_scattering_c(n_points)

! fractional large scale cloud cover used to compute total
! visibility in precipitation
real    (wp)                    :: precip_fractional_ls(n_points)

! fractional convective cloud cover used to compute total
! visibility in precipitation
real    (wp)                    :: precip_fractional_c(n_points)

! configuration to use - could be the default value in the global data
! module, or could use a value input to this routine
integer (wi)                    :: vera_config_use

! number of synthetic noises to produce - i.e. the number passed
! to this routine, or default number
integer (wi)                    :: n_noise_use

! counter to loop over the input sets of (P, T, q, qcl, [am])
integer (wi)                    :: ii_point

! counter to loop over the centiles required
integer (wi)                    :: ii_centile

! size of the random_number seed on entering Vera
integer (wi)                    :: seed_size

! store for random number generator random_number seed on entering Vera
integer (wi), allocatable       :: seed_store(:)

! synthetic noise to perturb the input (P, T, q, qcl, [am])
real    (wp), allocatable       :: noise(:, :)

! computed noisy visibilities, using just aerosol and Rayleigh scattering,
! i.e. no contribution from precipitation
real    (wp), allocatable       :: noisy_vis_no_precip(:)

! computed noisy visibilities, including the contribution from precipitation
real    (wp), allocatable       :: noisy_vis_precip(:)

! blank array of zero noise
real    (wp), allocatable       :: no_noise(:)

! a single element array to store the median visibility computed from a
! noisy ensemble
real    (wp)                    :: median_visibility(1)

! buffer for the call to umPrint
character(len=50000) :: lineBuffer

! Switch to determine whether we initialise parameters in vera_phantom_list().
! We should only need to do this on the first call to vera_phantom_list().
logical, parameter :: initialise_vera_parameters_on = .true.

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start of the executable code for vera
!===========================================================================

! grab the meteorological inputs
p   = p_in
t   = t_in
q   = q_in
qcl = qcl_in

!---------------------------------------------------------------------------
! gather the optional inputs to this subroutine

! if a configuration was passed to this routine, then use that
if (present(vera_config)) then
  vera_config_use            = vera_config
  vera_phantom%vera_config   = vera_config_use
else
  ! use the default configuration
  vera_config_use = vera_phantom%vera_config
end if

! check to see if aerosol mass mixing ratios have been input
if ( present(aerosol_mmr) .and.                                                &
     (vera_aerosol%aerosol_source == vera_aerosol%aerosol_source_murk) ) then
  ! use the input aerosol mass mixing ratios
  am = aerosol_mmr
else
  ! no aerosol mass mixing ratios have been passed to Vera, or the
  ! aerosol_source switch is not set to use MURK, so set the aerosol
  ! mass mixing ratio field to the background level
  am = vera_aerosol%am_background
end if

! check to see if large scale precipitation scattering coefficients
! have been input
if (present(scattering_ls)) then
  ! use the input precipitation scattering coefficient
  precip_scattering_ls = scattering_ls
else
  ! set the precipitation scattering coefficient to zero
  precip_scattering_ls = vera_constants%zero
end if

! check to see if convective precipitation scattering coefficients
! have been input
if (present(scattering_c)) then
  ! use the input precipitation scattering coefficient
  precip_scattering_c = scattering_c
else
  ! set the precipitation scattering coefficient to zero
  precip_scattering_c = vera_constants%zero
end if

! check for fractional large scale cloud coverage
if (present(fractional_ls)) then
  ! use the input fractional large scale cloud coverage
  precip_fractional_ls = fractional_ls
else
  ! set the large scale cloud coverage to zero
  precip_fractional_ls = vera_constants%zero
end if

! check for fractional convective cloud coverage
if (present(fractional_c)) then
  ! use the input fractional convective cloud coverage
  precip_fractional_c = fractional_c
else
  ! set the convective cloud coverage to zero
  precip_fractional_c = vera_constants%zero
end if

! if n_noise was passed to this routine, then use that
if (present(n_noise)) then
  n_noise_use                = n_noise
  vera_noise_control%n_noise = n_noise_use
else
  ! use the default number of synthetic noises
  n_noise_use = vera_noise_control%n_noise
end if

! check if any centiles were passed to this procedure
if (present(centiles_use)) then

  ! use the centiles passed to this procedure
  vera_noise_control%n_centiles = size( centiles_use )
  vera_noise_control%centiles(vera_constants%one_i         :                   &
                              vera_noise_control%n_centiles) = centiles_use

else

  ! use the default values specified in the Vera global data module
  vera_noise_control%n_centiles =                                              &
                       size(vera_noise_control % centiles_default)
  vera_noise_control%centiles(vera_constants%one_i         :                   &
                              vera_noise_control%n_centiles) =                 &
  vera_noise_control%centiles_default

end if

! check if any ranges were passed to this procedure
if (present(ranges_use)) then

  ! use the ranges passed to this procedure
  vera_noise_control%n_ranges = size( ranges_use )
  vera_noise_control%ranges(vera_constants%one_i       :                       &
                            vera_noise_control%n_ranges) = ranges_use

else

  ! use the default values specified in the Vera global data module
  vera_noise_control%n_ranges = size(vera_noise_control%ranges_default)
  vera_noise_control%ranges(vera_constants%one_i       :                       &
                            vera_noise_control%n_ranges) =                     &
  vera_noise_control%ranges_default

end if

! end of gathering optional inputs to this routine
!---------------------------------------------------------------------------

! Write some output to the log file.  The UM calls this routine once per
! timestep, but here it is called once per column from inside an OpenMP
! parallel do, so the parameters are reported on the first call only.
!$OMP CRITICAL(vera_report_parameters)
if ( report_vera_parameters ) then

  write(lineBuffer,'(A)')' ------------- Vera parameters -------------'
  call umPrint(lineBuffer,src='vera_mod')

  write(lineBuffer,'(A,I10)')' vera_config = ',    vera_phantom%vera_config
  call umPrint(lineBuffer,src='vera_mod')

  write(lineBuffer,'(A,I10)')' aerosol_source = ', vera_aerosol%aerosol_source
  call umPrint(lineBuffer,src='vera_mod')

  write(lineBuffer,'(A,I10)')' n_rd = ',           vera_phantom%n_rd
  call umPrint(lineBuffer,src='vera_mod')

  write(lineBuffer,'(A,I10)')' n_b0 = ',           vera_phantom%n_b0
  call umPrint(lineBuffer,src='vera_mod')

  write(lineBuffer,'(A,I10)')' n_noise_use = ',    n_noise_use
  call umPrint(lineBuffer,src='vera_mod')

  write(lineBuffer,'(A)')' ------ end of Vera parameters -------------'
  call umPrint(lineBuffer,src='vera_mod')

  report_vera_parameters = .false.

end if
!$OMP end CRITICAL(vera_report_parameters)

! push the meteorological inputs (P, T, q, qcl, [am]) through the Vera scheme,
! checking to see if any synthetic noise is required on the inputs
if ( n_noise_use == vera_noise_control%noise_switch_off ) then

  ! no synthetic noise required, so just push the input (P, T, q, qcl, [am])
  ! through the Vera scheme, remembering to cast a new phantom aerosol
  ! population for each value of am
  vera_phantom % cast_switch = vera_phantom % cast_switch_on
  call vera_scheme( p, t, q, qcl, am                         ,                 &
                    vera_config       = vera_config_use      ,                 &
                    scattering_ls     = precip_scattering_ls ,                 &
                    scattering_c      = precip_scattering_c  ,                 &
                    fractional_ls     = precip_fractional_ls ,                 &
                    fractional_c      = precip_fractional_c  ,                 &
                    vis_no_precip     = vis_no_precip        ,                 &
                    vis_precip        = vis_precip             )

  ! compute the required optional outputs, starting with ...
  ! ... median visibility from just aerosol and Rayleigh scattering,
  ! i.e. not including precipitation
  if (present(vis_median)) then
    vis_median = vis_no_precip
  end if

  ! median visibility including contribution from precipitation
  if (present(vis_median_precip)) then
    vis_median_precip = vis_precip
  end if

  ! geometric_mean visibility from just aerosol and Rayleigh scattering,
  ! i.e. not including precipitation
  if (present(vis_geometric_mean)) then
    vis_geometric_mean = vis_no_precip
  end if

  ! geometric_mean visibility including contribution from precipitation
  if (present(vis_geometric_mean_precip)) then
    vis_geometric_mean_precip = vis_precip
  end if

  ! if required, put the median visibility into the centiles ouput,
  ! computed using just aerosol and Rayleigh scattering, i.e.
  ! without precipitation ...
  if (present(centiles)) then
    do ii_centile = vera_constants%one_i, vera_noise_control%n_centiles
      centiles(ii_centile, :) = vis_no_precip
    end do
  end if

  ! ... or with precipitation
  if (present(centiles_precip)) then
    do ii_centile = vera_constants%one_i, vera_noise_control%n_centiles
      centiles_precip(ii_centile, :) = vis_precip
    end do
  end if

else

  ! need to generate some synthetic noise for the input sets of
  ! (P, T, q, qcl, [am])
  allocate( noise(vera_noise_covariance%n_variables, n_noise_use) )

  ! allocate the arrays to store the computed noisy visibilities
  allocate( noisy_vis_no_precip(n_noise_use) )
  allocate( noisy_vis_precip(n_noise_use) )
  allocate( no_noise(n_noise_use) )

  ! make sure that the no_noise array really is set to zero
  no_noise = vera_constants%zero

  ! Store the value of the current random number generator seed, so that
  ! afterwards this can be put back so as not to alter the prognostic
  ! evolution of the model when random numbers are changed by Vera but used
  ! elsewhere in stochastic physics.
  !
  ! The generator state and the Vera seed switch are both shared, so the
  ! seeding and the noise generation are done in a critical region: the UM
  ! calls this routine once per timestep from a serial region, whereas here
  ! it is called once per column from inside an OpenMP parallel do.  Every
  ! column seeds the generator in the same way and so sees the same noise
  ! realisation, which keeps the diagnostic independent of the domain
  ! decomposition.
!$OMP CRITICAL(vera_random_seed)
  call random_seed( size = seed_size )
  allocate ( seed_store(seed_size) )
  call random_seed( get  = seed_store )

  if ( (present(covariance)) .or.                                              &
       (vera_noise_covariance%switch_covariance ==                             &
        vera_noise_covariance%switch_covariance_on) ) then

    ! generate correlated noise
    if (present(covariance)) then
      ! use the input covariance matrix
      call vera_generate_noise( noise, Covariance = covariance )
    else
      ! use the default covariance matrix
      call vera_generate_noise( noise,                                         &
                             Covariance = vera_noise_covariance%cov_matrix )
    end if

  else

    ! generate uncorrelated noise, by calling the noise generation routine
    ! without passing a noise covariance matrix
    call vera_generate_noise( noise )

  end if

  ! set the switch to seed the random number generator
  ! so that the next call to Vera will seed the generator
  vera_noise_control % seed_flag = vera_noise_control % seed_flag_reset

  ! put back the stored value of the random number generator seed
  call random_seed( Put = seed_store )
  deallocate ( seed_store )
!$OMP end CRITICAL(vera_random_seed)

  ! scale the synthetic input noise for P and T, i.e. the two meteorological
  ! inputs that have noise that scales as fixed width Gaussian and fixed
  ! fixed bias, i.e. delta_x = N(gaussian) * st_x + bias_x
  noise(1,:) = ( noise(1,:)*vera_ns%width_p  ) + vera_ns%bias_p
  noise(2,:) = ( noise(2,:)*vera_ns%width_t  ) + vera_ns%bias_t

  ! check to see if the aerosol input is noisy
  if ( vera_noise_control % aerosol_noise ==                                   &
       vera_noise_control % aerosol_noise_off ) then

    ! no noise on the aerosol input, only need to cast the phantom
    ! aerosol population just the once
    call vera_phantom_list(vera_config_use, initialise_vera_parameters_on)
    vera_phantom % cast_switch = vera_phantom % cast_switch_off

  else

    ! noisy aerosol, so cast the phantom aerosol population for
    ! every noisy input
    vera_phantom % cast_switch = vera_phantom % cast_switch_on

  end if

  ! loop over the sets of input (P, T, q, qcl, [am])
  !
  ! Note: the OpenMP parallel region that wraps this loop in the UM is
  ! omitted here.  This routine is called from a kernel that is itself
  ! executed inside an OpenMP parallel do over columns, so the region would
  ! be a nested one.  With n_points = 1 it can win nothing, and were nesting
  ! ever enabled the inner threads would see their own, unallocated copies
  ! of the threadprivate aerosol population set up by the caller.
  do ii_point = vera_constants%one_i, n_points

    ! check to see if the aerosol input is noisy
    if ( vera_noise_control % aerosol_noise ==                                 &
         vera_noise_control % aerosol_noise_off ) then

      ! no noise on the aerosol input, so need to cast the phantom
      ! aerosol population once per gridpoint
      call vera_murk_cast( am(ii_point), vera_config_use )

    end if

    ! push the noisy (P, T, q, qcl, am) through the Vera scheme,
    ! scaling the inputs q, qcl, am using relative width Gaussian noise,
    ! i.e. delta_x = ( N(gaussian) * width_x * x ) + bias_x
    call vera_scheme( p(  ii_point) + noise(1, :)                       ,      &
                      t(  ii_point) + noise(2, :)                       ,      &
      (   q(  ii_point) * ( one_r + (noise(3, :)*vera_ns%width_q  ) ) ) +      &
                                                         vera_ns%bias_q ,      &
      ( qcl(  ii_point) * ( one_r + (noise(4, :)*vera_ns%width_qcl) ) ) +      &
                                                       vera_ns%bias_qcl ,      &
      (  am(  ii_point) * ( one_r + (noise(5, :)*vera_ns%width_am ) ) ) +      &
                                                        vera_ns%bias_am ,      &
                      vera_config   = vera_config_use                   ,      &
                      scattering_ls =                                          &
                               precip_scattering_ls(ii_point) + no_noise,      &
                      scattering_c  =                                          &
                               precip_scattering_c(ii_point)  + no_noise,      &
                      fractional_ls = precip_fractional_ls(ii_point) +         &
                                      no_noise                          ,      &
                      fractional_c  = precip_fractional_c(ii_point)  +         &
                                      no_noise                          ,      &
                      vis_no_precip = noisy_vis_no_precip               ,      &
                      vis_precip    = noisy_vis_precip                    )

    ! sort the list of noisy visibilities
    call vera_sort( noisy_vis_no_precip )
    call vera_sort( noisy_vis_precip )

    ! compute the required optional outputs, starting with ...
    ! ... median visibility from just aerosol and Rayleigh scattering,
    ! i.e. not including precipitation
    if (present(vis_median)) then
      ! compute the median of the noisy visibilities
      call vera_centile( noisy_vis_no_precip, median_visibility  ,             &
                         centile = vera_outputs % median_centile ,             &
                         sort    =  vera_noise_control % sort_off  )
      vis_median(ii_point) = median_visibility(vera_constants%one_i)
    end if

    ! median visibility including contribution from precipitation
    if (present(vis_median_precip)) then
      ! compute the median of the noisy visibilities, including the
      ! contribution of precipitation
      call vera_centile( noisy_vis_precip, median_visibility     ,             &
                         centile = vera_outputs % median_centile ,             &
                         sort    =  vera_noise_control % sort_off  )
      vis_median_precip(ii_point) = median_visibility(vera_constants%one_i)
    end if

    ! geometric mean visibility from just aerosol and Rayleigh scattering,
    ! i.e. not including precipitation
    if (present(vis_geometric_mean)) then
      vis_geometric_mean(ii_point) =                                           &
                                   exp( sum( log(noisy_vis_no_precip) ) /      &
                                        size(noisy_vis_no_precip,              &
                                        dim  = one_i            ,              &
                                        kind = wp                 )       )
    end if

    ! geometric mean visibility including contribution from precipitation
    if (present(vis_geometric_mean_precip)) then
      vis_geometric_mean_precip(ii_point) =                                    &
                                   exp( sum( log(noisy_vis_precip) ) /         &
                                        size(noisy_vis_precip,                 &
                                        dim  = one_i         ,                 &
                                        kind = wp              )       )
    end if

    ! if required, compute centiles from the noisy visibilities
    ! computed using just aerosol and Rayleigh scattering, i.e.
    ! no contribution from precipitation
    if (present(centiles)) then
      call vera_centile( noisy_vis_no_precip                   ,               &
                         centiles(:, ii_point)                 ,               &
                         sort =  vera_noise_control % sort_off   )
    end if

    ! if required, compute centiles from the noisy visibilities
    ! computed by including the contribution from precipitation
    if (present(centiles_precip)) then
      call vera_centile( noisy_vis_precip                      ,               &
                         centiles_precip(:, ii_point)          ,               &
                         sort =  vera_noise_control % sort_off   )
    end if

    ! if required, compute probabilities of visible ranges from
    ! the noisy visibilities computed using just aerosol and
    ! Rayleigh scattering, i.e. no contribution from precipitation
    if (present(ranges)) then
      call vera_range( noisy_vis_no_precip                   ,                 &
                       ranges(:, ii_point)                   ,                 &
                       sort =  vera_noise_control % sort_off   )
    end if

    ! if required, compute probabilities of visible ranges from
    ! the noisy visibilities computed by including the contribution
    ! from precipitation
    if (present(ranges_precip)) then
      call vera_range( noisy_vis_precip                      ,                 &
                       ranges_precip(:, ii_point)            ,                 &
                       sort =  vera_noise_control % sort_off   )
    end if

    ! end of looping over the input sets of (P, T, q, qcl, [am])
  end do

  ! tidy up
  if ( allocated( no_noise            ) ) deallocate( no_noise            )
  if ( allocated( noisy_vis_precip    ) ) deallocate( noisy_vis_precip    )
  if ( allocated( noisy_vis_no_precip ) ) deallocate( noisy_vis_no_precip )
  if ( allocated( noise               ) ) deallocate( noise               )
  ! tidy up the phantom aerosol population
  if ( allocated( vera_aerosol_population ) )                                  &
    deallocate( vera_aerosol_population )

  ! end of checking if any synthetic noise is required on the inputs
end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera

end module vera_mod
