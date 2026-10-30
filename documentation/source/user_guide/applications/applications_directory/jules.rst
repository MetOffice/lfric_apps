.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _jules:

The standalone land surface model - jules
=========================================

The `jules` application runs the JULES land surface model inside LFRic,
without the rest of the atmosphere model. It is meant to reproduce a
standalone JULES point (site) run. Observed near-surface meteorology drives
JULES, and the results can be checked directly against standalone JULES. It
also gives a cheap way to test the LFRic JULES interface.

The application uses the Gungho driver with ``method='jules'`` in the
``formulation`` namelist. Each timestep runs only the JULES surface
calculation (``jules_timestep_alg``). The dynamical core, transport, boundary
layer and radiation schemes are not run. The model is a single column (a
1x1 planar mesh) with one atmospheric layer. That layer represents the
forcing level.

The site forcing data is read through XIOS, and the code that applies it is
in ``interfaces/jules_interface``. Other applications could use it too.

Site forcing
------------

Switching it on
~~~~~~~~~~~~~~~

Two namelists switch the forcing on:

.. code-block:: fortran

   &files
     jules_forcing_path = '/path/to/forcing_file'   ! no .nc suffix
   /
   &initialization
     jules_forcing       = .true.
     jules_forcing_z1_tq = 10.0   ! height of the T and q data (m)
     jules_forcing_z1_uv = 10.0   ! height of the wind data (m)
   /

The forcing file
~~~~~~~~~~~~~~~~

The forcing file is a NetCDF file with one point, using the dimensions
``(time, land)``. It must contain the following variables:

=========  ==================================  ==============
Variable   Description                         Units
=========  ==================================  ==============
``Wind``   Wind speed                          m s\ :sup:`-1`
``Tair``   Air temperature                     K
``Qair``   Specific humidity                   kg kg\ :sup:`-1`
``PSurf``  Surface pressure                    Pa
``LWdown`` Downward longwave flux              W m\ :sup:`-2`
``SWdown`` Downward shortwave flux             W m\ :sup:`-2`
``Rainf``  Total precipitation rate            kg m\ :sup:`-2` s\ :sup:`-1`
=========  ==================================  ==============

It must also have a CF ``time`` coordinate variable with units of the form
``seconds since <date>``. The ``time`` dimension must be fixed, not
unlimited, because its length sets the size of the XIOS time axis. That
length is read when the model starts
(``get_jules_forcing_dim`` in ``time_dimensions_mod``). The file does not
need any UGRID mesh information.

The XIOS file definition is in ``iodef.xml``. It reads the variables onto an
axis-only grid (``jules_forcing_grid``), which avoids needing a mesh:

.. code-block:: xml

   <file id="jules_forcing" name="jules_forcing" mode="read"
         output_freq="30mi" convention="CF" cyclic="true" enabled=".FALSE.">
     <field id="jules_forcing_time" name="time" unit="seconds"
            operation="once" axis_ref="jules_forcing_time_axis" />
     <field id="forcing_wind_data" name="Wind" operation="instant"
            grid_ref="jules_forcing_grid" />
     ...
   </file>

The file is registered with XIOS in ``gungho_setup_io_mod`` when
``jules_forcing`` is true. Each XIOS field id is the model field name with
``_data`` added, for example ``forcing_wind`` and ``forcing_wind_data``. The
model fields are held in the ``jules_forcing_fields`` collection, which
``create_jules_forcing_fields`` in ``init_jules_forcing_mod`` creates.

.. important::

   The ``output_freq`` of the ``jules_forcing`` file must equal the data
   interval of the forcing file. It must also equal the model timestep.
   The data is not interpolated in time.

Timing
~~~~~~

The forcing time axis is set up without interpolation and with a window size
of one. Each data time is the start of the period that the data applies to,
and the value is held until the next data time. This is the ``'nf'``
interpolation option in standalone JULES. The forcing fields are updated in
``gungho_driver_mod`` after the clock ticks and before ``gungho_step``. As a
result, the first record in the file is used for the first timestep.

How the forcing drives JULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each timestep, before the surface calculations, ``apply_jules_forcing_alg``
calls ``jules_forcing_kernel`` to set the model state from the forcing:

* The surface Exner pressure is calculated from ``PSurf``. The Exner pressure
  at the forcing level comes from hydrostatic balance at height
  ``jules_forcing_z1_tq``.
* Potential temperature in the first layer is ``Tair`` divided by the Exner
  pressure there. The vapour mixing ratio is set to ``Qair``.
* The wind ``u`` is set to ``Wind`` and ``v`` is set to zero. Density comes
  from the equation of state.
* The downward SW and LW surface fluxes come from ``SWdown`` and ``LWdown``.
  The photosynthetically active radiation (PAR) is ``visible_fraction``
  (0.5) of the SW. The direct part of PAR is ``1 - diffuse_fraction`` (0.6),
  which matches ``diff_frac = 0.4`` in the JULES kernels.
* ``Rainf`` falls as rain when ``Tair`` is at least 274 K and as snow
  otherwise, as standalone JULES does.

With forcing on, the JULES explicit and implicit kernels use the namelist
forcing heights for ``z1_tq`` and ``z1_uv``, and set ``z1_tq_top`` and
``z1_uv_top`` to zero. ``jules_physics_init`` sets ``i_modiscopt`` to 0,
because the forcing is at a fixed height rather than a model level.

Radiation
~~~~~~~~~

No radiation scheme is run. Instead, ``jules_forcing_init`` (called from
``gungho_model_mod``) sets the SW bands normally read by ``socrates_init``.
It sets one broadband band (wavelength limits of -1, which is
``spec_band_bb`` in standalone JULES), with ``sw_weight_blue`` set to
``visible_fraction``. With this setup, the existing ``sw_rad_tile`` call in
``rad_tile_alg`` calculates broadband tile albedos weighted equally between
visible and near-infrared. These are the default ``wght_alb`` values in
standalone JULES. After ``rad_tile_alg``, ``jules_forcing_sw_up_kernel`` sets
the upward SW flux on each tile to the diffuse albedo multiplied by
``SWdown``.

The downward LW comes straight from the forcing. The tile emissivities come
from the existing LW tile calculation.

The solar zenith angle comes from ``illuminate_alg``, which is called in
``jules_timestep_alg`` each timestep. The radiation, orbit and star namelists
are not needed: ``illuminate`` uses default orbital parameters and a
radiation step of one when the radiation scheme is not Socrates.

Initial state and surface properties
------------------------------------

There are no ancillary files. The surface is set from the ``ideal_surface``
namelist, and the soil from the analytic branch of ``init_soil_fields_alg``.
To support site runs, these ``ideal_surface`` items were added. They were
previously hard-coded:

* soil hydraulic and thermal properties: ``clapp_horn_b``, ``soil_albedo``,
  ``soil_cond_sat``, ``soil_moist_crit``, ``soil_moist_sat``,
  ``soil_moist_wilt``, ``soil_suction_sat``, ``soil_thermal_cap`` and
  ``soil_thermal_cond``;
* the TOPMODEL topographic index mean and standard deviation,
  ``mean_topog_index`` and ``stdev_topog_index``;
* the initial ``water_table`` and ``wetness_under_soil``.

The ``vn32`` upgrade macros in ``jules-lsm`` and ``lfric-gungho`` set the
new items to the old hard-coded values. This leaves existing configurations
unchanged.

The atmospheric profile namelists only use ``profile_size=1`` placeholder
values, because the forcing overwrites the first layer before JULES runs.
``map_physics_fields_alg`` is skipped when forcing is on.

Configuration
-------------

Required namelists
~~~~~~~~~~~~~~~~~~

``applications/jules/source/jules_mod.f90`` defines
``jules_required_namelists``, which ``jules.f90`` uses in place of the Gungho
list. Only these namelists are required: ``finite_element``,
``formulation``, ``base_mesh``, ``initial_wind``, ``planet``,
``timestepping``, ``extrusion``, ``orography`` and ``energy_correction``.

The ``energy_correction`` namelist is required because a missing
``encorr_usage`` would switch energy correction on. That then fails on a
planar mesh.

The ``lfric-jules`` metadata imports ``lfric-gungho`` and makes the unused
namelists optional (``compulsory=false``). These are ``boundaries``,
``checks``, ``convection``, ``cosp``, ``departure_points``, ``esm_couple``,
``helmholtz_solver``, ``iau``, ``initial_density``, ``mixed_solver``,
``mixing``, ``solver``, ``specified_surface`` and ``transport``. The
transport checks in ``check_configuration_mod`` only run when the
``transport`` namelist has been loaded. The model is run with
``boundary_layer='none'``.

The canned example
~~~~~~~~~~~~~~~~~~

``applications/jules/example`` has a ``configuration.nml`` and
``iodef.xml`` for the Cardington site (52.103 N, 359.576 E). The setup is:

* start time 2005-01-01 00:00;
* 143 timesteps of 1800 s;
* 100% C3 grass;
* 10 m forcing heights.

The initial state, soil properties and topographic index match the
standalone JULES suite u-dz368 (RAL3). Before running, set
``jules_forcing_path`` to the location of the Cardington forcing file.

Rose stem
~~~~~~~~~

The rose-stem app is ``rose-stem/app/jules``. The Cardington test runs as
the task ``cardington-1x1-50000x50000`` in ``tasks_jules.cylc``. It is
included in the group
``jules_cardington-1x1-50000x50000_azspice_gnu_fast-debug-64bit``. The main
app is kept consistent with ``lfric_atm``. Three optional configurations
apply, in this order:

``rose-app-1x1-50000x50000.conf``
   Mesh and ``base_mesh`` settings.

``rose-app-ral4.conf``
   The JULES settings from the ``lfric_atm`` RAL3 optional configuration
   (``rose-stem/app/lfric_atm/opt/rose-app-ral3.conf``). It covers the lake,
   plant functional type and snow parameters, plus the ``jules_surface``
   switches. The ``lfric_atm`` settings for other schemes are left out,
   because those schemes are not run in this app.

``rose-app-cardington.conf``
   The settings where u-dz368 differs from ``ral4``.
   It is applied after ``ral4`` and overrides it. The settings are:

   * plant functional type albedo parameters (``alnir``, ``alpar``,
     ``omega``, ``omnir``);
   * shrub roughness (``z0hm_pft_io=0.05``, ``z0v_io=0.2``), which matches
     the RAL3 configuration of u-dz368;
   * snow canopy parameters (``can_clump`` and ``n_lai_exposed``);
   * ``formdrag='none'``;
   * ``l_mo_buoyancy_calc``, ``l_var_rainfrac`` and ``l_sea_alb_var_chl``
     switched off;
   * a CO\ :sub:`2` mixing ratio of 5.941e-4.

The ``[file:configuration.nml]`` source list holds every namelist, including
the ones that are ignored. Optional ones are marked ``(optional)``. This means
a namelist that is trigger-ignored is still output if a setting brings it
back into use.

Comparison with standalone JULES
--------------------------------

The Cardington test has been compared with the equivalent standalone JULES
run (u-dz368) for January 2005. The two agree closely:

* soil moisture, soil temperature and runoff (including TOPMODEL baseflow)
  are effectively identical;
* surface energy-balance terms agree to within about 0.3 W m\ :sup:`-2`.

There is one known remaining difference. Soil evaporation is about 2% higher
in LFRic, which offsets the latent and sensible heat fluxes by
+0.3 W m\ :sup:`-2` and -0.2 W m\ :sup:`-2`. This difference is present from
the first timestep, and its cause has not yet been found. The 1.5 m specific
humidity is also slightly lower, by about 0.02 g kg\ :sup:`-1`. This is
expected, because the screen diagnostic method (``iscrntdiag``) cannot be
chosen in LFRic.

Limitations
-----------

* Only single-point runs on one MPI rank are supported.
* The model timestep must equal the interval of the forcing data, and the
  ``output_freq`` of the ``jules_forcing`` file must match it. The data is
  not interpolated in time.
* All precipitation is supplied as ``Rainf`` and split into rain and snow at
  274 K. Separate rain, snow, or convective and large-scale precipitation
  are not supported.
* Some standalone JULES options cannot be set in LFRic. These include
  ``iscrntdiag`` and the urban schemes (for example ``l_urban2t``).
* Diffuse fraction and visible fraction are fixed at 0.4 and 0.5.
* Some JULES inputs are deliberately left unset, because the LFRic build of
  JULES does not use them. Examples are ``ti_sig_soilt``, ``a_step``, and the
  momentum terms of ``surf_couple_implicit``.
* Many fields that are only used by the rest of the atmosphere model are not
  initialised. They should not be written as diagnostics.
