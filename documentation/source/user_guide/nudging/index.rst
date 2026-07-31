.. -----------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _nudging_user_index:

Nudging
=======

Nudging is a form of external forcing that relaxes model fields towards a
specified reference state, typically taken from analysis or reanalysis data.
It can be used, for example, to keep a forecast close to observed large-scale
conditions, or to spin up a model state consistent with an external data
source.

A full mathematical description of the nudging scheme is given in the
:ref:`Science Guide <nudging_science_index>`. This page describes how to
configure and use nudging.

Enabling nudging
----------------

Nudging is enabled per-field via the ``namelist:external_forcing`` options
``theta_forcing`` and ``wind_forcing``. Setting either (or both) of these to
``'nudging'`` activates the nudging scheme for potential temperature or the
horizontal wind components respectively:

.. code-block:: ini

    [namelist:external_forcing]
    theta_forcing='nudging'
    wind_forcing='nudging'

The ``namelist:nudging`` namelist must then be configured as described
below.

Reference data files
---------------------

The external reference data used for nudging is read from a file specified
in the ``namelist:files`` namelist:

.. list-table::
    :header-rows: 1
    :widths: 25 75

    * - Option
      - Description
    * - ``nudging_directory``
      - Path to the directory containing the nudging reference data file.
    * - ``nudging_filename``
      - Name of the file (within ``nudging_directory``) containing the
        reference temperature, wind and surface pressure fields.

Configuring the nudging scheme
-------------------------------

The behaviour of the nudging scheme is controlled by the
``namelist:nudging`` namelist:

.. list-table::
    :header-rows: 1
    :widths: 25 75

    * - Option
      - Description
    * - ``nudging_method``
      - Selects the nudging algorithm: ``'newtonian'`` for pointwise
        relaxation towards the reference state, or ``'convolution'`` to
        relax only the large scales of the field, using a convolution in
        physical space that emulates a low-pass top-hat filter in spectral
        space (see the :ref:`Science Guide <nudging_science_index>`).
    * - ``nudging_relax_time``
      - Timescale, in seconds, over which the field (or its large-scale
        component) is relaxed towards the reference state. A value of zero
        sets the field to the reference state in a single time step.
    * - ``nudging_spinup_start`` / ``nudging_spinup_end``
      - Times, in seconds since the start of the run, over which the
        strength of nudging is ramped up linearly from zero to full
        strength. This allows nudging to be introduced gradually rather
        than as a step change.
    * - ``num_ref_data_levels``
      - Number of vertical levels on which the reference data for temperature
        and wind components are provided.
    * - ``nudging_level_bottom`` / ``nudging_level_top``
      - Model levels defining the vertical extent over which nudging is
        applied at full strength.
    * - ``nudging_level_taper_bottom`` / ``nudging_level_taper_top``
      - Model levels below/above which nudging is zero. Nudging is tapered
        linearly between ``nudging_level_taper_bottom`` and
        ``nudging_level_bottom``, and between ``nudging_level_top`` and
        ``nudging_level_taper_top`` (see the :ref:`Science Guide
        <nudging_science_vertical_taper>`).
    * - ``spectral_kmin`` / ``spectral_kmax``
      - (``nudging_method = 'convolution'`` only) Minimum and maximum
        wavenumbers retained by the top-hat spectral filter.
    * - ``spectral_stencil_extent``
      - (``nudging_method = 'convolution'`` only) Extent of the stencil
        used to construct the physical-space convolution.
    * - ``spectral_envelope_width``
      - (``nudging_method = 'convolution'`` only) Width, in radians, of the
        Gaussian envelope applied to the convolution filter (see the
        :ref:`Science Guide <nudging_science_convolution>`).

Choosing the stencil extent and envelope width
-----------------------------------------------

The ``spectral_stencil_extent`` and ``spectral_envelope_width`` options must
be chosen together, and some care is needed:

* The Legendre-polynomial filter that underlies the convolution does not
  have compact support, so truncating it to a finite stencil is itself an
  approximation. The Gaussian envelope is what makes this truncation
  acceptable, by bringing the filter smoothly down to a small value by the
  edge of the stencil.
* If ``spectral_envelope_width`` is too large relative to
  ``spectral_stencil_extent``, the filter has not decayed sufficiently by
  the edge of the stencil. The abrupt truncation this causes can amplify,
  rather than damp, some wavenumbers, leading to noise or numerical
  instability rather than the intended scale-selective nudging.
* Conversely, if ``spectral_envelope_width`` is too small, the envelope
  suppresses part of the range of wavenumbers between
  ``spectral_kmin`` and ``spectral_kmax`` that was intended to be retained,
  making the filter less selective than requested.
* Increasing ``spectral_stencil_extent`` allows a larger
  ``spectral_envelope_width`` to be used safely, but at the cost of a wider
  halo exchange and increased computational expense.

Because of this trade-off, it is recommended to test any new combination of
``spectral_kmin``, ``spectral_kmax``, ``spectral_stencil_extent`` and
``spectral_envelope_width`` for a given mesh resolution, checking that the
resulting weights decay smoothly across the stencil, before using it in
a forecast or research configuration.

Choosing the vertical extent of nudging
-----------------------------------------

The four vertical-level options must satisfy
``nudging_level_taper_bottom`` :math:`\le` ``nudging_level_bottom``
:math:`\le` ``nudging_level_top`` :math:`\le` ``nudging_level_taper_top``;
this ordering is enforced automatically when the namelist is checked.
Nudging is at full strength between ``nudging_level_bottom`` and
``nudging_level_top``, tapers linearly to zero outside this range, and is
exactly zero at and beyond ``nudging_level_taper_bottom`` /
``nudging_level_taper_top``. Setting a taper level equal to the
corresponding full-strength level removes the ramp on that side, giving an
abrupt step instead of a gradual transition.

Note that these levels are counted differently for the two sets of model
levels that nudging acts on: levels holding potential temperature are counted
from 0 at the surface, while levels holding density or the horizontal wind
components are counted from 1/2 in the lowest level, since they are
vertically staggered by half a level relative to the potential temperature
levels (see the :ref:`Science Guide <nudging_science_vertical_taper>`). The
same integer level values are used for both, so the two sets of levels see
a consistent, vertically-aligned taper profile.

Coarse-mesh nudging
--------------------

For efficiency, nudging may be computed on a mesh coarser than
the one used by the dynamical core, with the resulting increment mapped
back onto the dynamics mesh (see the :ref:`Science Guide
<nudging_science_multiresolution>`).

Coarse-mesh nudging is controlled from the ``namelist:multires_coupling``
namelist:

.. list-table::
    :header-rows: 1
    :widths: 25 75

    * - Option
      - Description
    * - ``coarse_nudging``
      - If ``.true.``, nudging increments are computed on a coarser mesh and
        then mapped back onto the dynamics mesh. If ``.false.`` (the
        default), nudging is computed directly on the dynamics mesh and no
        coarse mesh needs to be specified.
    * - ``nudging_mesh_name``
      - Tag-name of the coarser mesh to use when ``coarse_nudging`` is
        ``.true.``. This mesh must also be listed in
        ``multires_coupling_mesh_tags``.
