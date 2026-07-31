.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _nudging_science_formulation:

Formulation
===========

Let :math:`X` denote a field being nudged (potential temperature
:math:`\theta`, or one of the horizontal wind components :math:`u`, :math:`v`),
and let :math:`X_{ref}` denote the corresponding reference field, having
already been vertically interpolated onto the model's levels (see
:ref:`nudging_science_vertical_treatment`). At each call to the nudging
scheme, an increment :math:`\Delta X` is calculated and added to the
model's tendency for :math:`X`.

Two methods are available for calculating :math:`\Delta X`.

.. _nudging_science_formulation_newtonian:

Newtonian relaxation
---------------------

In Newtonian relaxation, the increment is simply proportional to the local
difference between the field and its reference value:

.. math:: :label: eq:nudging_newtonian

   \Delta X = w \left( X_{ref} - X \right)

where :math:`w` is a weight described in :ref:`nudging_science_weights`
below. This relaxes :math:`X` directly towards :math:`X_{ref}` at every
horizontal scale present in the reference data, with no spectral
selectivity.

.. _nudging_science_formulation_convolution:

Convolution
-----------

In the convolution method, the local difference between the field and its
reference value is first filtered so that only the large scales are
retained, before being scaled by the same weight :math:`w`:

.. math:: :label: eq:nudging_convolution

   \Delta X = w \, \mathcal{L}\left[ X_{ref} - X \right]

where :math:`\mathcal{L}` is a scale-selective filter that removes the
short wavelengths. Rather than performing this filtering with a global
spectral transform, :math:`\mathcal{L}` is applied as a convolution in
physical space. The motivation for this choice, and the derivation of the
convolution kernel on the sphere, are described in
:ref:`nudging_science_convolution`.

.. _nudging_science_weights:

Nudging weights
----------------

In both of the above methods, the increment is scaled by a single weight
:math:`w`, which combines a vertical tapering factor, a spin-up ramp, and a
relaxation-timescale factor:

.. math:: :label: eq:nudging_weight

   w = w_{taper} \; r_{spinup} \; \min\left(1, \frac{\Delta t}{\tau}\right)

where:

* :math:`w_{taper}` is a height-dependent factor, computed as described in
  :ref:`nudging_science_vertical_treatment`, that ramps the nudging on and
  off over a chosen range of vertical levels;
* :math:`r_{spinup}` is a factor that ramps linearly from 0 to 1 over a
  configured spin-up period, described below;
* :math:`\Delta t` is the model timestep and :math:`\tau` is a configured
  relaxation timescale. If :math:`\tau` is set to zero, this factor is
  instead set to 1, so that the field is relaxed fully towards the
  reference state within a single timestep.

Spin-up ramp
~~~~~~~~~~~~

The nudging weight can be ramped up gradually over a configurable spin-up
period :math:`[t_{start}, t_{end}]`:, allowing forecasts over short time scales
to be unperturbed.

.. math:: :label: eq:nudging_spinup

   r_{spinup} =
   \begin{cases}
     0 & t \le t_{start} \\
     \dfrac{t - t_{start}}{t_{end} - t_{start}} & t_{start} < t < t_{end} \\
     1 & t \ge t_{end}
   \end{cases}

where :math:`t` is the elapsed model time since the start of the run. If
:math:`t_{start}` and :math:`t_{end}` are equal, the ramp is skipped and
:math:`r_{spinup} = 1` once :math:`t` exceeds this value.
