.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.

    Some of the content of this file has been produced with the assistance of
    Met Office GitHub Copilot Enterprise.
   -----------------------------------------------------------------------------
.. _nudging_science_vertical_treatment:

Vertical Treatment
===================

The external reference data used for nudging may be supplied on its own set of
hybrid-pressure levels, which do not in general coincide with the model's
levels. Before the reference fields can be used in the nudging increment,
they must therefore be interpolated onto the model's levels, and the
temperature field converted to potential temperature. A height-dependent
weight is also computed, to allow nudging to be restricted to a chosen
range of vertical levels.

Reference pressure profile
---------------------------

The reference dataset's hybrid-pressure levels are defined by coefficients
:math:`a_k`, :math:`b_k`. Given the reference surface pressure,
:math:`p_{s,ref}`, the pressure on reference half-level :math:`k` is

.. math:: :label: eq:nudging_ref_pressure

   p_{ref,k} = a_k + b_k \, p_{s,ref} ,

and full-level values are formed by averaging :math:`a_k` and :math:`b_k`
between adjacent half levels.

Model pressure profile
------------------------

The corresponding pressure on the model's levels is obtained from the
model's Exner pressure, :math:`\Pi`, evaluated on both sets of model
levels (those on which pressure is held, and those on which potential
temperature is held):

.. math:: :label: eq:nudging_model_pressure

   p = p_0 \, \Pi^{1/\kappa}

where :math:`p_0` is a reference pressure and
:math:`\kappa = R_d / c_p`, with :math:`R_d` the gas constant for dry air
and :math:`c_p` the specific heat capacity at constant pressure.

Temperature to potential temperature
--------------------------------------

The reference dataset supplies temperature, :math:`T_{ref}`, rather than
potential temperature. This is converted using the reference pressure
profile computed above:

.. math:: :label: eq:nudging_T_to_theta

   \theta_{ref} = T_{ref} \left( \frac{p_0}{p_{ref}} \right)^{\kappa} .

Vertical interpolation
------------------------

The reference fields (:math:`\theta_{ref}`, :math:`u_{ref}`,
:math:`v_{ref}`) are interpolated from the reference pressure levels onto
the model's pressure levels by linear interpolation in log-pressure:

.. math:: :label: eq:nudging_log_p_interp

   X_{ref}(p) = X_{ref,l} + \left( X_{ref,l+1} - X_{ref,l} \right)
   \frac{\ln p - \ln p_{ref,l}}{\ln p_{ref,l+1} - \ln p_{ref,l}}

where :math:`l` is the reference level immediately below the model level
being interpolated to (in pressure). To avoid unbounded extrapolation
above the top or below the bottom of the reference data, the interpolation
factor is limited to the range :math:`[-1, 2]`. The lowest model
:math:`\theta` level is set equal to the value interpolated for the level
immediately above it.

.. _nudging_science_vertical_taper:

Vertical tapering
-------------------

To allow nudging to be targeted at a chosen range of vertical levels, a
height-dependent weight, :math:`w_{taper}` (introduced in
:eq:`eq:nudging_weight`), is computed for each function space
(:math:`\mathbb{W}_3` and :math:`\mathbb{W}_\theta`). The weight ramps
linearly from 0 to 1 between a lower taper level and a lower full-strength
level, remains at 1 between the full-strength levels, and ramps back down
to 0 between an upper full-strength level and an upper taper level:

.. math:: :label: eq:nudging_vertical_taper

   w_{taper}(\tilde{k}) =
   \begin{cases}
     0 &
       \tilde{k} \le k_{taper,bot} \\
     \dfrac{\tilde{k} - k_{taper,bot}}{k_{bot} - k_{taper,bot}} &
       k_{taper,bot} < \tilde{k} < k_{bot} \\
     1 &
       k_{bot} \le \tilde{k} \le k_{top} \\
     1 - \dfrac{\tilde{k} - k_{top}}{k_{taper,top} - k_{top}} &
       k_{top} < \tilde{k} < k_{taper,top} \\
     0 &
       \tilde{k} \ge k_{taper,top}
   \end{cases}

where :math:`k_{bot}` and :math:`k_{top}` are the model levels marking the
bottom and top of the region over which nudging is applied at full
strength, and :math:`k_{taper,bot}` and :math:`k_{taper,top}` are the
levels below and above which nudging is zero. The quantity
:math:`\tilde{k}` is the :math:`\mathbb{W}_\theta` level index starting from 0,
while :math:`\mathbb{W}_3` points are offset by half a level. Using this common
coordinate ensures that both sets of levels see a consistent, vertically-aligned
taper profile.

.. _nudging_science_tropopause_cap:

Tropopause capping
-------------------

Nudging towards external reference data is generally undesirable well above
the tropopause: the stratosphere is only loosely constrained (if at all) by
the coarse vertical and temporal sampling of typical reference datasets, and
imposing large-scale forcing there can disrupt the model's own
stratospheric transport and composition. Rather than relying solely on a
fixed upper taper level, chosen in advance to sit safely below the
tropopause for an expected range of conditions, the upper extent of the
tapering profile described above may instead be capped dynamically by the
model's own, evolving tropopause height.

When this option is used, the upper full-strength and taper levels,
:math:`k_{top}` and :math:`k_{taper,top}`, in :eq:`eq:nudging_vertical_taper`
are replaced at each application of the scheme by the lesser of their
configured values and the current tropopause level, diagnosed
independently for each atmospheric column from the model's evolving
temperature structure. This allows the vertical extent of nudging to track
a rising or falling tropopause, rather than remaining fixed, and ensures
nudging is switched off in the stratosphere even where the diagnosed
tropopause sits below the configured upper taper level. A fixed minimum
level may also be enforced as a floor for the diagnosed tropopause level,
guarding against spuriously low diagnoses; setting this floor high enough
disables the dynamic capping altogether, recovering the fixed taper profile
described above.

