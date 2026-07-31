.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _nudging_science_overview:

Overview
========

Spectral nudging relaxes the potential temperature, :math:`\theta`, and the
horizontal wind components, :math:`(u,v)`, towards reference fields derived
from an external dataset (for example a reanalysis, or the driving model in
a regional or limited-area configuration).

At each application of the scheme, an increment is calculated for a given
field :math:`X \in \{\theta, u, v\}` and added to the model's prognostic
tendencies. Two methods are available for calculating this increment:

* **Newtonian relaxation:** a pointwise relaxation of the field towards the
  reference field, described in
  :ref:`nudging_science_formulation_newtonian`.
* **Convolution:** a physical-space convolution that emulates a spectral
  low-pass (or band-pass) filter, so that only the large scales are nudged.
  The rationale for this approach, and its derivation on the sphere, are
  described in :ref:`nudging_science_convolution`.

In both cases the increment is scaled by a set of weights which:

* ramp the nudging increment smoothly up and down with height, targeting
  nudging at a chosen range of vertical levels (see
  :ref:`nudging_science_vertical_treatment`);
* ramp the nudging on gradually at the start of a run, over a configurable
  spin-up period;
* scale the increment according to the model timestep and a configurable
  relaxation timescale.

These weights are described fully in :ref:`nudging_science_formulation`.

The external reference data is supplied on its own set of pressure levels,
and must be vertically interpolated onto the model's levels before it can be
used; this process is described in
:ref:`nudging_science_vertical_treatment`.

Finally, when the convolution method is used, nudging can optionally be
performed on a coarser mesh than the mesh used for the rest of the model,
with the resulting increment prolongated back onto the model's mesh. This
improves efficiency and is described in
:ref:`nudging_science_multiresolution`.
