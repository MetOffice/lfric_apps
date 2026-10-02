.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _lipschitz_damping_doc:

=================
Lipschitz Damping
=================

Introduction
============

Lipschitz damping is a targeted scheme for damping a wind field wherever it is
locally unstable, as measured by its Lipschitz number. It is applied to the
wind increment from fast physics, before this increment is added to the
dynamics wind field, and may optionally also be applied to the initial wind
field.

Lipschitz numbers
==================

The Lipschitz numbers measure the divergence of a particular component of the
wind field, so that for the i-th direction:

.. math:: :label: lipschitz_definition

   \ell_i = \Delta t \nabla_i \cdot u_i


Let the computational :math:`\mathbb{W}_2` wind field (with values on the
faces, scaled by area) be given as :math:`u(E)`, :math:`u(W)`, :math:`u(N)`,
:math:`u(S)`, :math:`u(T)` and :math:`u(B)` on its east, west, north, south,
top and bottom faces respectively.
For a cell with volume :math:`V`, timestep :math:`\Delta t`, the three 1D
Lipschitz numbers can be simply computed as

.. math:: :label: lipschitz_1d

   \ell_x = \left( u(E) - u(W) \right) \frac{\Delta t}{V}, \quad
   \ell_y = \left( u(S) - u(N) \right) \frac{\Delta t}{V}, \quad
   \ell_z = \left( u(T) - u(B) \right) \frac{\Delta t}{V}

and the 3D Lipschitz number is their sum,

.. math:: :label: lipschitz_3d

   \ell_{3D} = \ell_x + \ell_y + \ell_z


These give a measure of whether the flow's divergence in that cell, over one
timestep, is large enough to cause numerical instability for mass-conserving
transport schemes. Positive values larger than 1 correspond to outflow that
can generate a negative mass in a cell. Negative values correspond to inflow
so do not directly cause this instability. Therefore this damping algorithm
only acts upon positive values that exceed a threshold.

Damping algorithm
==================

For each cell, the following steps are applied:

1. For each of the three 1D Lipschitz numbers, if it exceeds a threshold of
   1, the outflowing wind component(s) responsible are scaled back so that
   the Lipschitz number falls to exactly 1. Inflowing components are left
   untouched. Where both components on a pair of opposite faces are
   outflowing, the correction is shared between them in proportion to their
   outflow magnitude.
2. The 3D Lipschitz number is then recomputed from the (possibly already
   damped) wind components. If it exceeds a threshold of 0.5, the same
   procedure is applied across all six outflowing components to bring it
   back down to 0.5.

Because only the outflowing side of a face is ever adjusted, neighbouring
cells cannot both attempt to modify the same shared wind component, so the
scheme can be applied independently, cell by cell.