.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _nudging_science_multiresolution:

Multi-resolution Nudging
==========================

When the convolution method is used, the nudging increment can optionally be
computed on a mesh that is coarser than the one used for the rest of the
model. Since the nudging is usually desired only on large scales, this gives
similar scientific behaviour but at considerably reduced cost..

The convolution (see :ref:`nudging_science_convolution`) must span the range
of scales that are being nudged. Nudging is typically applied at horizontal
scales of the order of :math:`2000\,\text{km}`, which at a global
NWP resolution corresponds to a stencil spanning many cells and, in a
distributed implementation, many processor partitions. Computing the
convolution directly on the model mesh would therefore require a very large
halo depth and correspondingly large halo exchanges, which could be as
expensive as the global spectral transforms the convolution is intended to
avoid.

The cost is greatly reduced by applying the convolution only to a
coarse-grained form of the field. If the coarse mesh has of the order of a
single cell per partition, then nudging over a given physical distance
requires halo exchanges of only a modest depth, which is far more efficient.

Coarse-grained convolution
--------------------------

To describe this, let :math:`\mathcal{A}` be a restriction operator that
performs a local averaging to compute a coarse-grained field, and let
:math:`\mathcal{R}` be the corresponding prolongation operator, which
reconstructs a field on the fine mesh from a coarse-grained field using a
linear reconstruction. The filter :math:`\mathcal{L}` of
:eq:`eq:nudging_convolution` is then applied as

.. math:: :label: eq:nudging_coarse_convolution

   \mathcal{L}[X_{ref} - X] =
   \mathcal{R}\Big[\, \overline{\mathcal{C}}\big[\,
   \mathcal{A}[X_{ref} - X], \; \overline{c} \,\big] \Big] ,

where :math:`\overline{\mathcal{C}}` denotes the convolution computed on the
coarse mesh, and :math:`\overline{c}` is the corresponding coarse-mesh
kernel. In words, the increment field is first restricted to the coarse
mesh, the convolution is performed there, and the result is prolongated back
onto the model's mesh, where it is applied to the model's tendency for
:math:`X`.

When coarse-mesh nudging is not used, no restriction or prolongation takes
place, and the convolution is computed directly on the model's own mesh.

The external reference fields used to compute the increment (see
:ref:`nudging_science_vertical_treatment`) are themselves prepared on the
coarse mesh when coarse-mesh nudging is used: in particular, the pressure
fields used to construct the model pressure profile are first mapped onto
the coarse mesh before being used in the vertical interpolation.

For details of how to enable and configure coarse-mesh nudging, see the
:ref:`nudging_user_index` section of the User Guide.
