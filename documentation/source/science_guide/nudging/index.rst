.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.

    Some of the content of this file has been produced with the assistance of
    Met Office GitHub Copilot Enterprise.
   -----------------------------------------------------------------------------
.. _nudging_science_index:

Nudging
================

Nudging is a technique for relaxing a free-running simulation towards an
external reference state. It can be targeted at a chosen
range of horizontal scales and a chosen range of vertical levels. This
allows the large-scale state of the simulation to be kept close to the
reference state, while still permitting the model to generate its own
small-scale variability. Two approaches are available for this nudging:
selecting only certain horizontal scales (known as spectral nudging), or
relaxing the model state directly towards the reference state at every
horizontal scale (known as Newtonian relaxation).

This section describes the scientific formulation of nudging in GungHo. It
begins with an overview of the method and the relaxation increment, before
describing how the large scales can be selected using a convolution, how this
may be computed efficiently on a coarse mesh, and finally how the reference
data is prepared in the vertical. For details of how to configure and
enable nudging, see the :ref:`nudging_user_index` section of the User
Guide.

.. toctree::
    :maxdepth: 1

    overview
    formulation
    convolution
    multiresolution
    vertical_treatment
