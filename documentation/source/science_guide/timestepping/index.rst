.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
    Some of the content of this file has been produced with the assistance of
    Met Office GitHub Copilot Enterprise.
   -----------------------------------------------------------------------------

.. _science_guide_timestepping:

Timestepping
============

GungHo's dynamical core can be advanced in time using one of several
timestepping schemes (see the :ref:`User Guide <user_guide_timestepping>` for
how these are configured). This section describes
the numerical formulation of each scheme:

.. toctree::
    :maxdepth: 2

    semi_implicit
    runge_kutta
    tr_bdf2

.. attention::

   The formulations below use the notation of a state vector
   :math:`\boldsymbol{X}`, an operator :math:`\mathcal{F}` representing the
   evaluation of forces (such as the pressure gradient, Coriolis and
   gravitational forces), and a transport operator
   :math:`\mathcal{T}^{\Delta t}_{\boldsymbol{u}}`, relating to a velocity
   :math:`\boldsymbol{u}`, over a time step :math:`\Delta t`. The transport
   operator returns a transported field (rather than a transport increment).
   Increments from physics parametrisation schemes are represented by
   :math:`\mathcal{P}`.
