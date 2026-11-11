.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
    Some of the content of this file has been produced with the assistance of
    Met Office GitHub Copilot Enterprise.
   -----------------------------------------------------------------------------

.. _science_guide_timestepping_runge_kutta:

Runge-Kutta
===========

The Runge-Kutta timestepping method is an explicit,
multi-stage scheme. It can be used either to timestep
the whole dynamical core, or in combination with the Method-of-Lines
transport scheme. The Strong Stability Preserving (SSP) variants
are explicit multi-stage methods
which combine the values computed over multiple stages as a weighted linear
sum, in order to provide stability. A simple forward Euler option is also
available.

Runge-Kutta timestepping of the whole dynamical core has not been implemented
in combination with physics parametrization schemes.

Rather than splitting terms as in the Quasi-Newton schemes, the Runge-Kutta
treats the whole equation set for the state vector :math:`\boldsymbol{X}` as one:

.. math::

   \frac{\partial\boldsymbol{X}}{\partial t} + \mathcal{G}(\boldsymbol{X}) = \boldsymbol{0},

where :math:`\mathcal{G}` encompasses forcing and transport terms.
An explicit Runge-Kutta scheme advances one timestep by evaluating a sequence of
:math:`M` stages,

.. math::

   \begin{aligned}
   \boldsymbol{X}^{(0)} &= \boldsymbol{X}^{n}, \\
   \boldsymbol{X}^{(m)} &= \boldsymbol{X}^{(0)}
   - \Delta t\sum_{i=1}^{m} c_{m,i}\mathcal{G}[\boldsymbol{X}^{(i-1)}], \qquad \forall m = 1, \ldots, M, \\
   \boldsymbol{X}^{n+1} &= \boldsymbol{X}^{(M)}.
   \end{aligned}

The scheme is explicit because each stage :math:`i` depends only on
previous stages so no implicit solve is required.
In practice this means a
Runge-Kutta step is a sequence of tendency evaluations and weighted sums,
with higher-order SSP variants using more stages than Forward Euler to
improve stability and accuracy.

In GungHo, this same stage-update pattern is used whether Runge-Kutta is
applied to the whole dynamical core, or used inside Method-of-Lines
transport where :math:`\mathcal{G}` is the transport tendency operator.
