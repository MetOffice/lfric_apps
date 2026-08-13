.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
    Some of the content of this file has been produced with the assistance of
    Met Office GitHub Copilot Enterprise.
   -----------------------------------------------------------------------------

.. _science_guide_timestepping_tr_bdf2:

TR-BDF2
=======

Overview
--------

The Trapezoidal Backwards-Difference-2 (TR-BDF2) timestepper
is a two-stage semi-implicit method, using a
Quasi-Newton iterative structure like the :ref:`Semi-Implicit
<science_guide_timestepping_semi_implicit>` scheme. The Trapezoidal (TR) and
Backwards Difference (BDF2) stages each take the form of a Semi-Implicit
step. Unlike the Semi-Implicit Quasi-Newton (SIQN) scheme, TR-BDF2 has no
off-centering parameter.

A TR-BDF2 timestepper may have several advantages over the SIQN scheme:

* **Numerical stability**: the SIQN scheme does not damp wave modes unless it
  uses an implicit off-centering. The solution can then become polluted by
  high-frequency waves (typically acoustic modes) triggered by small-scale
  features such as those from physics parametrizations or orography. The
  growth of these waves can cause numerical instabilities. In contrast, the
  TR-BDF2 scheme is inherently damping of high-frequency waves, and in fact
  is more damping of the highest frequencies than an off-centred SIQN scheme.
* **Accuracy**: the SIQN scheme requires implicit off-centering for
  stability, but this off-centering reduces the formal order of accuracy of
  the scheme to first order in time. The TR-BDF2 scheme has no off-centering
  parameter, and damps high-frequency waves while keeping second-order
  accuracy.
* **Computational efficiency**: the TR-BDF2 scheme is a three-level time
  discretisation, compared with the two-level discretisation of the SIQN
  scheme. A single timestep with the TR-BDF2 scheme involves approximately
  twice as much computational work as a SIQN timestep, so it is natural to
  take twice the timestep length with TR-BDF2. The TR-BDF2 formulation also
  provides the mathematical basis for calling some schemes (physics
  parametrizations, or the transport of some variables) less frequently,
  reducing the computational cost. Further, the improved stability of the
  scheme may allow it to run stably with fewer Quasi-Newton iterations,
  again improving the computational efficiency.

The inspiration for much of this formulation comes from [Tumolo2015]_,
who presented a semi-Lagrangian formulation of the compressible Euler
equations using a TR-BDF2 timestepper with Discontinuous Galerkin spatial
discretisations.

Details
-------

Basic formulation
~~~~~~~~~~~~~~~~~~

The TR-BDF2 scheme is a three-level scheme with a "midpoint" level denoted by
:math:`m`. Using the off-centering-free implicit weight :math:`\gamma`, it
follows:

.. math:: :label: trbdf2_tr_stage

   \boldsymbol{X}^{m} + \gamma\Delta t \mathcal{F}(\boldsymbol{X}^m)
   = \mathcal{T}^{2\gamma\Delta t}_{\overline{\boldsymbol{u}}_n^m}
   \left[\boldsymbol{X}^n - \gamma\Delta t \mathcal{F}(\boldsymbol{X}^n)\right],

.. math:: :label: trbdf2_bdf2_stage

   \boldsymbol{X}^{n+1} + \gamma_2\Delta t \mathcal{F}(\boldsymbol{X}^{n+1})
   = (1-\gamma_3)\mathcal{T}^{\Delta t}_{\overline{\boldsymbol{u}}_n^{n+1}}
   \left[\boldsymbol{X}^n\right]
   + \gamma_3\mathcal{T}^{(1-2\gamma)\Delta t}_{\overline{\boldsymbol{u}}_m^{n+1}}
   \left[\boldsymbol{X}^m\right],

with transporting velocities:

.. math::

   \overline{\boldsymbol{u}}_n^{m} = \tfrac{1}{2}\left(\boldsymbol{u}^n + \boldsymbol{u}^{m}\right), \qquad
   \overline{\boldsymbol{u}}_n^{n+1} = \boldsymbol{u}^{n+1}, \qquad
   \overline{\boldsymbol{u}}_m^{n+1} = \boldsymbol{u}^{n+1}.

The parameters :math:`\gamma_2` and :math:`\gamma_3` are given by

.. math::

   \gamma_2 = \frac{1-2\gamma}{2(1-\gamma)}, \qquad
   \gamma_3 = \frac{1-\gamma_2}{2\gamma}.

GungHo takes :math:`\gamma=1-\sqrt{2}/2\approx 0.293`, which gives
:math:`\gamma_2=1-\sqrt{2}/2\approx 0.293` and
:math:`\gamma_3=(1+\sqrt{2})/2\approx 1.207`. The TR stage has an effective
timestep length of :math:`2\gamma\Delta t\approx 0.586\Delta t`, while the
BDF2 stage has an effective timestep length of
:math:`(1-2\gamma)\Delta t \approx 0.414\Delta t`.
Both stages take the same basic form as the Semi-Implicit scheme.


Combining BDF2 Transport Steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both stages of the TR-BDF2 scheme involve nested outer and inner iterative
loops (as for the SIQN scheme). If these follow the SIQN scheme, there would
be 2 iterations each of the outer and inner loops in both stages. The
transport operator is called in the outer loop, which means that even though
the TR-BDF2 timestep is twice the length of that used by the SIQN scheme, the
TR-BDF2 scheme would be at a disadvantage: it would call the transport
operator 6 times per step, compared with the 4 times it is called by SIQN.

This can be reduced to 4 calls per step by introducing intermediate
variables:

.. math::

   \boldsymbol{X}^\ast=\boldsymbol{X}^n - \gamma\Delta t \mathcal{F}(\boldsymbol{X}^n), \qquad
   \boldsymbol{X}^p = \mathcal{T}^{2\gamma\Delta t}_{\overline{\boldsymbol{u}}_n^m}\left[\boldsymbol{X}^\ast\right], \qquad
   \boldsymbol{X}^q = \boldsymbol{X}^n + \boldsymbol{X}^p - \boldsymbol{X}^\ast.

The TR step gives an approximation to the transport of :math:`\boldsymbol{X}^n`
over :math:`2\gamma\Delta t`:

.. math::

   \boldsymbol{X}^q
   = \mathcal{T}^{2\gamma\Delta t}_{\overline{\boldsymbol{u}}_n^m}\left[\boldsymbol{X}^n - \gamma\Delta t \mathcal{F}(\boldsymbol{X}^n)\right] + \gamma\Delta t \mathcal{F}(\boldsymbol{X}^n)
   \approx \mathcal{T}^{2\gamma\Delta t}_{\overline{\boldsymbol{u}}_n^m}\left[\boldsymbol{X}^n\right].

In the BDF2 stage, rather than transporting :math:`\boldsymbol{X}^n` over
:math:`\Delta t` from time level :math:`n` to :math:`n+1`, the approximation
from the TR stage can be used and transported over the remaining
:math:`(1-2\gamma)\Delta t`. Since the transporting velocity and length of
timestep is then the same as that used for :math:`\boldsymbol{X}^m`, these
fields can be added together and transported once:

.. math::

   \mathcal{T}^{\Delta t}_{\overline{\boldsymbol{u}}_n^{n+1}}\left[\boldsymbol{X}^n\right]
   \approx \mathcal{T}^{(1-2\gamma)\Delta t}_{\overline{\boldsymbol{u}}_m^{n+1}}
   \left[\mathcal{T}^{2\gamma\Delta t}_{\overline{\boldsymbol{u}}_n^{m}}\left[\boldsymbol{X}^n \right]\right],

and

.. math::

   (1-\gamma_3)\mathcal{T}^{(1-2\gamma)\Delta t}_{\overline{\boldsymbol{u}}_m^{n+1}}\left[\boldsymbol{X}^q\right]
   + \gamma_3 \mathcal{T}^{(1-2\gamma)\Delta t}_{\overline{\boldsymbol{u}}_m^{n+1}}\left[\boldsymbol{X}^m\right]
   \approx
   \mathcal{T}^{(1-2\gamma)\Delta t}_{\overline{\boldsymbol{u}}_m^{n+1}}\left[(1-\gamma_3)\boldsymbol{X}^q+ \gamma_3 \boldsymbol{X}^m\right].

The whole scheme is then written as:

.. math:: :label: cheaper_trbdf2

   \boldsymbol{X}^{m} + \gamma\Delta t \mathcal{F}(\boldsymbol{X}^m) = \boldsymbol{X}^p, \qquad
   \boldsymbol{X}^{n+1} + \gamma_2\Delta t \mathcal{F}(\boldsymbol{X}^{n+1})
   = \mathcal{T}^{(1-2\gamma)\Delta t}_{\overline{\boldsymbol{u}}_m^{n+1}}\left[(1-\gamma_3)\boldsymbol{X}^q+ \gamma_3 \boldsymbol{X}^m\right].

This is the form of the TR-BDF2 scheme implemented in GungHo. Pseudo-code for
this "cheaper" version of the whole scheme is given in the
:ref:`Algorithm <science_guide_timestepping_tr_bdf2_algorithm>` section
below.

Transport
~~~~~~~~~

* The transport operator for dynamical prognostics in the TR stage is
  :math:`\mathcal{T}_{\overline{u}_{n}^{m}}^{2\gamma\Delta t}`, where
  :math:`\overline{u}_{n}^{m}=\tfrac{1}{2}(u^n+u^m)`;
* The transport operator for dynamical prognostics in the BDF2 stage is
  :math:`\mathcal{T}_{\overline{u}_{m}^{n+1}}^{(1-2\gamma)\Delta t}`, where
  :math:`\overline{u}_m^{n+1}=u^{n+1}`;
* The transport operator for tracers in the BDF2 stage is
  :math:`\mathcal{T}_{\overline{u}_{n}^{n+1}}^{\Delta t}`, where
  :math:`\overline{u}_{n}^{n+1}=\tfrac{1}{2}(u^n+u^{n+1})`.

Linear solver
~~~~~~~~~~~~~

As for the Semi-Implicit scheme, each stage's implicit problem is solved
using the linear problem of the form

.. math::

   \boldsymbol{X}' + \sigma \mathcal{I}[\boldsymbol{X}'] + \tau_i \Delta t \mathcal{L}[\boldsymbol{X}'] = \boldsymbol{R},

with relaxation parameters :math:`\sigma` and :math:`\tau_i`. For the TR
stage, these are

.. math::

   \tau_u = \gamma, \qquad \tau_\rho=\tau_\theta=2\gamma, \qquad \sigma=2\gamma,

while for the BDF2 stage they are

.. math::

   \tau_u = \gamma_2, \qquad \tau_\rho=\tau_\theta=1, \qquad \sigma=1.

Physics coupling
~~~~~~~~~~~~~~~~~

Physics parametrisation schemes are included using the same slow physics
:math:`\mathcal{P}_S`, fast physics :math:`\mathcal{P}_F`, and end-of-step
physics :math:`\mathcal{P}_E` operators defined for the
:ref:`Semi-Implicit <science_guide_timestepping_semi_implicit>` scheme.

1. **Slow Physics**:  This is evaluated on the
   start-of-timestep state :math:`\boldsymbol{X}^n` before the TR step.
2. **Fast Physics**: This is evaluated within the outer loop of the BDF2 step.
3. **End-of-step Physics**: denoted :math:`\mathcal{P}_E`. Evaluated after the
   Quasi-Newton iteration has completed for the BDF2 step.

All three physics operators can be included simultaneously, giving intermediate
states :math:`\boldsymbol{X}^\dagger`, :math:`\boldsymbol{X}^{q+\dagger}` and
:math:`\boldsymbol{X}^{n+\dagger}` in addition to the
:math:`\boldsymbol{X}^\ast`, :math:`\boldsymbol{X}^p` and
:math:`\boldsymbol{X}^q` defined above,

.. math:: :label: trbdf2_physics_coupling

   \begin{aligned}
   \boldsymbol{X}^\dagger &= \boldsymbol{X}^n - \Delta t\,\mathcal{P}_S(\boldsymbol{X}^n), \\
   \boldsymbol{X}^\ast &= \boldsymbol{X}^\dagger - \gamma\Delta t\,\mathcal{F}(\boldsymbol{X}^n), \\
   \boldsymbol{X}^p &= \mathcal{T}^{2\gamma\Delta t}_{\overline{\boldsymbol{u}}_n^m}\left[\boldsymbol{X}^\ast\right], \\
   \boldsymbol{X}^{m} + \gamma\Delta t\,\mathcal{F}(\boldsymbol{X}^m) &= \boldsymbol{X}^p, \\
   \boldsymbol{X}^q &= \boldsymbol{X}^\dagger + \boldsymbol{X}^p-\boldsymbol{X}^\ast, \\
   \boldsymbol{X}^{q+\dagger} &= \mathcal{T}^{(1-2\gamma)\Delta t}_{\overline{\boldsymbol{u}}_m^{n+1}}
   \left[(1-\gamma_3)\boldsymbol{X}^q+ \gamma_3 \boldsymbol{X}^m\right], \\
   \boldsymbol{X}^{n+\dagger} + \gamma_2\Delta t\,\mathcal{F}(\boldsymbol{X}^{n+\dagger})
   + \Delta t\,\mathcal{P}_F(\boldsymbol{X}^{n+\dagger})
   &= \boldsymbol{X}^{q+\dagger}, \\
   \boldsymbol{X}^{n+1} &= \boldsymbol{X}^{n+\dagger} - \Delta t\,\mathcal{P}_E(\boldsymbol{X}^{n+\dagger}).
   \end{aligned}

Moisture
~~~~~~~~~

To ensure correlation between water vapour and potential temperature,
moisture needs careful handling through the TR-BDF2 timestep. Post-solver
conservation corrections to moisture are included, which can be re-used to
update the moisture to be transported in the BDF2 stage.

It is easiest to understand the evolution of moisture through the dynamical
core alongside the evolution of dry density,

.. math:: :label: trbdf2_moisture_evolution

   \begin{array}{ll}
   \rho^p = \mathcal{T}^{2\gamma\Delta t}_{\overline{u}^m_n}[\rho^n],
   &
   m^p = \mathcal{T}^{2\gamma \Delta t}_{\overline{u}^m_n}[m^n;\rho^n], \\[2mm]
   \rho^m=\rho^p - \boldsymbol{\nabla\cdot}\boldsymbol{S}^d_{TR},
   &
   m^m = \dfrac{1}{\rho^m}\left[\rho^p m^p - \boldsymbol{\nabla\cdot}\boldsymbol{S}^m_{TR}\right], \\[2mm]
   \rho^{BDF}=(1-\gamma_3)\rho^p + \gamma_3\rho^m \equiv \rho^p -\gamma_3\boldsymbol{\nabla\cdot}\boldsymbol{F}^d_{TR},
   &
   m^{BDF} = \dfrac{1}{\rho^{BDF}}\left[\rho^p m^p - \gamma_3\boldsymbol{\nabla\cdot}\boldsymbol{F}^m_{TR} \right], \\[2mm]
   \rho^{q+\dagger} = \mathcal{T}^{(1-2\gamma)\Delta t}_{u^{n+1}}[\rho^n],
   &
   m^{q+\dagger} = \mathcal{T}^{(1-2\gamma)\Delta t}_{u^{n+1}}[m^{BDF};\rho^{BDF}], \\[2mm]
   \rho^{n+1}=\rho^{q+\dagger} - \boldsymbol{\nabla\cdot}\boldsymbol{F}^d_{BDF},
   &
   m^{n+1} = \dfrac{1}{\rho^{n+1}}\left[\rho^{q+\dagger} m^{q+\dagger} - \boldsymbol{\nabla\cdot}\boldsymbol{F}^m_{BDF}\right],
   \end{array}

where :math:`\boldsymbol{F}^d_{TR}` is the dry density flux from the TR
solver, and :math:`\boldsymbol{F}^d_{BDF}` is the corresponding flux from the
BDF2 stage. The superscript :math:`m` values indicate the values at the end of
the TR stage, and the fluxes :math:`\boldsymbol{S}^d_{TR}` and
:math:`\boldsymbol{S}^m_{TR}` are the fluxes corresponding to the increment
from the linear solver on the final outer iteration of the TR step.
If :math:`m^{BDF}` is calculated as described above, then the moisture mixing
ratios will be consistent with the dry density throughout the timestep.
The moist fluxes :math:`\boldsymbol{F}^m_{TR}` and
:math:`\boldsymbol{F}^m_{BDF}` are evaluated from the corresponding dry
fluxes, using the upwind values of :math:`m^p` and :math:`m^{q+\dagger}`
respectively. Note that :math:`m^m` is never used, so does not need to be
evaluated.

.. _science_guide_timestepping_tr_bdf2_algorithm:

Algorithm
---------

Equation :eq:`cheaper_trbdf2` above summarises the whole TR-BDF2
scheme, as implemented in GungHo. We now describe this in more detail as a
sequence of equations, using the same intermediate states
:math:`\boldsymbol{X}^\ast`, :math:`\boldsymbol{X}^p` and
:math:`\boldsymbol{X}^q` introduced above.

Both of the TR and BDF2 stages are solved using the outer/inner Quasi-Newton
structure described in the :ref:`Semi-Implicit <science_guide_timestepping_semi_implicit>` section.
The outer loop includes transport and fast physics terms, while the inner loop
calls the implicit forcing terms and the linear solver.

The explicit forcing contribution is evaluated once, using forces at time
level :math:`n`, before the TR stage begins,

.. math:: :label: trbdf2_alg_explicit_forcing

   \boldsymbol{X}^\ast = \boldsymbol{X}^n - \gamma\Delta t\,\mathcal{F}(\boldsymbol{X}^n).

This is transported over :math:`2\gamma\Delta t`, using the velocity
:math:`\overline{\boldsymbol{u}}_n^m`, to give the TR-stage transported
state,

.. math:: :label: trbdf2_alg_tr_transport

   \boldsymbol{X}^p = \mathcal{T}^{2\gamma\Delta t}_{\overline{\boldsymbol{u}}_n^m}\left[\boldsymbol{X}^\ast\right].

The TR-stage implicit problem is then solved for :math:`\boldsymbol{X}^m`,
using the outer/inner Quasi-Newton iteration with relaxation parameters
:math:`\tau_u=\gamma`, :math:`\tau_\rho=\tau_\theta=2\gamma`,
:math:`\sigma=2\gamma`,

.. math:: :label: trbdf2_alg_tr_solve

   \boldsymbol{X}^{m} + \gamma\Delta t\,\mathcal{F}(\boldsymbol{X}^m) = \boldsymbol{X}^p.

An approximation to the transport of :math:`\boldsymbol{X}^n` over
:math:`2\gamma\Delta t` is then formed,

.. math:: :label: trbdf2_alg_q

   \boldsymbol{X}^q = \boldsymbol{X}^n + \boldsymbol{X}^p-\boldsymbol{X}^\ast.

Convergence of the BDF2 stage can be accelerated by extrapolating from the
TR solution and the previous timestep, to give the first guess for the BDF2
stage, via

.. math:: :label: trbdf2_alg_bdf_first_guess

   \boldsymbol{X}^{n+1}_{(0)} = (1-1/(2\gamma))\boldsymbol{X}^n + 1/(2\gamma)\boldsymbol{X}^m.

Finally, the combination :math:`(1-\gamma_3)\boldsymbol{X}^q+ \gamma_3
\boldsymbol{X}^m` is transported over the remaining
:math:`(1-2\gamma)\Delta t`, using the velocity
:math:`\overline{\boldsymbol{u}}_m^{n+1}`, and the BDF2-stage implicit
problem is solved for :math:`\boldsymbol{X}^{n+1}`, using the outer/inner
Quasi-Newton iteration with relaxation parameters :math:`\tau_u=\gamma_2`,
:math:`\tau_\rho=\tau_\theta=1`, :math:`\sigma=1`,

.. math:: :label: trbdf2_alg_bdf_solve

   \boldsymbol{X}^{n+1} + \gamma_2\Delta t\,\mathcal{F}(\boldsymbol{X}^{n+1})
   = \mathcal{T}^{(1-2\gamma)\Delta t}_{\overline{\boldsymbol{u}}_m^{n+1}}
   \left[(1-\gamma_3)\boldsymbol{X}^q+ \gamma_3 \boldsymbol{X}^m\right].

This process is summarised by the pseudo-code below:

.. code-block:: text

   Slow physics:
      X^dag = X^n - dt * P_S(X^n)

   # TR Stage ------------------------------------------------------------------
   Explicit forcing:
      X_FE = X^dag - gamma * dt * F(X^n)
   Calculate transport predictor (X_FE)

   Set the first guess:
      X^m_(0) = X^dag

   # Outer Loop ----------------------------------------------------------------
   for j = 1 to n_o_TR:
       Update transporting velocity:
           u_bar_n_m^(j) = 0.5 * (u_n + u_m^(j-1))
       Transport:
           X^p_(j) = T_{u_bar_n_m^(j)}(X_FE)
       Implicit forcing:
           X_FI^(j,1) = X^{n+1}_(j-1) + alpha * dt * F(X^{n+1}_(j-1))

       X^m_(j,0) = X^m_(j-1)

       # Inner Loop ------------------------------------------------------------
       for k = 1 to n_i_TR^j:
           if k > 1:
              Implicit forcing:
                  X^m_FI^(j,k) = X^m_(j,k-1) + gamma * dt * F(X^m_(j,k-1))
           Form the residual:
               R_(j,k) = -( X^m_FI^(j,k) - X^p_(j) )
           Solve the linear system for the increment:
               S_TR[X'] = R_(j,k)
           Increment the state:
               X^m_(j,k) = X^m_(j,k-1) + X'

       Set final estimate:
           X^m_(j) = X^m_(j,n_i_TR^j)

   Advance the TR stage:
      X^m = X^m_(n_o_TR)
      X^p = X^p_(n_o_TR)

   # BDF2 Stage ----------------------------------------------------------------
   Form the approximate transport of X^n:
      X^q = X^dag + X^p - X_FE
   Set BDF2 pre-transport state:
      X_BDF = (1 - gamma3) * X^q + gamma3 * X^m

   Set the first guess:
      X^{n+1}_(0) = (1-1/(2*gamma)) * X^n + 1/(2*gamma) * X^m

   # Outer Loop ----------------------------------------------------------------
   for j = 1 to n_o_BDF:
       Update transporting velocity:
           u_bar_m_np1^(j) = u_{n+1}^(j-1)
       Transport:
           X^q_dag_(j) = T_{u_bar_m_np1^(j)}(X_BDF)
       Implicit forcing:
           X^{n+1}_FI^(j,1) = X^{n+1}_(j-1) + gamma2 * dt * F(X^{n+1}_(j-1))
       Fast physics predictor:
           X^star_(j) = X^q_dag_(j) + (X^{n+1}_FI^(j,1) - X^{n+1}_(j-1))
       Fast physics increment:
           dX_F^(j) = dt * P_F(X^star_(j))

       X^{n+1}_(j,0) = X^{n+1}_(j-1)

       # Inner Loop ------------------------------------------------------------
       for k = 1 to n_i_BDF^j:
           if k > 1:
               Implicit forcing:
                   X^{n+1}_FI^(j,k) = \
                      X^{n+1}_(j,k) + gamma2 * dt * F(X^{n+1}_(j,k-1))
           Form the residual:
               R_(j,k) = -( X^{n+1}_FI^(j,k) - X^q_dag_(j) - dX_F^(j) )
           Solve the linear system for the increment:
               S_BDF[X'] = R_(j,k)
           Increment the state:
               X^{n+1}_(j,k) = X^{n+1}_(j,k-1) + X'

       Set final estimate:
           X^{n+1}_(j) = X^{n+1}_(j,n_i_BDF^j)

   End-of-step physics:
      X^{n+1} = X^{n+1}_(n_o_BDF) - dt * P_E(X^{n+1}_(n_o_BDF))
   Advance the timestep:
      X^n = X^{n+1}

Tracer Transport
-----------------

Tracers that do not feed back onto the dynamical core can be transported in
just the final outer iteration of the BDF2 stage. If following the advective
form of the transport equation, tracers can be transported using the
time-centered wind
:math:`\overline{\boldsymbol{u}}_n^{n+1} = \tfrac{1}{2}(\boldsymbol{u}^n + \boldsymbol{u}^{n+1})`.
Tracers that are transported conservatively can use the mass flux used to
transport the dry density over the combination of the TR and BDF2 steps.

References
----------

.. [Tumolo2015] Tumolo, G. and Bonaventura, L. (2015). A semi-implicit,
   semi-Lagrangian discontinuous Galerkin framework for adaptive numerical
   weather prediction.
   *Quarterly Journal of the Royal Meteorological Society*, 141(692), 2582-2601.
   https://doi.org/10.1002/qj.2544
