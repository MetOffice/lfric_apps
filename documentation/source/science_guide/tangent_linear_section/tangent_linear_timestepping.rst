.. _tangent_linear_timestepping:

Iterated (Semi-) Implicit Timestepping
====================

Nonlinear Equations
-------------------

The continous nonlinear equations are

.. math::

   \begin{aligned}
     \frac{\partial u}{\partial t} &= - u . \nabla u - 2 \Omega \times u - \nabla \Phi - c_p \theta \nabla \Pi \\
     \frac{\partial \rho}{\partial t} & = - \nabla . (\rho u) \\
     \frac{\partial \theta}{\partial t} & = - u. \nabla \theta
     \end{aligned}

together with the nonlinear equation of state:

.. math:: \Pi ^{\left( \frac{1-\kappa}{\kappa} \right)} = \frac{R}{p_0} \rho \theta

where the prognostic variables are :math:`u` the velocity vector, :math:`\rho` the density, :math:`\theta` the potential temperature, and :math:`\Pi` the Exner pressure.  The constants are :math:`\Omega` the rotation vector, :math:`\Phi` the geopotential and :math:`c_p` the specific heat capacity.

Discretisation
--------------------
In the iterated-implicit timestepping, the equations are solved iteratively
using a Newton method. If the full nonlinear discrete equations are
written as

.. math:: \mathcal{R} (\mathbf{x}^{n+1}) =0

then this can be solved iteratively as

.. math::

   \begin{aligned}
     \mathbf{x}^{n+1}_{(k+1)} & =  \mathbf{x}^{n+1}_{(k)} +  \mathbf{x}'_{(k)} \\
     \mathcal{L} (\mathbf{x}^*) \mathbf{x}'_{(k)} & = - \mathcal{R} ( \mathbf{x}^{n+1}_{(k)})
     \end{aligned}

where :math:`\mathcal{L}` is an approximation to the Jacobian and
:math:`\mathbf{x}^*` is the reference state.

The reference state :math:`u^*` is zero whilst the other reference state
variables are set to :math:`\mathbf{x}^n`. Relaxation parameters
:math:`\tau <1` are also introduced to give extra weight to the diagonal
terms (and hence similar to solving the resulting linear equation using
the Jacobi method or a damped Newton method).

This gives the approximation to the Jacobian as

.. math::

   \mathcal{L}(\mathbf{x}^*) \mathbf{x}' = \left\{
     \begin{aligned}
       u' & + \tau_u \Delta t c_P (\theta' \nabla \Pi^* + \theta^* \nabla \Pi ') \\
       \rho' & + \tau_\rho \Delta t \nabla.(\rho^* u') \\
       \theta' & + \tau_\theta \Delta t u' \nabla. \nabla \theta^* \\
   & \frac{1-\kappa}{\kappa}  \frac{\Pi'}{\Pi^*}  - \frac{1}{E} \frac{\rho'}{\rho^*} - \frac{1}{E} \frac{\theta'}{\theta^*}
     \end{aligned}
   \right.

where

.. math:: E = \frac{p_0}{R} \frac{\Pi^{*(1-\kappa)/\kappa}}{\theta^* \rho^*}

The procedure is:

- | Let :math:`\mathbf{x}^{n+1}_{(0)} = \mathbf{x}^n`
  | The first iteration is given by the state at the previous timestep.

- | Set
    :math:`\mathbf{x}* \equiv (\rho*, \theta*, \Pi* ) = (\rho^n, \theta^n,\Pi^n)`
  | The reference state is given by the state at the previous timestep
    (for density, potential temperature and exner pressure only).

- Loop k=1, K

  - | Compute :math:`\mathcal{R}(\mathbf{x}^n, \mathbf{x}^{n+1}_{(k)})`

  - Solve the linear equation for :math:`\mathbf{x}'_{(k+1)}` and
    hence for :math:`\mathbf{x}^{n+1}_{(k+1)}`.

Linear Equations
-----------------------

We can proceed in a similar way with the linear model. Using the notation of :math:`\delta \mathbf{x}` is the perturbation and :math:`\mathbf{x}` is the linearisation state,

.. math::

   \begin{aligned}
     \delta \mathbf{x}^{n+1}_{(k+1)} & =  \delta \mathbf{x}^{n+1}_{(k)} +  \delta \mathbf{x}'_{(k)} \\
     \mathcal{L} (\mathbf{x}^*) \delta \mathbf{x}'_{(k)} & = - R ( \mathbf{x}_{(k)}) \delta \mathbf{x}_{(k)}
     \end{aligned}

where :math:`R` is the linearized version of :math:`\mathcal{R}`. As the operator :math:`\mathcal{L}` is already linear, this does not need to be linearized.

The procedure is:

- | Let :math:`\delta \mathbf{x}^{n+1}_{(0)} = \delta \mathbf{x}^n`
  | The first iteration is given by the perturbation at the previous
    timestep.

- | Set
    :math:`\mathbf{x}* \equiv (\rho*,\theta*, \Pi* ) = (\rho^n, \theta^n, \Pi^n)`
  | The reference state is given by the linearisation state at the
    previous timestep - so that the linear operator is identical to the
    that in the nonlinear model.

- Loop k=1, K

  - Compute
    :math:`R ( \mathbf{x}_{(k)}) \delta \mathbf{x}_{(k)}`

  - Solve for :math:`\delta \mathbf{x}'_{(k)}` using :math:`\mathcal{L} (\mathbf{x}^*) \delta \mathbf{x}'_{(k)} = R ( \mathbf{x}_{(k)}) \delta \mathbf{x}_{(k)}` and hence for
    :math:`\delta \mathbf{x}^{n+1}_{(k+1)}`.

The linear operator :math:`\mathcal{L} ( \bar{\mathbf{x}}*)` is exactly
the same as for the nonlinear model. The solution of the equation (often known as the solver) uses the same method as the nonlinear model. This could be multi-grid or Krylov iterative methods, or a combination of the both. A very exact line-by=line procedure would linearise every iteration of the solver procedure, using the nonlinear solution at every step. However, this is not necessary and the more linearise then discretize approach is taken, to form the linear equation and then solve.

The right hand side term :math:`R ( \mathbf{x}_{(k)}) \delta \mathbf{x}_{(k)}` require both the perturbation :math:`\delta \mathbf{x}_{(k)}` and the linearisation state :math:`\mathbf{x}_{(k)}` at the same iteration. Therefore, for an exact tangent linear semi-implicit algorithm, the
nonlinear semi-implicit step needs to be run first, and the intermediate
``ls_state`` values need to be stored and used in the following tangent
linear semi-implicit step.

In practice the impact of the variation of the linearisation state through the timestep is very small. Furthermore, the evolution of the linearisation state in the linear model timestepping does not include any physics and will be at a different resolution to the outer loop nonlinear model, and hence will not be very accurate. Therefore, it is advised for both efficiency and accuracy reasons to use the ``fixed_ls`` option. This option applies the ``ls_state`` values set at the beginning of the timestep to all iterations, without recalculation.

