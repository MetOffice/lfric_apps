.. _tangent_linear_coding:

General principles
=====================

To obtain the tangent linear model, we need to differentiate the
nonlinear model. There are two main ways to achieve this: differentiate the nonlinear continuous equations and then discretize, or differentiate the nonlinear discretized equations.

Mainly the approach of differentiating the discretized equations is used here - the so-called line-by-line approach, as each line of code can be treated individually. However, some parts such as the semi-implicit solve take the approach of differentiating the nonlinear equations and then discretizing (although in fact as the solver is already linear, in practice this just means reusing the nonlinear code).

The following examples demonstrate how the line-by-line approach is implemented, to differentiate the discretized equations.

First derivative
----------------

For the most simple example, if the nonlinear model is :math:`y = f(x)`,
then the tangent linear model is

.. math:: \delta y = \frac{df}{dx} \delta x

For example, if the nonlinear model code is

::

     y = x**2

then :math:`f(x)=x^2`, and the first derivative is :math:`f'(x) = 2x`.
So the linear model code is

::

     y = 2 * ls_x * x

where now ``y``, and ``x`` are the perturbations, and ``ls_x`` is the
linearisation state.

Note that we could rewrite this in other variable names as e.g.
``delta_y = 2 * x * delta_x``, but in the LFRic code we have chosen to
keep the variable names for the perturbations identical to the nonlinear
model variable names, and to introduce ``ls_`` to indicate the
linearisation state. This linearisation state is produced from a
separate nonlinear model run.

The chain rule
--------------

For a more complicated expression, then the chain rule of
differentiation is required

.. math::

   \begin{aligned} 
     y & = f(g(x)) \\
     \delta y & = \frac{df}{dx} \delta x \\
              & = \frac{df}{dg} \frac{dg}{dx} \delta x
     \end{aligned}

For example, if the nonlinear model code is

::

     y = (sin v)**2

so

.. math::

   \begin{aligned}
   f(g) & = g^2 \\
   g(x) & = sin(x)\\
   f'(g) & = 2 g = 2 sin(x) \\ 
   g'(x) & = cos(x)
   \end{aligned}

then the linear model code is

::

     y = 2 * sin(ls_v) * cos(ls_v) * v

Multivariable calculus
----------------------

When dealing with more than one variable, then the rules need to be
extended to multi-variable calculus.

If the nonlinear model is

.. math:: y = f(x,z)

the linear model is

.. math:: \delta y = \frac{\partial f}{\partial x} \delta x + \frac{\partial f}{\partial z} \delta z

and the multivariate chain rule gives

.. math::

   \begin{aligned} 
     y & = f(g_1(x_1), ..., g_n(x_n)) \\
     \delta y & = \sum_i \frac{\partial f}{\partial g_i} \frac{d g_i}{d x_i} \delta x_i
     \end{aligned}
