.. _tangent_linear_components:

Right-hand-side Components
=====================

Kinetic energy gradient
-----------------------

The nonlinear model is

.. math::

   \begin{aligned}
     f & = \nabla (\frac{1}{2} u^2) \\
       & = g(z(u)) \\
     \text{where} \quad z(u) & = \frac{1}{2} u^2 \\
                        g(z) & = \nabla(z)
     \end{aligned}

The linear model is,

.. math::

   \begin{aligned}
     \delta f & = \frac{d g}{d z} \frac{ dz}{du} \delta u\\
              & = \nabla (u \delta u)
     \end{aligned}

Note that an alternative way to derive this is to form the difference
between the perturbed and unperturbed, and remove the small terms

.. math::

   \begin{aligned}
    \delta f & = f(\bar{u}+ \delta u) - f(\bar{u}) \\
       & = \frac{1}{2} \nabla (\bar{u} + \delta u)^2 - \frac{1}{2} \nabla \bar{u}^2 \\
       & = \frac{1}{2} \nabla (\bar{u}^2 + 2\bar{u} \delta u + \delta u^2) - \frac{1}{2} \nabla \bar{u}^2 \\
       & \approx \nabla (\bar{u} \delta u)
     \end{aligned}

Potential temperature
---------------------

The nonlinear model is

.. math:: \theta_e  = \theta (m_g/m_t)

The linear model is,

.. math::

   \begin{aligned}
     \delta \theta_e & = \frac{\partial \theta_e}{\partial \theta} \delta \theta
                       + \frac{\partial \theta_e}{\partial m_g} \delta m_g
                       + \frac{\partial \theta_e}{\partial m_t} \delta m_t \\
                     & = \frac{m_g}{m_t} \delta \theta
                       + \frac{\theta}{m_t} \delta m_g
                       -  \frac{\theta m_g}{m_t^2} \delta m_t
      \end{aligned}

using

.. math:: \frac{d (1/m_t) }{dx} = \frac{d (m_t^{-1}) }{dx} = -1 (m_t^{-2}) = - \frac{1}{m_t^2}

This can be simplified, by substituting in the equation for the
nonlinear model :math:`\theta_e  = \theta (m_g/m_t)`

.. math:: \delta \theta_e =  \theta_e \left( \frac{ \delta \theta}{\theta} + \frac{\delta m_g}{m_g} - \frac{\delta m_t}{m_t} \right)

Logspace
--------

The nonlinear model is

.. math::

   \begin{aligned}
     p   & = \Pi_i |\rho_i|^{c_i} \\
         & = \Pi_i u_i(y_i(z_i(\rho_i))) \\
         & = \Pi_i u_i \\
    \text{where} \quad u_i & = y_i ^{c_i} \\
     y_i & = | \rho_i | \\
         & = (z_i) ^{1/2} \\
     z_i & = \rho_i^2
     \end{aligned}

The derivatives are

.. math::

   \begin{aligned}
     \frac{ \partial p}{\partial u_i} & =  \Pi_{j \neq i} u_j =  \Pi_{j \neq i} | \rho_j |^{c_j} \\
     \frac{ \partial u}{\partial y_i} & = c_i y_i ^{c_i -1} = c_i | \rho_i | ^{c_i -1} \\
     \frac{ \partial y}{\partial z_i} & = \frac{1}{2} z_i^{-1/2} = \frac{1}{2} \frac{1}{| \rho_i | } \\
     \frac{ \partial z}{\partial \rho_i} & = 2 \rho_i
     \end{aligned}

Therefore the linear model is,

.. math::

   \begin{aligned}
     \delta p & =  \sum_i \frac{ \partial p}{\partial u_i}  \frac{ \partial u_i}{\partial y_i}\frac{ \partial y_i}{\partial z_i} \frac{ \partial z_i}{\partial \rho_i} \delta \rho_i \\
     & = \sum_i \left( \left( \Pi_{j \neq i} | \rho_j |^{c_j} \right) c_i  | \rho_i | ^{c_i -1} \frac{\rho_i}{| \rho_i |} \delta \rho_i \right)
     \end{aligned}

This can be simplified, by using a rearrangement of the nonlinear model

.. math::

   \begin{aligned}
     p   & = \Pi_i |\rho_i|^{c_i} \\
         & = |\rho_i|^{c_i}  \Pi_{j \neq i} | \rho_j|^{c_j} \\
   \text{so } \quad \Pi_{j \neq i} | \rho_j|^{c_j} &  = \frac{p}{|\rho_i|^{c_i}} \\
     \end{aligned}

Substituting this into the linear model, gives

.. math::

   \begin{aligned}
     \delta p & = p \sum_i c_i \frac{\rho_i }{ | \rho_i |^2} \delta \rho_i \\
              & = p \sum_i \frac{c_i }{ \rho_i} \delta \rho_i \\
     \end{aligned}

Newton iteration of the Equation of State
-----------------------------------------

The equation of state is applied iteratively in the semi-implicit
timestepping. The nonlinear model rhs is

.. math::

   \begin{aligned}
     N & = 1/E -1
     \end{aligned}

where

.. math::

   \begin{aligned}
   E & = \frac{p_0}{R} \frac{\Pi^\gamma}{\theta \rho} \\ 
   \gamma & = \frac{1-\kappa}{\kappa}
   \end{aligned}

The derivatives are:

.. math::

   \begin{aligned}
     \frac{\partial N}{\partial E} & = - E^{-2} \\
     \frac{\partial E}{\partial \Pi} & = \frac{p_0}{R} \gamma \Pi^{\gamma-1}\theta^{-1} \rho^{-1} = E \gamma \Pi^{-1} \\
   \frac{\partial E}{\partial \rho} & = - \frac{p_0}{R} \Pi^\gamma \rho^{-2} \theta^{-1} = - E \rho^{-1} \\
   \frac{\partial E}{\partial \theta} & = - \frac{p_0}{R} \Pi^\gamma \rho^{-1} \theta^{-2} = - E \theta^{-1}
   \end{aligned}

The linear model is

.. math::

   \begin{aligned}
   \delta N & = \frac{\partial N}{\partial E} \left( \frac{ \partial E}{\partial \Pi} \delta \Pi + \frac{ \partial E}{\partial \rho} \delta \rho + \frac{ \partial E}{\partial \theta} \delta \theta \right) \\
    &  = - \frac{1}{E} \left( \frac{1-\kappa}{\kappa} \frac{\delta \Pi}{\Pi} - \frac{\delta \rho}{\rho} - \frac{\delta \theta}{\theta} \right)
   \end{aligned}
