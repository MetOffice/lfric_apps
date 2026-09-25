.. _tangent_linear_mol_transport:

MoL  Transport
====================

The Method of Lines (MoL) transport scheme is a Finite Volume Scheme that uses Runge-Kutta timestepping. For the linear model, we use the notation that :math:`\bar{u}` is the linearisation state and :math:`u'` is the perturbation.

Nonlinear MoL
----------

The nonlinear continuous equation is

.. math:: \frac{\partial \theta}{\partial t} = - u. \nabla \theta

which is discretised using an upwind scheme. On each face, 

.. math::  \begin{aligned} \frac{\partial \theta}{\partial t}
	                  & = - u. \nabla_L \theta & \text{   if } u>0 \\ 
                          & = - u. \nabla_R \theta & \text{   if } u<0 \\
			  & =   0                  & \text{   if } u=0 \\
	  \end{aligned}

where :math:`\nabla_L` uses the upwind reconstruction calculated from the left and :math:`\nabla_R` uses the upwind reconstruction calculated from the right.
As required, :math:`\frac{\partial \theta}{\partial t}` tends to 0 as u tends to 0.
			    
This is implemented by initialising
:math:`\frac{\partial \theta}{\partial t}` as :math:`0` and then adding
on an increment if :math:`u>0` or :math:`u<0`.

Linearised MoL
----------

The linearised continuous equation is

.. math:: \frac{\partial \theta'}{\partial t} = - \bar{u}. \nabla \theta' - u'. \nabla \bar{\theta}

This is discretised as

.. math::  \begin{aligned} \frac{\partial \theta'}{\partial t}
	                  & = - \bar{u}. \nabla_L \theta' - u'. \frac{1}{2} (\nabla_L + \nabla_R) \bar{\theta}
			  & \text{   if } \bar{u}>0 \\ 
                          & = - \bar{u}. \nabla_R \theta' - u'. \frac{1}{2} (\nabla_L + \nabla_R) \bar{\theta}
			  & \text{   if } \bar{u}<0 \\
			  & = - u'. \frac{1}{2} (\nabla_L + \nabla_R) \bar{\theta}
			  & \text{   if } \bar{u}=0 \\
	  \end{aligned}

which, as required, tends to
:math:`- u'. \frac{1}{2} (\nabla_L + \nabla_R) \bar{\theta}` as
:math:`\bar{u}` tends to :math:`0`.

This is implemented by initialising
:math:`\frac{\partial \theta'}{\partial t}` as
:math:`- u'. \frac{1}{2} (\nabla_L + \nabla_R) \bar{\theta}` and then
adding an increment if :math:`\bar{u}>0` or :math:`\bar{u}<0`.

Runge-Kutta
-----------

The nonlinear Runge Kutta scheme is:

.. math::

   \begin{aligned}
     u^{(i)} & = u^n + \Delta t \sum_{k=0}^{i-1} \beta_{i,k} f(u^{(k)}) \quad i=1,\ldots,m \\
     u^{n+1} & = u^{(m)}
     \end{aligned}

The tangent linear Runge Kutta scheme comprises two steps. First the linearisation
state is computed:

.. math:: \bar{u}^{(i)} = \bar{u}^n + \Delta t \sum_{k=0}^{i-1} \beta_{i,k} f(\bar{u}^{(k)}) \quad i=1,\ldots,m

then the perturbations are computed:

.. math::

   \begin{aligned}
     u'^{(i)} & = u'^n + \Delta t \sum_{k=0}^{i-1} \beta_{i,k} f'(u'^{(k)},\bar{u}^{(k)}) \quad i=1,\ldots,m \\
     u'^{n+1} & = u'^{(m)}
     \end{aligned}


If ``transport_efficiency`` is set to true, then the ``ls_state`` is not updated, so that the same linearisation state is used for all iterations.
