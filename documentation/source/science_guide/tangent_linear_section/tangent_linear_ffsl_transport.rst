.. _tl_ffsl_transport:

FFSL Transport
================

| The non-linear flux-form semi-Lagrangian (FFSL) transport scheme solves the conservative form transport
  equation in the form

  .. math:: q^{n+1} = q^{n} - \Delta{t} F(u,q^n)

  where :math:`q` is the variable we are transporting, :math:`u` is the
  advecting wind, :math:`\Delta{t}` is the time step, and :math:`F` is
  the divergence of the mass flux

  .. math:: F(u,q) = \frac{1}{V} \nabla \cdot \tilde{f}(u,q)

  where the flux :math:`\tilde{f}`, is computed using the chosen FFSL
  scheme.

The non-linear scheme solves the advective form transport equation in
the form

.. math::

       q^{n+1}_{adv} = \frac{q^{n} - \Delta{t} F(u,q^n)}{\sigma - \Delta{t} F(u,\sigma)}

where :math:`\sigma=1` is the unity field, and therefore
:math:`F(u,\sigma)` is the wind divergence. The unity field plays the
role of a psuedo-density.

Tangent Linear Transport
-------------------------

For the linear model we split the variables into a linear state (LS)
denoted with a bar, and a perturbation denoted with a prime

.. math:: q = \bar{q} + q' \hspace{1cm} u = \bar{u}+u'

The conservative transport equation becomes

.. math:: \frac{\partial \bar{q}}{\partial t} + \frac{\partial q'}{\partial t} + \nabla \cdot( (\bar{u}+u')(\bar{q}+q') ) = 0

The product of perturbations is assumed to be zero, and so multiplying
out the divergence term gives

.. math:: \frac{\partial \bar{q}}{\partial t} + \frac{\partial q'}{\partial t} + \nabla \cdot( \bar{u}\bar{q} ) + \nabla \cdot( \bar{u}q' ) + \nabla \cdot( u'\bar{q} ) = 0

The linear state conservative transport equation is

.. math:: \frac{\partial \bar{q}}{\partial t} + \nabla \cdot( \bar{u}\bar{q} ) = 0

Leaving the conservative perturbation transport equation as

.. math::

       \frac{\partial q'}{\partial t} + \nabla \cdot( \bar{u}q' ) + \nabla \cdot( u'\bar{q} ) = 0

Flux Computation
----------------------

The flux computation can be split into two parts.

The first term, corresponding to :math:`\bar{u}q'`, can be computed
using standard FFSL operators, provided the scheme is linear (i.e. no
flux limiting). This is the reconstruction of the perturbation field
with the ls wind, and we denote this :math:`\tilde{f}_1`.

The second term, corresponding to :math:`u'\bar{q}`, uses the
perturbation departure distance from the ls departure point. This is
computed as, where the subscript :math:`d` signifies the ls wind
departure cell,

.. math:: \tilde{f}_2(u',\bar{q}) = u' \bar{q}_d

The constant reconstruction is to ensure that this term is linear in the
perturbation wind.

We can now use the notation

.. math::

   \begin{aligned}
       F_1(u,q) & = \frac{1}{V} \nabla \cdot \tilde{f}_1(u,q) \\
       F_2(u,q) & = \frac{1}{V} \nabla \cdot \tilde{f}_2(u,q)
   \end{aligned}

Advective Update
--------------------

If we expand the variables
into linearisation state and perturbation parts we find the flux
divergences can be written as

.. math:: F(u,q) = F(\bar{u},\bar{q}) + F_1(\bar{u},q') + F_2(u',\bar{q})

and, as :math:`\sigma = \bar{\sigma}`,

.. math:: F(u,\sigma) = F(\bar{u},\sigma) + F_2(u',\sigma)

and the advective form becomes

.. math::

       q^{n+1}_{adv} = \frac{\bar{q} + q' - \Delta{t} F(\bar{u},\bar{q}) - \Delta{t}F_1(\bar{u},q') - \Delta{t}F_2(u',\bar{q}) }{\sigma - \Delta{t} F(\bar{u},\sigma)- \Delta{t} F_2(u',\sigma)}

Let

.. math::

   \begin{aligned}
       \bar{N} &= \bar{q} - \Delta{t} F(\bar{u},\bar{q}) \\
       N' &= q' - \Delta{t}F_1(\bar{u},q') - \Delta{t}F_2(u',\bar{q}) \\
       \bar{D} &= \sigma - \Delta{t} F(\bar{u},\sigma) \\
       D' = - \Delta{t} F_2(u',\sigma)
   \end{aligned}

and then this can be written as

.. math:: q^{n+1}_{adv} = \frac{\bar{N}+N'}{\bar{D}+D'}

Noting that

.. math::

       \bar{q}^{n+1}_{adv} = \frac{\bar{N}}{\bar{D}}

then the advective perturbation field can be found from
:math:`q^{n+1}_{adv} - \bar{q}^{n+1}_{adv}`, i.e.

.. math:: q'^{n+1}_{adv} = \frac{\bar{N}+N'}{\bar{D}+D'} - \frac{\bar{N}}{\bar{D}}

Using the first terms of the expansion
:math:`(1+D'/\bar{D})^{-1} \approx 1 - D'/\bar{D}` we get

.. math::

   \begin{aligned}
       q'^{n+1}_{adv}=& \frac{\bar{N}+N'}{\bar{D}}(1-D'/\bar{D}) - \frac{\bar{N}}{\bar{D}} \\
                     =& \frac{N'}{\bar{D}} - \frac{\bar{N}D'}{\bar{D}\bar{D}} - \frac{N'D'}{\bar{D}\bar{D}}
   \end{aligned}

Using the fact products of prime
variables equals zero, this simplifies to

.. math:: q'^{n+1}_{adv} = \frac{N' - D' \bar{q}^{n+1}_{adv}}{\bar{D}}

Horizontal Splitting
-----------------------

The horizontal FFSL scheme uses COSMIC/Lin-Rood splitting which has, for
the non-linear scheme, the form

.. math::

   \begin{aligned}
       q^x &= q^{n} - \frac{\Delta}{2}F^x_{adv}(u,q^n) \\
       q^x &= q^{n} - \frac{\Delta}{2}F^y_{adv}(v,q^n) \\
       q^{n+1} &= q^{n} - \Delta{t} F^x(u,q^y)- \Delta{t} F^y(v,q^x)
   \end{aligned}

Transport Efficiency
------------------------

To improve efficiency, we can make the approximation that
:math:`\bar{q}^{n+1}_{adv} \approx \bar{q}^{n}`. This removes the need
to recompute :math:`\bar{q}^{n+1}_{adv}` at every step.
