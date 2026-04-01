====================
The PC2 Cloud Scheme
====================

:Author: D. Wilson, A. Bushell, C. Morcrette, V. Varma\ :math:`^{1}`,
         M. Whitall

.. role:: raw-latex(raw)
   :format: latex
..

PC2 developers
==============

We would like to acknowledge those who developed the PC2 cloud scheme:
Damian Wilson, Andrew Bushell, David Gregory, Amanda Kerr-Munslow, John
Edwards, Jeremy Price, Cyril Morcrette, Martin Sharpe, Thomas Mirfield,
Ian Boutle. Many others offered considerable help, advice and analysis,
including Roy Kershaw, Malcolm Brooks, Richard Forbes and Alejandro
Bodas-Salcedo, and we would like to thank them all for their input.

Introduction
============

This document describes the PC2 *(prognostic cloud, prognostic
condensate)* cloud scheme. It should be seen as a complete reference
source for the scheme’s physical assumptions, numerical techniques,
application to the Unified Model and coding within the Unified Model. It
does not describe results from the scheme, please refer to the various
reports and papers written on this. Except where commented on
explicitly, the description applies to the PC2:66 version of the PC2
scheme, which is the version that will be available at UM6.5. The
version available at 6.4 is PC2:64.

This paper will first introduce the concepts that underlie cloud
schemes, before developing a study of the theoretical behaviour of the
prognostic PC2 scheme under certain, well-defined, situations. The next
sections shows how the theory can be applied to the physical and
dynamical processes represented in the Unified Model. Finally, we
outline the way in which the PC2 scheme is implemented within the code
of the Unified Model.

Cloud schemes
-------------

The basic requirements of any cloud scheme within a large-scale model
are to:

- calculate the amount of condensation (from water vapour to liquid
  water or vice-versa) within each gridbox each timestep

- to calculate or update the cloud fractions for use by the radiation
  and large-scale precipitation schemes (or any other physics scheme).

Depending on the model involved, cloud schemes may also treat the
deposition / sublimation process from vapour to ice. The problem is
straightforward to solve if one is allowed to assume that there is no
variability of moisture or temperature on a scale of a model gridbox. In
this case the cloud fraction scheme is redundant and only the
condensation part remains, which may be solved diagnostically using the
instantaneous condensation assumption in section `2.2 <#sec:s_dist>`__.
However, the ‘no-variability’ assumption is poor until very high
resolutions close to, or maybe exceeding, 1 km in the horizontal are
reached. Although we may eventually assume that computer power will
enable such resolutions to be reached globally, for many years we will
need a subgrid-scale cloud scheme to properly account for the
variability in the atmosphere. This is the principal challenge of cloud
parametrization.

There are several approaches to take to the solution of the problem,
although they are not as independent as often portrayed, since they
nearly all require the same instantaneous condensation assumption
(discussed in section `2.2 <#sec:s_dist>`__). Hence there are
mathematical links between all the approaches. *The following are all
valid structures to use in this respect.*

- One may diagnose cloud fractions and condensate contents from
  knowledge of gridbox mean variables. This forms the basis of the
  :raw-latex:`\cite{smith90}` scheme, which is described in .

- A mixed scheme, such as :raw-latex:`\cite{sundqvist1978}` uses a
  prediction of condensate contents, but a diagnostic cloud fraction.

- Alternatively, one may predict cloud fraction and condensate content
  changes as a result of each modelled process. This forms the basis of
  the :raw-latex:`\cite{t93}` scheme and the PC2 scheme.

- Hybrid schemes, such as :raw-latex:`\cite{t02}`, will predict various
  moments of the subgrid-scale variability, and use this knowledge to
  diagnose the cloud fraction and condensate contents.

Many years of experience of the results from the
:raw-latex:`\cite{smith90}` scheme have highlighted deficiences in the
diagnosis of cloud from this scheme, which we feel can only be tackled
by adding the memory of cloud history available by using a prognostic
based scheme. We chose to develop a scheme that directly specified the
impacts on observable prognostics (condensates and cloud fractions, as
in :raw-latex:`\cite{t93}`) rather than on moments of a probability
density function (as in :raw-latex:`\cite{t02}`). This is because we
believe it is easier to physically relate (and hence parametrize) the
effect processes to quantities such as cloud fraction and condensate
rather than to the more abstract quantities of moments of a probability
density function of moisture. However, although the PC2 scheme is
similar to :raw-latex:`\cite{t93}` in its very basic prognostic variable
structure, the assumptions behind the formulation of the prognostic
terms in PC2 are very different and much improved. The PC2 scheme should
not be considered to be merely an extension of :raw-latex:`\cite{t93}`.

In particular, we wish to use a prognostic formulation in order to link
the detraiment of moisture from convection directly to cloud fraction,
and to break the hard diagnostic link between cloud fraction and
condensate. These major features of the :raw-latex:`\cite{t93}` scheme
provide the motivation to develop the PC2 cloud scheme.

.. _`sec:s_dist`:

The ‘s’ distribution
--------------------

Most cloud schemes are based on the concept of a distribution of
fluctuations of moisture and temperature in the gridbox. Here we
mathematically formalize this concept, since it is used both in the PC2
scheme and the :raw-latex:`\cite{smith90}` scheme.

This method was first formulated by :raw-latex:`\cite{m77}` and
:raw-latex:`\cite{sommeria_deardorff_1977}` for large-eddy simulations.
It can also be applied to larger scale models. It allows us to calculate
vapour and liquid contents and liquid cloud fraction from knowledge only
of the combined vapour+liquid content, :math:`\overline{q_T}`, and the
liquid temperature, :math:`\overline{T_L}`. These variables are
unchanged during condensation processes, so it is useful to write the
cloud scheme in terms of these variables.

In this derivation we will consider only liquid condensate. We assume,
as above, that, locally, the water content in a cloud is such as to
remove any supersaturation. This gives the equation

.. math::

   q_{cl} = q_T - q_{sat}(T,p)
   \label{eq:basic_qcl}

assuming that :math:`q_T > q_{sat}(T,p)` (:math:`q_{cl}` will be zero
otherwise). :math:`q_T` is the local total water content, equal to the
sum of the condensate :math:`(q_{cl})` plus the vapour :math:`(q)`,
:math:`T` is the temperature, :math:`p` is the pressure and
:math:`q_{sat}(T,p)` is the saturation specific humidity at temperature
T and pressure p *with respect to liquid water*. (Many earlier
diagnostic cloud schemes use a similar instantaneous condensation
assumption for ice, which would mean that :math:`q_{sat}` must be taken
with respect to ice when :math:`T < 0 ^{\circ} C`, but the Unified Model
does not). We now introduce the liquid temperature (:math:`T_L`), where
:math:`T_L` is given by

.. math::

   T_L = T - \frac{L}{c_p} q_{cl} , 
   \label{eq:tl}

and :math:`L` is the latent heat of vaporization and :math:`c_p` is the
heat capacity of air. Note that :math:`T_L` is unaffected by changes of
phase between vapour and liquid. We now write
(`[eq:basic_qcl] <#eq:basic_qcl>`__) as an *equality*

.. math::

   q_{cl} = q_T - \left( q_{sat}(T_L,p) + \alpha (T - T_L) \right)
   \label{eq:alpha_t_tl}

where

.. math::

   \alpha = \frac{ q_{sat}(T,p) - q_{sat}(T_L,p) }{T - T_L } .
   \label{eq:alpha}

Using (`[eq:tl] <#eq:tl>`__) in (`[eq:alpha_t_tl] <#eq:alpha_t_tl>`__)
gives the expression

.. math:: q_{cl} = q_T - q_{sat} (T_L) - \alpha \frac{L}{c_p} q_{cl}

or

.. math::

   q_{cl} = a_L \left( q_T - q_{sat}(T_L,p) \right)
   \label{eq:l_eq_al}

where :math:`a_L` is given by

.. math::

   a_L = \left( 1 + \alpha \frac{L}{c_p} \right) ^{-1} .
   \label{eq:a_L}

Thus (`[eq:basic_qcl] <#eq:basic_qcl>`__) has been rewritten *exactly*
in terms of the conserved variables, :math:`q_T` and :math:`T_L`,
although the temperature, :math:`T`, does remain in the definition of
:math:`a_L`. We will need to consider variations across a gridbox for a
parametrization scheme, so we expand the expression for condensate
(`[eq:l_eq_al] <#eq:l_eq_al>`__) into terms relating to the gridbox mean
and variation from the gridbox mean.

.. math::

   q_{cl} = \overline{ a_L \left( q_T - q_{sat}(T_L,p) \right)}
   + [ a_L \left( q_T - q_{sat}(T_L,p) \right) ]'
   \label{eq:bar_plus_pri1}

where :math:`\overline{\phi}` represents the mean of a distibution of
:math:`\phi` and :math:`\phi = \overline{\phi} + {\phi}'`. The
expression (`[eq:bar_plus_pri1] <#eq:bar_plus_pri1>`__) is *exact* when
using the definition of :math:`\alpha` given in
(`[eq:alpha] <#eq:alpha>`__).

The idea of a PDF scheme is to calculate the first (mean) term,
:math:`\overline{\phi}`, from the known gridbox mean parameters,
:math:`q_T`, :math:`T_L` and :math:`p`, and to parametrize the
distribution of the second, variable term, :math:`{\phi}'`.
Unfortunately, the mean term is difficult to write in terms of the
gridbox mean variables :math:`\overline{q_T}` and :math:`\overline{T_L}`
because :math:`q_{sat}(T_L,p)` is not a linear function of :math:`T_L`
(or of :math:`p`). In order to proceed, we will now make an
*approximation* that :math:`q_{sat}(T,p)` is a linear function of
:math:`T_L` and :math:`p`. This equivalently implies that :math:`a_L`
and :math:`\alpha` are approximated as being constant across the
gridbox. The expression now becomes more tractable,
(`[eq:bar_plus_pri1] <#eq:bar_plus_pri1>`__) becoming:

.. math::

   q_{cl} = a_L \left( \overline{q_T} - q_{sat}(\overline{T_L},\overline{p})
   \right)  + a_L \left(  {q_T}' - \alpha {T_L}' - \beta {p}' \right)
   \label{eq:l_eq_bar_plus_pri}

where :math:`\beta = {\frac{\partial q_{sat}}{\partial p}}` at constant
temperature. The first term is connected with the mean properties of the
gridbox, and is written as :math:`Q_c`, the second term is connected
with the deviation of the local conditions from the mean and is written
as :math:`s`.

.. math::

   Q_c = a_L \left( \overline{q_T} - q_{sat}(\overline{T_L},\overline{p}) \right)
   \label{eq:qc_eq_qt-qs}

.. math::

   s = a_L \left(  {q_T}' - \alpha {T_L}' - \beta {p}' \right)
   \label{eq:s}

This gives the equation

.. math::

   q_{cl} = Q_c + s
   \label{eq:l_qc_s}

with the assumption that :math:`s \ge -Q_c` (i.e. :math:`q_{cl} \ge 0`).
If :math:`s < -Q_c` then :math:`q_{cl} = 0`. The term :math:`a_L` can be
calculated using (`[eq:a_L] <#eq:a_L>`__) from
(`[eq:alpha] <#eq:alpha>`__) with gridbox mean temperatures, i.e.

.. math::

   \alpha = \frac{ q_{sat}(\overline{T},\overline{p})
   - q_{sat}(\overline{T_L},\overline{p}) }{\overline{T} -
   \overline{T_L} } .
   \label{eq:alpha_mean}

This definition of :math:`\alpha` and :math:`a_L` will retrieve an
*exact* value for the gridbox mean :math:`\overline{q_{cl}}` *if* the
distribution is monodispersed. Hence it is the sensible form to use for
a purely diagnostic representation such as :raw-latex:`\cite{smith90}`
where we explicitly consider distributions of :math:`s`. Strictly, the
linear approximation implies that other approximations for
:math:`\alpha` are valid: PC2 will do this (see section
`3.2.3 <#sec:homog_num_app>`__) since we are concerned in PC2 with the
best estimate of the *changes* to :math:`\overline{q_{cl}}`, not the
best estimate of :math:`\overline{q_{cl}}` itself.

We now assume that within any particular gridbox a distribution
:math:`G` of :math:`s` occurs (with mean, by definition, of zero).
Considering cloud to be where the water content is greater than zero
(i.e. where :math:`s > -Q_c`) gives an expression for the liquid cloud
*volume* fraction, :math:`C_l`, within the gridbox as

.. math::

   C_l = \int_{s=-Q_c}^{\infty} G(s) ds
   \label{eq:int_gs_ds}

and the expression for mean condensate, :math:`{\overline{q_{cl}}}`,
using (`[eq:l_qc_s] <#eq:l_qc_s>`__) to expand :math:`q_{cl}`, is

.. math::

   \overline{q_{cl}} = \int_{s=-Q_c}^{\infty} (Q_c + s) G(s) ds .
   \label{eq:qclbar=int}

If we know (parametrize) the PDF given by :math:`G(s)` then we can solve
for :math:`C_l` and :math:`\overline{q_{cl}}`. Note that this
distribution is in terms of :math:`s`, there is no need to know the
three-dimensional distribution in terms of three separate variables
:math:`q_T`, :math:`T_L` and :math:`p`. This is the method used by
:raw-latex:`\cite{smith90}`, where a symmetric triangular distribution
function is used. For further information on the
:raw-latex:`\cite{smith90}` scheme, please refer to . Physics and
dynamics schemes hence only need to provide increments to
:math:`\overline{q_T}` and :math:`\overline{T_L}`, provided that a
diagnostic scheme (such as :raw-latex:`\cite{smith90}`) is called at
some point in the timestep to partition :math:`\overline{q_T}` into
:math:`\overline{q}` and :math:`\overline{q_{cl}}`, to calculate the dry
bulb temperature :math:`\overline{T}` (from :math:`\overline{T_L}` and
:math:`\overline{q_{cl}}`) and to calculate the liquid cloud fraction,
:math:`C_l`. The diagnostic scheme effectively allows a calculation of
condensation associated with any physical process. However, its results
remain tied to the distribution of :math:`G(s)` that is chosen in
(`[eq:int_gs_ds] <#eq:int_gs_ds>`__) and
(`[eq:qclbar=int] <#eq:qclbar=int>`__) and it is this tie that we seek
to break by the use of a prognostic scheme.

Concept of PC2
--------------

The PC2 scheme develops prognostic expressions for the rates of change
of cloud fraction and condensate contents as a result of each process
that acts in the model. We consider ice and liquid condensate as two
distinct aspects of clouds, which may or may not overlap Figure
`1 <#fig:schematic>`__ provides a schematic summary of the PC2 scheme.
The equations for the five prognostic cloud variables can be written
schematically:

.. math::

   \begin{aligned}
   \frac{\partial \overline{q_{cl}}}{\partial t} =
   \frac{\partial \overline{q_{cl}}}{\partial t} |_{advection} +
   \frac{\partial \overline{q_{cl}}}{\partial t} |_{convection} +
   \frac{\partial \overline{q_{cl}}}{\partial t} |_{boundary \, layer} +
   \frac{\partial \overline{q_{cl}}}{\partial t} |_{precipitation} + ... \nonumber \\
   \frac{\partial \overline{q_{cf}}}{\partial t} =
   \frac{\partial \overline{q_{cf}}}{\partial t} |_{advection} +
   \frac{\partial \overline{q_{cf}}}{\partial t} |_{convection} +
   \frac{\partial \overline{q_{cf}}}{\partial t} |_{boundary \, layer} +
   \frac{\partial \overline{q_{cf}}}{\partial t} |_{precipitation} + ... \nonumber \\
   \frac{\partial C_l}{\partial t} =
   \frac{\partial C_l}{\partial t} |_{advection} +
   \frac{\partial C_l}{\partial t} |_{convection} +
   \frac{\partial C_l}{\partial t} |_{boundary \, layer} +
   \frac{\partial C_l}{\partial t} |_{precipitation} + ... \nonumber \\
   \frac{\partial C_i}{\partial t} =
   \frac{\partial C_i}{\partial t} |_{advection} +
   \frac{\partial C_i}{\partial t} |_{convection} +
   \frac{\partial C_i}{\partial t} |_{boundary \, layer} +
   \frac{\partial C_i}{\partial t} |_{precipitation} + ... \nonumber \\
   \frac{\partial C_t}{\partial t} =
   \frac{\partial C_t}{\partial t} |_{advection} +
   \frac{\partial C_t}{\partial t} |_{convection} +
   \frac{\partial C_t}{\partial t} |_{boundary \, layer} +
   \frac{\partial C_t}{\partial t} |_{precipitation} + ... ,
   \label{eq:dqcldt_and_dcdt}
   \end{aligned}

where :math:`\overline{q_{cf}}` is the ice water specfic humidity,
:math:`C_l` is the liquid cloud *volume* fraction, :math:`C_i` is the
ice cloud volume fraction, and :math:`C_t` is the combined ice or liquid
cloud volume fraction. The amount of mixed phase cloud, :math:`C_{mp}`,
can be calculated by the overlap of the ice and liquid fractions:

.. math::

   C_{mp} = C_i + C_l - C_t.
   \label{eq:mp}

The idea is to parametrize each of the terms in the above equations.
This approach removes the diagnostic method, hence it will be critical
that we can write expressions for
:math:`\frac{\partial \overline{q_{cl}}}{\partial t}` and
:math:`\frac{\partial C_l}{\partial t}` for *each process that alters
:math:`\overline{T}`, :math:`\overline{p}`, :math:`\overline{q}`, or
:math:`\overline{q_{cl}}` in the model* (and similarly for the ice
terms). In doing so, we will not lose sight of underlying PDF approach
given by (`[eq:int_gs_ds] <#eq:int_gs_ds>`__) and
(`[eq:qclbar=int] <#eq:qclbar=int>`__) since we will still use the
concept of instantaneous condensation for liquid clouds. Equations
`[eq:int_gs_ds] <#eq:int_gs_ds>`__ and
`[eq:qclbar=int] <#eq:qclbar=int>`__ will form the basis of the
homogeneous forcing methods discussed in section `3.2 <#sec:homog>`__.
We note in particular that the convective cloud fraction, previously a
quantity that is diagnosed separately from the large-scale cloud
fraction calculated by the :raw-latex:`\cite{smith90}` scheme, may, in
PC2, be included as part of the large-scale cloud fraction. This aspect
is similar to the :raw-latex:`\cite{t93}` approach.

The final aim of PC2 is that the parametrization of each term in
(`[eq:dqcldt_and_dcdt] <#eq:dqcldt_and_dcdt>`__) is performed by each
part of the model that alters :math:`\overline{T}`,
:math:`\overline{p}`, :math:`\overline{q}`, :math:`\overline{q_{cl}}` or
:math:`\overline{q_{cf}}` as an integral part of that physics or
dynamics scheme. However, in this PC2 scheme we acknowledge that this
will not be possible, at least, not to begin with. Hence we have
specifically developed generic approaches that can be used to calculate
expressions for :math:`\frac{\partial \overline{q_{cl}}}{\partial t}`
and :math:`\frac{\partial C_l}
{\partial t}` . These are referred to as Homogeneous forcing (section
`3.2 <#sec:homog>`__), Injection source (or inhomogeneous forcing,
section `3.5 <#sec:inhomog>`__) and Width Changing (section
`3.3 <#sec:width>`__). Two additional modules are available to assist
with PC2, liquid cloud initiaion (section `3.4 <#sec:init>`__) and the
calculation of total cloud fraction changes (section `3.6 <#sec:ct>`__).
At the present time, only the large-scale precipitation (section
`4.2 <#sec:precip>`__) scheme has been rewritten fully to use the PC2
concept of prognostic cloud fractions. The existing mass-flux convection
scheme has been modified to enable calculation of the detrained
condensate, but direct modification to the cloud fraction is not
included. All other physics schemes use one of the generic approaches
below.

A note on convective cloud fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It was the original intention that PC2 be able to replace the two
separate diagnostic cloud fractions (large-scale and convective) with a
single cloud fraction, as in :raw-latex:`\cite{t93}`. The hypothesis was
that by detraining cloud directly from the convection scheme we would no
longer need a separate representation of this cloud type. Our experience
with PC2 is that this is not necessarily the case. We suspect that the
basic reason is that we are unable to truely represent the extreme PDF
shapes that result from convective activity. Additionally, we only
create cloud associated with the detrainment part of the convection
scheme, assuming that cloud associated with the active updraughts in
convection is small. This assumption is not necessarily applicable.
Similar arguments, and model results, come from analysis of the
:raw-latex:`\cite{t93}` and :raw-latex:`\cite{t02}` scheme (Ben Johnson,
personal communication). We also note that with two cloud fraction types
and two different optical depths it is possible to have a basic degree
of representation of cloud inhomogeneity.

Hence the code still exists to enable PC2 to be run with or without a
diagnostic convective cloud fraction, although PC2:66 does not include a
diagnostic term. More details are in section `4.7 <#sec:convec>`__.

Physical basis of the PC2 prognostic cloud scheme
=================================================

In this section we will develop the physical models that PC2 uses in
order to calculate its prognostic increment terms. We will also consider
the numerical solution of the models. The way in which these are
incorporated into the Unifed Model will be discussed in section
`5 <#sec:um>`__

Instantaneous condensation
--------------------------

Liquid clouds in PC2 use the concept of instantaneous condensation.
Hence the ‘s’ distribution methods are fully applicable to the
development of the equations that govern the parametrization of liquid
cloud in PC2. We will start by looking at changes to
:math:`\overline{q_{cl}}` and :math:`C_l` when a uniform forcing is
applied to a gridbox, under the assumption of instantaneous
condensation.

.. _`sec:homog`:

Homogeneous forcing
-------------------

We define the expression *uniform forcing* (or *homogeneous forcing*) to
refer to changes in local values of :math:`T_L` and :math:`q_T` that
occur at a rate independent of the part of the gridbox in which they are
located. This implies that :math:`G(s)` will not alter due to such a
process. Uniform forcing simply alters :math:`Q_c` in
(`[eq:int_gs_ds] <#eq:int_gs_ds>`__) and
(`[eq:qclbar=int] <#eq:qclbar=int>`__). In the Unified Model, this
concept will be applied to several different sets of physics increments
in order to calculate the condensation and cloud fraction changes
associated with each one, where the physics routine does not allow the
explicit calculation of condensation and cloud fraction changes by
another method. Large-scale ascent may be considered a meteorological
example of such a process. By differentiating
(`[eq:int_gs_ds] <#eq:int_gs_ds>`__) and
(`[eq:qclbar=int] <#eq:qclbar=int>`__) with respect to time, assuming
uniform forcing (so :math:`{\frac{\partial G}{\partial t}}` terms are
zero), we obtain

.. math::

   {{\frac{\partial C_l}{\partial t}} = G(-Q_c) {\frac{\partial Q_c}
   {\partial t} }.}
   \label{dcdt}

.. math::

   {\frac{\partial \overline{q_{cl}}}{\partial t}} =
   C_l {\frac{\partial Q_c}{\partial t}}
   \label{dqcldt}

The quantity :math:`G(-Q_c)` is the value of the PDF of :math:`G` at
:math:`s=-Q_c`, which defines the boundary between the saturated and
unsaturated parts of the distribution.

If we wish to consider a prognostic cloud scheme with equations for the
rate of change of condensate and cloud fraction based upon
(`[dqcldt] <#dqcldt>`__) and (`[dcdt] <#dcdt>`__) then we need to close
(`[dcdt] <#dcdt>`__) by specifying the value of :math:`G(-Q_c)`. We will
choose to develop a parametrization for this quantity based upon the
quantities :math:`C_l`, :math:`\overline{q_{cl}}` and the saturation
deficit, :math:`SD`, rather than tie :math:`G(-Q_c)` to a process. The
saturation deficit is *defined* here in the ‘s’ framework to be the
first moment of the PDF for ‘s’ values less than :math:`-Q_c`. In this
way it is analogous to the liquid water content,
:math:`\overline{q_{cl}}`. Appendix A of :raw-latex:`\cite{wg03}` writes
this *definition* as

.. math::

   {SD = - \int_{-\infty}^{-Q_c} {( s+Q_c ) G(s) ds}}
   \label{SD}

and shows this is equivalent to

.. math::

   {SD = a_L  ( q_{sat}({\overline{T}},{\overline{p}}) - {\overline{q}} ) .}
   \label{SD2}

The basis behind the parametrization for :math:`G(-Q_c)` is to consider
an underlying form of the distribution :math:`G(s)` near the
:math:`+b_s` and :math:`-b_s` ends. We borrow the notation of
:raw-latex:`\cite{smith90}` and refer to a quantity :math:`b_s` that is
the value of :math:`s` when a monomodal distribution :math:`G(s)` just
equals zero. We suppose that the distribution G can be described as a
power law near :math:`s=b_s`.

.. math::

   G(s) ~ \propto ~ {(-s + b_s)}^n 
   \label{eqn19}

provided :math:`s<b_s`, where :math:`b_s` represents the ‘width’ of the
distribution (so :math:`G(-Q_c)` = 0 at :math:`-Q_c` = :math:`b_s`) and
n is a power. From (`[eqn19] <#eqn19>`__) it can be shown (see appendix
B of :raw-latex:`\cite{wg03}`) that

.. math::

   {G_1(-Q_c) = {\frac{(n+1)}{(n+2)}} {\frac{C_l^2}{\overline{q_{cl}}}} .}
   \label{eqn20}

An important feature is that the proportionality between :math:`G(-Q_c)`
and :math:`{\frac{C^2}{\overline{l}} }` holds for any power law
description. Also, this relationship is independent of the value of
:math:`b_s`. The triangular :raw-latex:`\cite{smith90}` scheme obeys
this relationship (for :math:`C_l` less than 0.5) with :math:`n`\ =1, as
does a ‘top hat’ function which is a limiting case of :math:`n` tending
to zero. This invariant functional form can be exploited in deriving a
generalized :math:`G(-Q_c`) closure. If we assume a similar power law
for the other end of the distribution we can write a second estimate of
:math:`G(-Q_c)` as:

.. math::

   {G_2(-Q_c) = {\frac{(n+1)}{(n+2)}} {\frac{{(1-C_l)}^2}{SD}} .}
   \label{eqn21}

We note that if n tends to zero then (`[eqn21] <#eqn21>`__) is identical
to the expression used by :raw-latex:`\cite{jgt99}`. This is because
:raw-latex:`\cite{jgt99}` also uses a similar description of a ‘top-hat’
PDF of fluctuations.

In order that the closure of :math:`G(-Q_c)` is reversible, we take a
linear combination of :math:`G_1(-Q_c)` and :math:`G_2(-Q_c)`. To close
our parameterisation, we must choose suitable weights to apply to the
two solutions, and there are currently 2 options for the choice of
weights, discussed below.

The equations (`[dqcldt] <#dqcldt>`__), (`[dcdt] <#dcdt>`__) and either
(`[eqn22] <#eqn22>`__) or (`[eq:gmqc_width] <#eq:gmqc_width>`__) below
form a complete mathematical set for the solution of :math:`C_l` and
:math:`\overline{q_{cl}}` under homogeneous forcing, provided that the
initial value of :math:`C_l` is not identically 0 or 1. If :math:`C_l`
is 0 or 1 then :math:`G(-Q_c)` remains at zero. The equation set then
needs to be initiated in some way. This is discussed further in
:raw-latex:`\cite{wg03}` and in section `3.4 <#sec:init>`__. If
:math:`G(-Q_c)` is defined, then application of the above equation set
may be used to trace out an underlying PDF for any input values of
:math:`C_l`, :math:`\overline{q_{cl}}` and :math:`SD`. Although we never
need to define the complete PDF in PC2 (just the value of
:math:`G(-Q_c)` from equation `[eqn22] <#eqn22>`__), a PDF can be
inferred off-line if required.

:raw-latex:`\cite{wg03}` analyse the performance of this parametrization
under idealised tests, and shows it performs well against other cloud
parametrizations used in large-scale models. After extensive analysis
including observed PDF shapes and moments from balloon, performance in
the Unified Model, and numerical stability tests, the value of the shape
parameter :math:`n` has been chosen to be 0.0, corresponding to a
top-hat distribution shape.

Weight as a function of cloud-fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This option is selected by setting **i_pc2_homog_g_method=1** in the UM
large-scale cloud namelist.

Under this closure, we choose the weights such that :math:`G_1(-Q_c)`
(`[eqn20] <#eqn20>`__) is used when cloud fractions are small, and
:math:`G_2(-Q_c)` (`[eqn21] <#eqn21>`__) is used when cloud fractions
are large (a small cloud fraction indicates that the saturation boundary
is close to the right-hand end of the PDF, which :math:`G_1(-Q_c)` is
based on). We choose relative weights of :math:`{{(1-C_l)}^{m}}` and
:math:`{C_l^m }` respectively, where :math:`m` is a power currently set
to 0.5. Hence the complete suggested closure of :math:`G(-Q_c)` is:

.. math::

   {G(-Q_c) = {\frac{(n+1)}{(n+2)}} {\frac{ ( {( 1-C_l )}^m {\frac{C_l^2}
   {\overline{q_{cl}}}} + C_l^m
   {\frac{{(1-C_l)}^2}{SD}} ) }{( {(1-C_l)}^m + C_l^m ) }} . }
   \label{eqn22}

A problematic property of equation `[eqn22] <#eqn22>`__ is that it goes
to infinity if either :math:`q_{cl}` or :math:`SD` goes to zero (and
:math:`C_l` is not zero or unity). There are realistic scenarios in
which this limit will be approached; e.g. if heavy rain falls through a
layer of cloud, nearly all of the cloud liquid water content maybe
removed by accretion, without reducing the cloud-fraction. In this
situation, any subsequent homogeneous forcing applied to the cloud will
result in a huge tendency in cloud-fraction in equation
`[dcdt] <#dcdt>`__, due to the term
:math:`\frac{C_l^2}{\overline{q_{cl}}}` in `[eqn22] <#eqn22>`__ becoming
huge.

Note that for a homogeneous forcing acting to dry the layer / reduce the
cloud, a huge negative tendency is the “right” answer; if there is only
an infinitessimally small amount of liquid water content left within the
cloud, then the cloud fraction should indeed vanish extremely rapidly
under a negative forcing. However, for a positive homogeneous forcing,
the cloud-fraction will very rapidly increase in this scenario, for no
physical reason. This might not be a problem if one exactly integrated
the differential equations `[dcdt] <#dcdt>`__ and
`[dqcldt] <#dqcldt>`__, since :math:`q_{cl}` would immediately increase
away from zero, so that `[eqn22] <#eqn22>`__ immediately becomes
well-defined (an infinitely large tendency maintained for an infinitely
small period of time can yield a finite, sensible increment!)
Unfortunately, PC2 uses an explicit numerical method, and a finite
(often quite large) model timestep, so an instantaneous very large
tendency will yield a very large increment, even if the continuous
differential equations would not have done.

In practice, the code that implements `[eqn22] <#eqn22>`__ simply sets
:math:`G(-Q_c)` to zero if either :math:`q_{cl}` or :math:`SD` falls
below 1.0E-10 kg kg\ :math:`^{-1}`, to avoid a floating-point error.
This in itself can be problematic, since leaving :math:`C_l` unmodified
under a homogeneous forcing can allow unrealistic states to develop.

Weight in proportion to PDF width
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This option is selected by setting **i_pc2_homog_g_method=2** in the UM
large-scale cloud namelist.

:raw-latex:`\cite{morcrette_2020}` proposed an alternative choice of
weights applied to :math:`G_1(-Q_c)` (`[eqn20] <#eqn20>`__) and
:math:`G_2(-Q_c)` (`[eqn21] <#eqn21>`__), so-as to make each one’s
weight go to zero in the limit that it goes to infinity, reliably
yielding a sensible, finite solution for :math:`G(-Q_c)`.

We choose the weights to be :math:`\frac{\overline{q_{cl}}}{C_l}`
applied to :math:`G_1(-Q_c)`, and :math:`\frac{SD}{1-C_l}` applied to
:math:`G_2(-Q_c)`, yielding:

.. math::

   G(-Q_c) = \frac{(n+1)}{(n+2)} \frac{
        \frac{\overline{q_{cl}}}{C_l} {\frac{C_l^2}{\overline{q_{cl}}}}
     + \frac{SD}{1-C_l} {\frac{{(1-C_l)}^2}{SD}}
             }{ \frac{\overline{q_{cl}}}{C_l} + \frac{SD}{1-C_l} }

This is equivalent to weighting the saturated and subsaturated solutions
for :math:`G(-Q_c)` by their respective PDF-widths. Note that everything
cancels-out in the numerator, so this reduces to:

.. math::

   G(-Q_c) = \frac{(n+1)}{(n+2)} \frac{ 1
             }{ \frac{\overline{q_{cl}}}{C_l} + \frac{SD}{1-C_l} }
   \label{eq:gmqc_width}

In full UM tests, :raw-latex:`\cite{morcrette_2020}` found that this
closure removed occasional spurious very large cloud-fraction increments
that occur when using `[eqn22] <#eqn22>`__, but did not significantly
impact the performance of the model forecast.

.. _`sec:homog_num_app`:

Numerical application
~~~~~~~~~~~~~~~~~~~~~

The timestepping methods that are used in the homogeneous forcing were
developed off-line using a single gridbox model, to ensure smooth,
accurate and convergent behaviour of the solution, rather than as a
result of a mathematical analysis of the problem.

We need to timestep forward (`[dqcldt] <#dqcldt>`__) and
(`[dcdt] <#dcdt>`__) using a timestep of :math:`\Delta{t}`, knowing
values of :math:`\Delta{\overline{T}}`, :math:`\Delta{p}`,
:math:`\Delta{\overline{q}}` and :math:`\Delta{\overline{q_{cl}}}`,
which are provided by an existing physics scheme in the model. We assume
that the physics scheme we are applying this to has not already
calculated its condensation and cloud fraction increments by another
means.

Firstly, we need to calculate the forcing term :math:`\Delta{Q_c}`. From
(`[eq:qc_eq_qt-qs] <#eq:qc_eq_qt-qs>`__) we can write

.. math::

   \Delta{Q_c} = a_L ( \Delta{\overline{q_T}} - \Delta{\overline{q_{sat}(T_L)}} )
   \label{eq:deltaqc}

assuming that :math:`a_L` does not change (see below). This can be
expanded, using a linear approximation for
:math:`q_{sat}(T+\Delta{T},p+\Delta{p})` in terms of
:math:`q_{sat}(T,p)` as

.. math::

   \Delta{Q_c} = a_L ( \Delta{\overline{q}} + \Delta{\overline{q_{cl}}}
   - \alpha \Delta{\overline{T_L}}
   - \beta \Delta{\overline{p}} )
   \label{eq:deltaqc_exp}

where :math:`\alpha` is the rate of change of :math:`q_{sat}` with
respect to temperature at constant pressure
(:math:`\frac{\partial q_{sat}}{\partial T}`), and :math:`\beta` is the
rate of change of :math:`q_{sat}` with respect to pressure at constant
temperature (:math:`\frac{\partial q_{sat}}{\partial p}`). Using
(`[eq:tl] <#eq:tl>`__) to expand :math:`T_L` in terms of :math:`T` and
:math:`q_{cl}` and gathering terms together we obtain

.. math::

   \Delta{Q_c} = a_L ( \Delta{\overline{q}} - \alpha \Delta{\overline{T}}
   - \beta \Delta{\overline{p}} ) +  \Delta{\overline{q_{cl}}} .
   \label{eq:deltaqc_exp2}

We must now note the method we use here to calculate :math:`a_L`,
defined in (`[eq:a_L] <#eq:a_L>`__), uses a gradient expansion value of
:math:`\alpha` that corresponds to :math:`\frac{\partial{q_{sat}}}
{\partial{T}}` at constant pressure and not the chord expression in
(`[eq:alpha] <#eq:alpha>`__). This is because we are trying to find the
best estimate of the increments, not the absolute value. We have seen
errors arise in simple numerical tests when the chord expression is
used. We need to define the temperature around which this calculation of
:math:`\frac{\partial{q_{sat}}}{\partial{T}}` is made. Throughout PC2 we
choose the dry-bulb temperature :math:`\overline{T}`, and not the
liquid-temperature (:math:`\overline{T_L}`), since :math:`q_{sat}`
locally is defined by the local dry-bulb temperature (:math:`T`) and we
need to consider *changes* in the condensate. This has been confirmed
using simulations using a single gridbox model. contains a longer
discussion of this issue, but we note here that the
:raw-latex:`\cite{smith90}` scheme performs best when it does not use
:math:`T` to calculate :math:`\alpha` but the gradient of the chord
between :math:`(\overline{T_L}, q_{sat}(\overline{T_L}))` and
:math:`(\overline{T}, q_{sat}(\overline{T}))`, as in
(`[eq:alpha] <#eq:alpha>`__). The best choice is dependent on the method
of implementation. PC2 uses (`[eq:a_L] <#eq:a_L>`__) with the standard
thermodynamic relationships (e.g. :raw-latex:`\cite{ry89}`, chapter 2)

.. math::

   \alpha = \frac{ \epsilon L q_{sat}(\overline{T}) } { R \overline{T}^2} ,
   \label{eq:alpha_exp}

where :math:`R` is the gas constant for dry air, and

.. math::

   \beta = \frac{-q_{sat}(\overline{T})}{\overline{p}} .
   \label{eq:beta}

The right hand side of (`[eq:deltaqc_exp2] <#eq:deltaqc_exp2>`__) now
contains forcing values which we know from the physics scheme we are
applying the homogeneous forcing to.

We next estimate :math:`G(-Q_c)` from the parametrization
(`[eqn22] <#eqn22>`__) and the expression for the saturation deficit
(`[SD2] <#SD2>`__). The change in :math:`C_l` is then estimated using a
simple forward step of (`[dcdt] <#dcdt>`__) using
(`[eq:deltaqc_exp2] <#eq:deltaqc_exp2>`__):

.. math::

   \Delta{C_l} = G(-Q_c) \Delta{Q_c} .
   \label{eq:deltac}

The final value of :math:`C_l` is then limited to lie between 0 and 1.

.. math::

   C_l^{[n+1]} = (0, ~ C_l^{[n]} + \Delta{C_l}, ~ 1)
   \label{eq:c_l^n+1}

where :math:`[n]` and :math:`[n+1]` label the timesteps. The
timestepping of :math:`\overline{q_{cl}}` is more involved. We will use
a mid-timestep estimate of :math:`C_l` in the discrete form of
(`[dqcldt] <#dqcldt>`__).

.. math::

   \overline{q_{cl}}^{[n+1]} = \overline{q_{cl}}^{[n]} + \frac{1}{2}
   (C_l^{[n]} + C_l^{[n+1]}) \Delta{Q_c}.
   \label{eq:qcl_l^n+1}

This completes the homogeneous forcing routine. We note that it is
possible for :math:`\overline{q_{cl}}^{[n+1]}` to be negative if the
forcing is strong enough. Originally, it was not deemed desirable to
prevent homogeneous forcing processes from doing this, in order not to
interfere with the possible cancellation of positive and negative
increments from different physics schemes. A checking routine is
applied, however, in the Unified Model to remove any negative values
that are generated, which is discussed in section
`4.10 <#sec:checks>`__. However, the checking routine (Q-Pos) involves a
lot of communication between processors and can significantly increase
the run-time of the model. The option to “Ensure consistent sinks of qcl
and CFL” performs a check at the end of the homogeneous forcing routines
to ensure that we are not trying to remove more condensate than was
there to start with.

Note that (`[eq:c_l^n+1] <#eq:c_l^n+1>`__) and
(`[eq:qcl_l^n+1] <#eq:qcl_l^n+1>`__) *include* the contribution of the
forcing itself, it is not just the reactionary condensation. If we wish
to isolate the condensation associated with the forcing, then we must
subtract any liquid forcing from the final solution. The net change in
:math:`\overline{q}` as a result of the homogeneous forcing is simply
the net change in :math:`\overline{q_T}` minus the net change in
:math:`\overline{q_{cl}}`. The net change in :math:`\overline{T}` is
calculated simply to account for the latent heat released due to the
condensation.

.. _`sec:width`:

Changing the width of the PDF - PC2 erosion
-------------------------------------------

Another basic change to the PDF that can be mathematically analysed is
if the PDF shape is kept constant but its width (and therefore height)
is altered. This is, perhaps, the simplest method of representing a
process that changes the shape of the PDF, and we will, in PC2, apply it
to represent mixing of air within a gridbox, although this is a
significant approximation of the process. Its application fulfils the
role of the “cloud erosion” term in the :raw-latex:`\cite{t93}` scheme.
By linking the term to the PDF shape we can place this term on a
stronger mathematical footing than the simple reduction term
parametrized by :raw-latex:`\cite{t93}`. Equivalent arguments enabled
:raw-latex:`\cite{ww99}` to retrieve the same result as presented here.

We can consider a change in the width of the PDF to alter its form
according to

.. math::

   G^{[n+1]}(s) = \xi G^{[n]} (\xi s)
   \label{eq:g_xi}

where :math:`G^{[n+1]}(s)` is the distribution after the change in
width, :math:`G^{[n]} (s)` is the distribution before the change in
width and :math:`\xi` is a scaling factor. If :math:`\xi > 1` then the
distribution is narrowed. For the liquid cloud fraction we therefore
have

.. math::

   C_l^{[n+1]} = \int_{s=-Q_c}^{\infty} \xi G(\xi s) ds .
   \label{eq:c_l_xi}

If we transform variables to :math:`s' = \xi s` we can rewrite this
integral as

.. math::

   C_l^{[n+1]} = \int_{s'=-Q_c \xi}^{\infty} G(s') ds' .
   \label{eq:c_l_xi2}

Hence the expression for :math:`C_l^{[n+1]}` is equivalent to using the
same distribution function :math:`G(s)` as for :math:`C_l^{[n]}` except
that the saturation boundary has been moved from :math:`-Q_c` to
:math:`-Q_c \xi`. The result is the same as applying a homogeneous
forcing (`[eq:deltac] <#eq:deltac>`__) with a modified forcing,

.. math::

   \Delta Q_c \equiv \xi Q_c - Q_c ,
   \label{eq:deltac_modified}

or the continuous version

.. math::

   \frac{\partial Q_c}{\partial t} \equiv
   Q_c \frac{\partial}{\partial t}(\xi - 1) .
   \label{eq:xi_equiv}

We can write :math:`\xi` in a slightly more informative way by linking
it to the relative change in width of the PDF
:math:`\frac{1}{b_s} \frac{\partial b_s}{\partial t}`. For a PDF that
changes its width, :math:`\xi` is defined as

.. math::

   \xi = \frac{b_s}{b_s + \delta b_s} = \frac{1}{1 + \frac{\delta b_s}{b_s}}.
   \label{eq:xi_equiv1}

For an infintessimal timestep :math:`\delta t` we therefore have

.. math::

   \xi = \frac{1}{1 + \frac{1}{b_s} \frac{\partial b_s}{\partial t} \delta t }
   \label{eq:xi}

and hence, by expanding (`[eq:xi] <#eq:xi>`__) to give
:math:`\xi = 1 -  \frac{1}{b_s}
\frac{\partial b_s}{\partial t} \delta t` and using the homogeneous
forcing expression (`[dcdt] <#dcdt>`__) with the modified forcing
(`[eq:xi_equiv] <#eq:xi_equiv>`__), we retrieve the continuous form

.. math::

   \frac{\partial C_l}{\partial t} = - G(-Q_c) Q_c \frac{1}{b_s}
   \frac{\partial b_s}{\partial t} .
   \label{eq:dcdt_width}

A similar analysis can be performed for
:math:`\frac{\partial \overline{q_{cl}}}
{\partial t}` from (`[dqcldt] <#dqcldt>`__) to give

.. math::

   \overline{q_{cl}}^{[n+1]} = \frac{1}{\xi} \int_{s'=-Q_c \xi}^{\infty}
   (- \xi Q_c + s') G(s') ds' .
   \label{eq:qcl_xi2}

Again, this is equivalent to using the homogeneous forcing with the
modified forcing (`[eq:xi_equiv] <#eq:xi_equiv>`__), but it also
includes a scaling term :math:`\frac{1}{\xi}`. In the infinitessimal
limit, this scaling gives a second term that is proportional to the
value of the integral (i.e. :math:`\overline{q_{cl}}`). Hence we obtain
the final continuous solution

.. math::

   \frac{\partial \overline{q_{cl}}}{\partial t} =
   (- C_l Q_c+\overline{q_{cl}}) \frac{1}{b_s} \frac{\partial b_s}{\partial t} .
   \label{eq:dqcldt_width}

To close the solution, we need to parametrize :math:`\frac{1}{b_s}
\frac{\partial b_s}{\partial t}` , which could be linked to the physics
of the process that is occuring. Note we don’t need to calculate
:math:`b_s` separately, just its *fractional* rate of change. Options
for the parameterisation of
:math:`\frac{1}{b_s}\frac{\partial b_s}{\partial t}` due to turbulent
“erosion” are described in section `4.3 <#sec:turb>`__, along with the
numerical methods used to integrate the equations.

.. _`sec:init`:

Initiation of cloud
-------------------

In section `3.2 <#sec:homog>`__ we commented that the closure
(`[eqn22] <#eqn22>`__) for :math:`G(-Qc)` is only valid if :math:`C_l`
is not identically 0 or 1. If :math:`C_l` is 0 or 1 we know that
:math:`G(-Q_c)` is equal to 0 but we have lost the information that will
tell us when :math:`G(-Q_c + \Delta Q_c)` starts differing from 0. Hence
the homogeneous forcing equation set (`[dcdt] <#dcdt>`__),
(`[dqcldt] <#dqcldt>`__) and (`[eqn22] <#eqn22>`__) is not complete if
we start from a position where :math:`C_l` is 0 or 1. To complete this
set, we will need to define a width, :math:`b_s`, to the PDF and provide
an initiation increment to :math:`C_l` and :math:`\overline{q_{cl}}`
when the value of :math:`-Q_c` crosses the limit of the distribution.
There is more discussion in :raw-latex:`\cite{wg03}`.

To initiate new partial cloud-cover (or new partial clear-sky), we
essentially call a diagnostic cloud scheme to initialise the prognostics
:math:`C_l` and :math:`\overline{q_{cl}}`. In the UM there is currently
a choice of 2 different diagnostic cloud schemes that can be used for
this; either a version of the Smith scheme (see UMDP 029), or the
bimodal scheme (see UMDP 039). These two options are described below...

Initiation using a “Smith-like” method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This option is selected by setting the UM namelist switch
**i_pc2_init_method = 1** (Smith).

We will assume the same form of the PDF at its boundaries as is assumed
in the derivation of the :math:`G(-Q_c)` closure. For the high
‘:math:`s`’ end of the PDF distribution we integrate the power law
description in (`[eqn19] <#eqn19>`__) to obtain the expressions

.. math::

   C_l = \frac{1}{2 b_s^{n+1}} (b_s + Q_c)^{n+1} ,
   \label{eq:initc}

.. math::

   \overline{q_{cl}} = \frac{1}{2 b_s^{n+1}} \frac{(b_s + Q_c)^{n+2}}{n+2} .
   \label{eq:initqcl}

We now need to parametrize the PDF width :math:`b_s`. Unlike the
:raw-latex:`\cite{smith90}` scheme, this is the only location in the PC2
cloud scheme where the width needs to be defined for the liquid cloud
(although see section `4.2.4 <#sec:mp_depsub>`__ for a discussion of an
equivalent width in the deposition / sublimation relationship for ice
cloud). We still choose to define :math:`b_s` in terms of a critical
relative humidity parameter, :math:`RH_{crit}`. Like the
:raw-latex:`\cite{smith90}` scheme (see ), we define the value of
:math:`b_s` as

.. math::

   b_s = a_L (1 - RH_{crit}) q_{sat} (\overline{T_L}) .
   \label{eq:bs}

Hence, if the parameter :math:`n` was the same in PC2 as the equivalent
in :raw-latex:`\cite{smith90}`, the initial creation of liquid cloud
would follow precisely that diagnosed by the :raw-latex:`\cite{smith90}`
scheme (assuming that the numerical implementation of the calculation is
the same). Its subsequent behaviour in PC2, though, would be different,
because the subsequent physical processes that act are parametrized in
different ways. Note: for some reason, the implementation in the UM uses
a fixed value of :math:`n = 0` (corresponding to a top-hat distribution)
if a constant :math:`RH_{crit}` profile is used, but instead sets
:math:`n = 1` (a triangular distribution) in the PC2 initiation
calculation if a TKE-based variable :math:`RH_{crit}` is used. In the
latter case, :math:`n = 0` is still hardwired in the PC2 homogeneous
forcing calculations, so it is not handled consistently.

An equivalent initiation scheme is required if :math:`C_l` is 1 and
:math:`Q_c` is being reduced - at some point we need to introduce clear
sky into the solution. Because we make the choice of symmetry (which
could be relaxed if we used different :math:`RH_{crit}` values for
:math:`C_l` of 1 and :math:`C_l` of 0), the problem is entirely
equivalent to that of initiating from :math:`C_l = 0`, with the
exception that :math:`\overline{q_{cl}}` is replaced by :math:`SD`,
:math:`C_l` is replaced by :math:`(1-C_l)`, and :math:`Q_c` is replaced
by :math:`-Q_c`. We hence have the solution

.. math::

   1 - C_l = \frac{1}{2 b_s^{n+1}} (b_s - Q_c)^{n+1} ,
   \label{eq:init1mc}

.. math::

   SD = \frac{1}{2 b_s^{n+1}} \frac{(b_s - Q_c)^{n+2}}{n+2} .
   \label{eq:initSD}

The conversion between :math:`SD` and :math:`\overline{q_{cl}}` follows
(`[SD2] <#SD2>`__). We will choose, as we do throughout PC2, to define
:math:`\alpha` (and hence :math:`a_L`) in terms of
:math:`\frac{\partial q_{sat}(\overline{T})}{\partial t}`, although
within this diagnostic calculation of SD it might actually be better to
use the representation (`[eq:alpha] <#eq:alpha>`__) used by the
diagnostic :raw-latex:`\cite{smith90}` scheme.

.. _`sec:numapp_init`:

Numerical Application of the Smith method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to calculate and compare the state of the model to :math:`b_s`,
we first calculate :math:`T_L`, :math:`q_{sat}(\overline{T_L})` and
calculate the mean relative total humidity, :math:`RH_T`, where

.. math::

   RH_T = \frac{ \overline{q} + \overline{q_{cl}} } {q_{sat}(\overline{T_L}) } .
   \label{eq:rht}

We then assess whether initiation is required. There are only two
circumstances in which we wish to proceed further:

- If the current cloud fraction :math:`C_l` is 0 and :math:`-Q_c < b_s`.
  By dividing the second condition by
  :math:`a_L q_{sat} (\overline{T_L})` we see, using the definitions
  (`[eq:qc_eq_qt-qs] <#eq:qc_eq_qt-qs>`__) and (`[eq:bs] <#eq:bs>`__),
  that this second condition is equivalent to :math:`RH_T > RH_{crit}`.

- If the current cloud fraction :math:`C_l` is 1 and :math:`-Q_c > -b_s`
  (or, equivalently, :math:`RH_T < 2 - RH_{crit})`.

Note: in the UM implementation, the actual conditions for when
initiation may occur are more complicated than this, and there are
several options depending on a namelist switch. See section
`4.9 <#sec:init2>`__ for details...

In the second case, we then make the temporary transformation of
variables in order to use the same solution set as in the first case:
:math:`C_l'` takes the value :math:`(1-C_l)` and :math:`RH_t'` takes the
value :math:`(2-RH_t)` (which is equivalent to the replacing of
:math:`Q_c` by :math:`-Q_c`). In the first case, :math:`C_l'` and
:math:`RH_t` take the same values as :math:`C_l` and :math:`RH_t`
respectively.

We then solve for the initiated cloud fraction :math:`C_l'`, using the
similar methods as described in , except that we allow the solution to
vary with the PDF shape :math:`n`. We first write :math:`Q_N` as

.. math::

   Q_N = \frac{Q_c}{b_s} = \frac{ a_L (\overline{q_T} - q_{sat}(\overline{T_L})) }
   { a_L (1 - RH_{crit}) q_{sat} (\overline{T_L})} = \frac{RH_T - 1}{1-RH_{crit}}
   \label{eq:qn_def}

and then use :math:`Q_N` to solve the initiated cloud fraction. We
assume a PDF described by a power law as in (`[eqn19] <#eqn19>`__) (and
the equivalent for the other end of the distribution, the two
expressions switching at :math:`Q_c=0`), which is normalized. The
solution to (`[eq:int_gs_ds] <#eq:int_gs_ds>`__) is hence

.. math::

   C_l^{init'} = \left\{  \begin{array}{ll}
                  0,                                              &  Q_N \le -1      \\
                  \frac{1}{2} {\left( 1 + Q_N \right)}^{n+1},     &  -1 < Q_N \le 0  \\
                  1 - \frac{1}{2} {\left( 1 - Q_N \right)}^{n+1}, &  0 < Q_N < 1     \\
                  1,                                              &  1 \le Q_N .
                \end{array} \right.
   \label{eq:c_qn}

where :math:`C_l^{init'}` is the initiated value of liquid cloud
fraction. If we had performed the variable transformation we then we
need to transform back, so :math:`C_l^{init} = 1 - C_l^{init'}`,
otherwise :math:`C_l^{init} = C_l^{init'}`.

In practice, it is likely to be only the second of the options in
(`[eq:c_qn] <#eq:c_qn>`__) that the scheme uses, since we will be at
that end of the distribution function, unless previous parts of the
model timestep have resulted in large forcings to :math:`Q_c`.

The solution for the initiated liquid water,
:math:`\overline{q_{cl}}^{init}` is more difficult, since it depends on
the width of the distribution :math:`b_s`, hence on :math:`a_L` and
:math:`\alpha`, and :math:`\alpha` is a function of the dry-bulb
temperature :math:`\overline{T}`, which is not known until we know the
amount of condensation. Hence we will need to iterate to a solution.

We first calculate :math:`q_{sat}(\overline{T})`, :math:`\alpha`,
:math:`a_L` and :math:`b_s`, using (`[eq:alpha_exp] <#eq:alpha_exp>`__),
(`[eq:a_L] <#eq:a_L>`__) and (`[eq:bs] <#eq:bs>`__). We then solve for
the liquid water content:

.. math::

   \frac{\overline{q_{cl}}^{init'}}{b_s} = \left\{  \begin{array}{ll}
                     0,                                            &  Q_N \le -1        \\
                     \frac{1}{2 (n+2)} {\left( 1 + Q_N \right)}^{n+2},       &  -1 < Q_N \le 0    \\
                     Q_N + \frac{1}{2 (n+2)} {\left( 1 - Q_N \right)}^{n+2}, &  0 < Q_N < 1       \\
                     Q_N,                                          &  1 \le Q_N .
                   \end{array} \right.
   \label{eq:l_bar}

If we have been working in transformed variables we now transform back,
so the initiated saturation deficit, :math:`SD^{init}`, takes the value
of :math:`\overline{q_{cl}}^{init'}`. We then use (`[SD2] <#SD2>`__) to
estimate :math:`\overline{q_{cl}}^{init}` using our initial estimates of
:math:`q_{sat}(\overline{T})` and :math:`a_L`. If we are not in
transformed variables, we have the first estimate
:math:`\overline{q_{cl}}^{init}=\overline{q_{cl}}^{init'}`.

We now use this estimate of :math:`\overline{q_{cl}}^{init}` to
calculate a more accurate estimate of :math:`a_L` etc. by iteration. In
order to achieve a faster convergence of the iteration, we do not use
(`[eq:l_bar] <#eq:l_bar>`__) directly in the estimation of :math:`a_L`
etc., but use a combination of this value and the one from the previous
iteration.

.. math::

   \overline{q_{cl}}^{init~[i+1]} = f  \overline{q_{cl}}^{init~[i]} +
   (1 - f) \overline{q_{cl}}^{init~[i-1]}
   \label{eq:iter}

where the superscript :math:`[i]` labels each iteration. We find that 10
iterations is effective for convergence, with the weighting :math:`f`
given by :math:`a_L^{[i]}`.

.. _`sec:bimodal_init`:

Initiation using the bimodal scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This option is selected by setting the UM namelist switch
**i_pc2_init_method = 2** (Bimodal).

First, the diagnosis of entrainment zones is performed at all
grid-points, as described in UMDP 039. The parameters of the moisture
PDF are then constructed, assuming either a sum of two Gaussian modes
from the top and bottom of an inversion layer (when within an
entrainment zone), or a single symmetric Gaussian mode (when not in an
entrainment zone). The variance of each Gaussian mode is estimated based
on the TKE and other information output by the boundary-layer scheme
(with a minimum limit applied to the PDF width, consistent with
:math:`RH_{crit}` = 99%. Crucially, each Gaussian mode is truncated to
zero at plus and minus 3 standard deviations; this sets the overall
width of the moisture PDF at each point.

The positions of the upper and lower truncated bounds of the moisture
PDF relative to the saturation threshold are expressed in terms of a
normalised :math:`Q_N` = :math:`Q_c` over PDF-width (see equation
`[eq:qn_def] <#eq:qn_def>`__). In entrainment zones, the sum of the two
Gaussian modes can lead to a highly skewed distribution; hence
:math:`Q_N` can have different values for the upper and lower bounds,
each normalised by the different widths on either side of the PDF. The
upper and lower values of :math:`Q_N` are then compared to -1 and 1
respectively, to determine whether the saturation boundary lies within
the PDF bounds. This is the basic condition for initiation to occur
(though there are additional conditions and various options for these in
the soure code; see section `4.9 <#sec:init2>`__).

If the initiation conditions are met, the diagnostic bimodal cloud
scheme code is then called (see UMDP 039), and the diagnosed :math:`C_l`
and :math:`q_{cl}` are used to set the prognostic :math:`C_l` and
:math:`q_{cl}`.

.. _`sec:inhomog`:

Injection forcing
-----------------

Injection forcing (sometimes referred to as inhomogeneous forcing) uses
another concept of how the underlying moisture PDF may change in order
to calculate a change in cloud fraction as a result of a known injection
of condensate into a gridbox. The term was developed in order to be
coupled with a modified mass-flux convection scheme, but is first
presented here in its basic form.

We will assume a physical model whereby saturated air, containing
condensate, randomly replaces already existing air in the gridbox. (Such
a formulation is designed to represent air detrained from convection
replacing pre-existing air when averaged over a large horizontal
domain). :raw-latex:`\cite{bwg03}` discusses the situation in more
detail. Briefly, we consider two parts to the distribution function
:math:`G(s)`. One part represents the background air. This maintains its
PDF shape (in terms of absolute :math:`q_T` and :math:`T_L` values)
because we assume it is *randomly* replaced, but will reduce in
amplitude as it is replaced by a second PDF representing the injected
air.

The fractional rate at which existing air is replaced by the injected
source air we will write as :math:`\frac{\partial{C_S}}{\partial{t}}`.
Provided that only the liquid phase exists (see section
`3.5.1 <#sec:multiple>`__ for the extention to multiple phases), we then
note that the rate of change of liquid cloud fraction and liquid water
content in the gridbox can be written in two parts: firstly the change
due to the background, and secondly the change due to the source.

.. math::

   \frac{\partial{C_l}}{\partial{t}} =
   - C_l \frac{\partial{C_S}}{\partial{t}} + \frac{\partial{C_S}}{\partial{t}}
   \label{eq:dcdt_inhom}

.. math::

   \frac{\partial{\overline{q_{cl}}}}{\partial{t}} =
   - \overline{q_{cl}} \frac{\partial{C_S}}{\partial{t}}
   + q_{cl}^S \frac{\partial{C_S}}{\partial{t}}
   \label{eq:dqcldt_inhom}

where :math:`q_{cl}^S` is the liquid water content of the injected air.
Eliminating :math:`\frac{\partial{C_S}}{\partial{t}}` gives the
relationship

.. math::

   \frac{\partial{C_l}}{\partial{t}} = \frac{1 - C_l}{q_{cl}^S
   - \overline{q_{cl}}} Q4_l,
   \label{eq:dcdt_inhom2}

where :math:`Q4_l` is the net (*including* the liquid water in the
background distribution that was randomally replaced) injection source
change of :math:`\overline{q_{cl}}`:

.. math::

   Q4_l = \frac{\partial{\overline{q_{cl}}}}{\partial{t}} |_{injection \, source}.
   \label{eq:q4}

We see that we do not need to know anything about the nature of the two
PDFs involved, except the assumption that the injected PDF contains
completely cloudy air. This equation allows one to calculate the change
in :math:`C_l` associated with an injection source change of
:math:`\overline{q_{cl}}` for the example of convection. Modifications
to the mass-flux convection scheme for PC2 (far from trivial and
discussed in depth in section `4.7 <#sec:convec>`__) allow :math:`Q4_l`
to be calculated (:math:`q_{cl}^S` is already available), and
(`[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__) can then be used to calculate
the equivalent :math:`C_l` change. We note at this stage that the
denominator in (`[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__), being the
difference in two terms that may be close to each other, may cause
problems when we attempt to numerically apply this equation.

It is reasonable to ask what happens to the air in the distribution that
was replaced. In this mathematical representation of a single gridbox we
need not know anything other than that the air is displaced into a
neighbouring gridbox. In practical use with a mass-flux convection
scheme we know more that this air is displaced downwards in the column.
We could reasonably calculate the change in cloud fraction following the
same methods as used to calculate the change in :math:`\overline{q}` or
the change in a tracer and we discuss this later.

.. _`sec:multiple`:

Multiple phases in the injection source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The injection source formulation can be extended to multiple phases of
condensate. In practice, this will simply be the two phases ice and
liquid, although we need to recognize that they can overlap with each
other. :raw-latex:`\cite{wilson2001}` provides the background to the
derivation and it is briefly presented below.

We firstly rewrite (`[eq:dqcldt_inhom] <#eq:dqcldt_inhom>`__) but use
the net condensate
(:math:`\overline{q_c} = \overline{q_{cl}} + \overline{q_{cf}}`) instead
of just the liquid water expression, and the net cloud amount
:math:`C_t`, instead of the liquid cloud amount :math:`C_l`. The same
argument as before leads to the expressions

.. math::

   \frac{\partial{C_t}}{\partial{t}} =
   - C_t \frac{\partial{C_S}}{\partial{t}} + \frac{\partial{C_S}}{\partial{t}}
   \label{eq:dctdt_inhom}

and

.. math::

   \frac{\partial{\overline{q_{c}}}}{\partial{t}} =
   - \overline{q_{c}} \frac{\partial{C_S}}{\partial{t}}
   + q_{c}^S \frac{\partial{C_S}}{\partial{t}} .
   \label{eq:dqcdt_inhom}

The left hand side of (`[eq:dqcdt_inhom] <#eq:dqcdt_inhom>`__) is
written as :math:`Q4_c`. :math:`q_{c}^S` is the in-cloud condensate
content (ice plus liquid) of the source.

Hence eliminating :math:`\frac{\partial{C_S}}{\partial{t}}` we obtain

.. math::

   \frac{\partial{C_t}}{\partial{t}} = \frac{(1-C_t)}{q_{c}^S -
   \overline{q_{c}}} Q4_c .
   \label{eq:dctdt_q4}

We will assume that the proportion of the injected volume that contains
liquid cloud can be written as :math:`g_l`, and the proportion that
contains ice cloud can be written as :math:`g_i`. Note that it is not
necessary to have :math:`g_l + g_i = 1` if there is mixed phase cloud
injected. We can write the change in *liquid* cloud fraction
equivalently to (`[eq:dcdt_inhom] <#eq:dcdt_inhom>`__) as

.. math::

   \frac{\partial{C_l}}{\partial{t}} =
   - C_l \frac{\partial{C_S}}{\partial{t}} + g_l \frac{\partial{C_S}}{\partial{t}} .
   \label{eq:dcldt_inhom}

Combining (`[eq:dcldt_inhom] <#eq:dcldt_inhom>`__) and
(`[eq:dctdt_inhom] <#eq:dctdt_inhom>`__) by eliminating
:math:`\frac{\partial{C_S}}{\partial{t}}` gives

.. math::

   \frac{\partial{C_l}}{\partial{t}} = \frac{g_l - C_l}{1 - C_t}
   \frac{\partial{C_t}}{\partial{t}}
   \label{eq:dctdt_dcdt}

and hence from (`[eq:dctdt_q4] <#eq:dctdt_q4>`__) we have the result

.. math::

   \frac{\partial{C_l}}{\partial{t}} =
   \frac{g_l - C_l}{q_c^S - \overline{q_{c}}} Q4_c .
   \label{eq:dcltdt_almost_final}

An equivalent expression holds for the ice cloud. Hence the change in
the amount of cloud for each phase may be calculated assuming we know
the volume proportions of the source term that contain each of the
phases and the net increase in the amount of condensate, :math:`Q4_c`
(regardless of phase). This expression is coded for use in a generically
available inhomogeneous forcing module. However, we can also write this
in a slightly more accessible form by noting the ratio of
(`[eq:dqcdt_inhom] <#eq:dqcdt_inhom>`__) and
(`[eq:dqcldt_inhom] <#eq:dqcldt_inhom>`__) with the :math:`Q4`
definitions following (`[eq:q4] <#eq:q4>`__).

.. math::

   \frac{Q4_c}{q_c^S - \overline{q_c}} =
   \frac{Q4_l}{q_{cl}^S - \overline{q_{cl}}} .
   \label{eq:q4_ratios}

Using (`[eq:q4_ratios] <#eq:q4_ratios>`__) in
(`[eq:dcltdt_almost_final] <#eq:dcltdt_almost_final>`__) gives the final
expression

.. math::

   \frac{\partial{C_l}}{\partial{t}} =
   \frac{g_l - C_l}{q_{cl}^S - \overline{q_{cl}}} Q4_l
   \label{eq:dctdt_final}

and similarly for the ice. Note that this expression accounts for the
possibility that liquid cloud is displaced from the gridbox by added ice
cloud. We can further write :math:`q_{cl}^S` as a fraction of
:math:`q_{c}^S`

.. math::

   q_{cl}^S = h_l q_{c}^S
   \label{eq:qcls_qcs}

where :math:`h_l` is the factor between them (i.e. the *mass* fraction
of the injected condensate that is liquid). It is not necessary in this
theory to have :math:`h_l` equal to :math:`g_l`: if a mixed phase plume
exists :math:`g_l + g_i` need not equal 1, but since :math:`h_l` and its
ice equivalent, :math:`h_i`, refer to mass, :math:`h_l + h_i` must equal
1. However, if we do not allow a mixed phase injection (which is the
case in the current mass-flux convection scheme, where only one phase
can be injected), :math:`h_l` and :math:`g_l` are equal (and either zero
or one in the current mass-flux convection scheme) and we can write
(`[eq:dctdt_final] <#eq:dctdt_final>`__) as

.. math::

   \frac{\partial{C_l}}{\partial{t}} =
   \frac{ (\delta_{xl} - C_l) }{ \delta_{xl} q_{c}^S - \overline{q_{cl}} }
   Q4_l 
   \label{eq:dctdt_xl}

where :math:`\delta_{xl} = h_l = g_l`. Equivalent expressions exist for
the ice cloud fraction and total cloud fraction.

.. math::

   \frac{\partial{C_i}}{\partial{t}} =
   \frac{ (\delta_{xi} - C_l) }{ \delta_{xi} q_{c}^S - \overline{q_{cf}} }
   Q4_i 
   \label{eq:dctdt_xi}

.. math::

   \frac{\partial{C_t}}{\partial{t}} =
   \frac{ (1 - C_t) }{ q_{c}^S - \overline{q_{c}} } Q4_c
   \label{eq:dctdt_xc}

with :math:`\delta_{xi} = h_i = g_i`. These are the expressions that are
used within the convection scheme. It still remains to parametrize
:math:`\delta_{xl}`, which is given by the convection scheme itself.
This is discussed in section `4.7.9 <#sec:plume_phase>`__.

.. _`sec:multi_numapp`:

Numerical application
~~~~~~~~~~~~~~~~~~~~~

The numerical application using
(`[eq:dcltdt_almost_final] <#eq:dcltdt_almost_final>`__) may be
performed with a basic forward timestep. Each of the three cloud
fractions can be incremented, assuming we know
:math:`\Delta{\overline{q_{cl}}}` and :math:`\Delta{\overline{q_{cf}}}`,
as

.. math::

   \Delta{C_t} = \frac{(1 - C_t)} {q_c^S - \overline{q_{cl}} - \overline{q_{cf}}}
   ( \Delta{\overline{q_{cl}}} + \Delta{\overline{q_{cf}}} ), 
   \label{eq:cft_ts}

.. math::

   \Delta{C_l} = \frac{ (g_l - C_l)}
   {q_c^S - \overline{q_{cl}} - \overline{q_{cf}}}
   ( \Delta{\overline{q_{cl}}} + \Delta{\overline{q_{cf}}} ),
   \label{eq:cfl_ts}

.. math::

   \Delta{C_i} = \frac{ (g_i - C_i)} 
   {q_c^S - \overline{q_{cl}} - \overline{q_{cf}}}
   ( \Delta{\overline{q_{cl}}} + \Delta{\overline{q_{cf}}} ).
   \label{eq:cff_ts}

The application from within the convection scheme is slightly different.
We start with (`[eq:dctdt_xl] <#eq:dctdt_xl>`__), but enforce two
numerical restrictions to avoid the equation set becoming
ill-conditioned. Firstly, we limit the denominator
:math:`q_c^S - \overline{q_{cl}}` to a minimum value if we are
considering changes of the same phase as the injected source.

.. math::

   \Delta C_l = \frac  {\delta_{xl} - C_l}  {\delta_{xl} \text{Max}( q_{c}^S
   - \overline{q_{cl}} , q_c^{S0} ) + ( 1 - \delta_{xl} ) (-\overline{q_{cl}}) }
   Q4_l
   \label{eq:delta_cl}

where :math:`q_c^{S0}` is specified as
:math:`5 \times 10^{-5} kg \, kg^{-1}`. The denominator also has an
additional check. If its absolute value is less than a tolerance value
of :math:`1 \times 10^{-10} kg \, kg^{-1}` then no change in cloud
fraction will be considered. A similar equation is used for the ice
cloud and the change in total cloud fraction

.. math::

   \Delta C_i = \frac  {\delta_{xi} - C_i}  {\delta_{xi} \text{Max}( q_{c}^S
   - \overline{q_{ci}} , q_c^{S0} ) + ( 1 - \delta_{xi} ) (-\overline{q_{cf}}) }
   Q4_i
   \label{eq:delta_ci}

.. math::

   \Delta C_t = \frac  {1 - C_t}
   { \text{Max}(q_c^S - \overline{q_c} , q_c^{S0} ) } Q4_c .
   \label{eq:delta_ct}

We now limit the change in cloud fraction to ensure that the cloud
fraction remains within its physical bounds.

.. math::

   C_l^{[n+1]} = ( 0,  C_l^{[n]} + \Delta C_l, 1)
   \label{eq:delta_cl_conv_final}

and similar equations are used for :math:`C_i^{[n+1]}` and
:math:`C_t^{[n+1]}`.

.. _`sec:conv_imp_note`:

A note on the implementation of the cloud fraction change
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Equation `[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__ has been derived
assuming that the only change in the cloud properties within the gridbox
comes from the detrainment of air from the convective plume (so that the
injection source is an appropriate model). Attention should be drawn to
the fact that this is not the only source of change from the convection
scheme. Two other terms require consideration, namely advection of the
environmental air downwards by compensating subsidence and the
condensation resulting from the adiabatic warming due to this
subsidence. The former is considered correctly in the calculation of
:math:`\frac{\partial \overline{q_{cl}}}{\partial t}`, which corresponds
to :math:`Q4`. However, the calculation of
:math:`\frac{\partial C_l}{\partial t}` is then performed using
(`[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__) and **incorrectly** assuming
that all the :math:`\overline{q_{cl}}` change comes from the
detrainment. It is possible to calculate directly the change in
:math:`C_l` that should occur due to the detrainment and compensating
subsidence treated together, in the same way that
:math:`\Delta \overline{q_{cl}}` is calculated (see section
`4.7.3 <#subsect:q4calculation>`__), and this is the way in which the
cloud fraction change **should** be done. It is an unfortunate
historical emphasis in the early development of PC2 on the derivation of
(`[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__) that has led to the treatment
used within the Unified Model for the change in cloud fractions due to
convection.

The change in :math:`\overline{q_{cl}}` and :math:`C_l` due to the
adiabatic warming associated with the compensating subsidence is
considered explicitly in the model implementation (see section
`4.7.4 <#sec:conv_homog>`__) for both :math:`\overline{q_{cl}}` and
:math:`C_l` after the rest of the convective process has been
calculated. It is perhaps arguable that if
(`[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__) is going to be applied then the
value of :math:`Q4` used in (`[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__)
should include this term.

Any major future developments of PC2 for a mass-flux convection scheme
would be advised to consider whether it is appropriate to use
(`[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__) at all.

.. _`sec:ct`:

Ice cloud and mixed phase regions
---------------------------------

The homogeneous forcing, initiation and PC2 erosion sections described
above have only considered the generation and dissipation of liquid
clouds. Although the forcing methods will not influence the generation
and dissipation of ice cloud (which is primarily performed in the
large-scale precipitation scheme, section `4.2 <#sec:precip>`__) we are
still left with the issue of how created or dissipated liquid cloud
overlaps with existing ice cloud in the gridbox. The opposite situation,
where changes in ice cloud are specified and changes in the overlap with
liquid cloud need to be calculated, is also possible in PC2 (e.g. in the
boundary layer, see section `4.6 <#sec:bl>`__).

Here we need a simple assumption to close the problem. The assumption
that we now choose is that liquid cloud fraction *changes* are
*minimally* overlapped with ice cloud fraction changes. This choice is
based upon observational evidence that mixed phase cloud is relatively
rare, and also on results from earlier PC2 development that indicated
less supercooled liquid water cloud than is observed from ground-based
lidar.

With this assumption, the equation set becomes straightforward to write
down. We firstly consider that a change in liquid cloud fraction
:math:`\Delta C_l` is known and we wish to estimate the resulting change
in the total cloud fraction. There is, of course, no change in the ice
cloud fraction :math:`C_i`, since, from our *definitions* in
(`[eq:dqcldt_and_dcdt] <#eq:dqcldt_and_dcdt>`__) and
(`[eq:mp] <#eq:mp>`__), this includes the mixed phase contribution.
Hence we write

.. math::

   \Delta C_i = 0 .
   \label{eq:deltaci_eq_0}

The change in the total cloud fraction, :math:`C_t` will depend upon the
sign of the change of the liquid cloud fraction. If
:math:`\Delta C_l > 0`, then :math:`\Delta C_t` is going to be the same
as :math:`\Delta C_l` (:math:`C_l` is being added with minimum overlap
to :math:`C_i`), unless the gridbox becomes completely covered in cloud,
when there is no choice but to generate mixed phase cloud. Hence we have

.. math::

   \Delta C_t = \text{Min} ( \Delta C_l , 1 - C_t ).
   \label{eq:deltact_min}

If :math:`\Delta C_l < 0`, then we still consider minimum overlap of the
*changes* (this is so that the solution is reversible as much as
possible). Hence :math:`\Delta C_t` is going to be the same as
:math:`\Delta C_l` unless :math:`C_l` is reduced below the existing
:math:`C_i`, in which case no more change to :math:`C_t` is possible.

.. math::

   \Delta C_t = \text{Max} ( \Delta C_l , C_i - C_t ) ,
   \label{eq:deltact_min2}

remembering that both quantities in the maximum expression in
(`[eq:deltact_min2] <#eq:deltact_min2>`__) have negative values.

We can write similar expressions if a known amount of ice cloud is added
or removed, and we need to calculate the effect on :math:`C_t`. Similar
to the results above we have:

.. math::

   \Delta C_l = 0 .
   \label{eq:deltacl_eq_0}

and

.. math::

   \Delta C_t = \left\{ \begin{array}{ll}
               \text{Max} ( \Delta C_i , C_l - C_t ),   &  \Delta C_i < 0  \\
               \text{Min} ( \Delta C_i , 1 - C_t ),     & \Delta C_i > 0 .
                \end{array} \right.
   \label{eq:deltact_min_array}

For completeness, we also present here the equation set for random
overlap of changes in liquid cloud with existing ice cloud. We have, as
before,

.. math::

   \Delta C_i = 0 .
   \label{eq:deltaci_eq_0_2}

For :math:`\Delta C_l > 0` additional liquid cloud is added randomly to
any location outside that of the current liquid cloud. A proportion
:math:`\frac{1-C_t}{1-C_l}` of this will be additionally outside that of
existing ice cloud. Hence the net change in total cloud fraction can be
written as

.. math::

   \Delta C_t = \Delta C_l \frac{1 - C_t}{1 - C_l} .
   \label{eq:deltact_ran1}

Similarly, if :math:`\Delta C_l < 0`, the liquid cloud is removed
randomly from the existing liquid cloud. A proportion
:math:`\frac{C_t - C_i}{C_l}` of this is from liquid cloud that does not
overlap with existing ice cloud. Hence,

.. math::

   \Delta C_t = \Delta C_l \frac{C_t - C_i}{C_l} .
   \label{eq:deltact_ran2}

Equivalent equations to (`[eq:deltact_ran1] <#eq:deltact_ran1>`__) and
(`[eq:deltact_ran2] <#eq:deltact_ran2>`__) but with :math:`C_l` and
:math:`C_i` swapped apply when we need to estimate changes in
:math:`C_t` from a known :math:`\Delta C_i`, when assuming random
overlap.

Numerical Implementation
~~~~~~~~~~~~~~~~~~~~~~~~

In general, although the situation does not occur within the current
implementation of PC2 , we might have increments to both :math:`C_l` and
:math:`C_i` simultaneously. Hence the implementation is to calculate
:math:`\Delta C_t` from the sum of that predicted by
(`[eq:deltact_min] <#eq:deltact_min>`__) or
(`[eq:deltact_min2] <#eq:deltact_min2>`__), and
(`[eq:deltact_min_array] <#eq:deltact_min_array>`__). For the random
overlap situation we also need to apply a check on the denominator in
(`[eq:deltact_ran1] <#eq:deltact_ran1>`__) and
(`[eq:deltact_ran2] <#eq:deltact_ran2>`__) before calculation, with the
result set to the limit :math:`\Delta C_t = 0` if the denominator is 0.
For the minimum overlap situation a final check is made that :math:`C_t`
lies between 0 and 1, with the value being reset to 0 or 1 if not.

Forced convective cloud
-----------------------

Forced convective clouds are clouds that form at the top of a convective
boundary layer but are too shallow to reach their level of free
convection (and become fully fledged cumulus clouds). These clouds
currently require special treatment because initiation in PC2 uses the
Smith scheme with a specified value of :math:`RH_{crit}` while the large
:math:`RH` variability associated with these clouds implies much lower
values than are typically used.

A profile of “forced cloud fraction”, :math:`C_{forced}`, is
parametrized as linearly varying with height between a cloud-base value,
at the lifting condensation level (LCL) from the convection diagnosis
parcel ascent, and a cloud-top value of 0.1 at the top of the capping
inversion, :math:`z_i^{top}`. The cloud-base value of :math:`C_{forced}`
varies linearly between 0.1 and 0.3 for cloud depths between 100 m and
300 m based loosely on SGP ARM site observations
:raw-latex:`\cite{zk13}`. The inversion top is taken to be the boundary
layer depth, :math:`z_h` plus the inversion thickness,
:math:`\Delta z_i` parametrized following :raw-latex:`\cite{rb08}` as:

.. math::

   \Delta z_i  = 6.3 \, w_m^2 /  \int_{z_h}^{z_h+\Delta z_i} b \, dz 
   \label{dz_param}

where :math:`w_m` is the boundary layer velocity scale
(:math:`w_m^3 = u_*^3 + 0.25 w_*^3`) and :math:`b` is the parcel
buoyancy that is integrated over the depth of the inversion assuming a
piece-wise linear variation between grid-levels. Note that the constant
in (`[dz_param] <#dz_param>`__) is the same as in
:raw-latex:`\cite{rb08}` because :math:`6.3 = 2.5 * 4^{2/3}` and
:math:`w_m^3` differs by a factor of 4.

The in-cloud water content at the top of the inversion is estimated
using the water content from the diagnostic parcel ascent (used to
diagnose boundary layer type and trigger convection), with linear
interpolation used between the lifting condensation level and inversion
top. To allow for sub-adiabatic water content (due to lateral mixing or
microphysical processes) the in-cloud water content can be reduced by a
factor, forced_cu_fac, that has been set to 0.5 in GA7.

These cloud fraction and water content profiles are then used as minimum
values and increments to :math:`C` and :math:`\overline{q_{cl}}`
calculated if necessary. This methodology can also optionally be applied
to cloud layers diagnosed as cumulus, if the boundary layer option to
mix across the lifting condensation level is selected that generates a
cloud base transition zone thickness which is then treated analgously to
the inversion thickness above.

Also, there is an option to treat the calculated forced cumulus cloud
fraction and water content as diagnostic quantities passed directly to
the radiation scheme as part of the “convective” cloud, instead of using
them to modify the prognostic “large-scale” cloud variables :math:`C`
and :math:`\overline{q_{cl}}`. If this option is used, the convective
cloud fraction :math:`CCA` and water content :math:`CCW` output by the
convection scheme are updated, by taking the forced cumulus profiles as
their minimum allowed values. Note that only the convective cloud fields
passed to radiation are updated (i.e. the versions of :math:`CCA` and
:math:`CCW` that are stored in the model dump / D1 array). The UM code
contains other copies of the convective cloud fields that are only used
for diagnostics; these are *not* updated.

The different options for how to treat forced cumulus cloud are
controlled by the cloud namelist input :math:`forced\_cu`, and are
summarised below:

- :math:`forced\_cu = 0`: No treatment of forced cumulus clouds.

- :math:`forced\_cu = 1`: Forced cumulus cloud applied to :math:`C` and
  :math:`\overline{q_{cl}}` only in dry-convective boundary-layers.

- :math:`forced\_cu = 2`: Forced cumulus cloud applied to :math:`C` and
  :math:`\overline{q_{cl}}` in both dry-convective and cumulus-capped
  boundary-layers.

- :math:`forced\_cu = 3`: Forced cumulus cloud applied to :math:`CCA`
  and :math:`CCW` in both dry-convective and cumulus-capped
  boundary-layers.

.. _`sec:turb_qcl_scheme`:

Turbulence-driven production of subgrid scale liquid cloud
----------------------------------------------------------

.. _`sec:sgt_intro`:

Introduction
~~~~~~~~~~~~

:raw-latex:`\cite{fhfk14}` developed a model for subgrid liquid water
production by turbulent motions. Their method uses an exactly soluble
stochastic process to describe subgrid relative humidity (RH)
fluctuations. The probability density function (PDF) of the fluctuations
can be diagnosed in terms of the local turbulent local state and any
pre-existing ice cloud. The liquid cloud properties (cloud fraction and
liquid water content) can be then be calculated as truncated moments of
the PDF.

:raw-latex:`\cite{fhfk14}` initially used their model to understand and
parametrize the results of Large Eddy Simulations (LES) of
shear-induced, Altostratus clouds. They obtained excellent agreement
between their theoretically predicted predicted mean cloud properties
and the bulk properties of the LES clouds. Subsequently, their model has
been used as the basis of subgrid cloud initiation method for use in the
Unified Model in conjunction with the PC2 prognostic cloud scheme. In
Section `3.8.2 <#sec:sgt_model_describe>`__ we outline the model of
:raw-latex:`\cite{fhfk14}`. In Section
`3.8.3 <#sec:sgt_model_implement>`__ we described its implementation in
the GCM.

.. _`sec:sgt_model_describe`:

Model description
~~~~~~~~~~~~~~~~~

:raw-latex:`\cite{fhfk14}` started from the equation for the dynamics of
ice supersaturation :math:`S_i=e_v/e_{sat\;ice}-1`:

.. math::

   \label{eqn:squires_eqn}
    \frac{D S_i}{D t} = -b_i B_0 {\cal M}_1 S_i
      -\left(\frac{\varepsilon}{L^2}\right)^{1/3}(S_i-S_E) + a_i w,

where :math:`{\cal M}_1` is the first moment of ice particle size
distribution (PSD), :math:`\varepsilon` is the turbulent dissipation
rate, :math:`L` is a prescribed mixing length for the turbulence,
:math:`S_{\rm E}` is the ice supersaturation of the environment
surrounding the cloud and :math:`b_i,B_0` and :math:`a_i` are function
of :math:`p` and :math:`T` given by

.. math::

   \begin{aligned}
    b_i &=& \frac{1}{q} + \frac{\epsilon L_s^2}{c_p R T^2}, \\
    B_0 &=& 4\pi C \left( \frac{\epsilon L_s^2}{K_a R T^2} + \frac{R T}{\epsilon e_{si} \psi} \right)^{-1}, \\
    a_i &=& \frac{g}{R T}\left( \frac{\epsilon L_s}{c_p T} - 1 \right), \\
   \end{aligned}

The first term on the right hand side of
Eq. `[eqn:squires_eqn] <#eqn:squires_eqn>`__ is the sink of vapor due to
depositional growth of ice crystals, the second term models entrainment
(mixing) of environmental air into the cloudy volume and the third term
is a source term due to vertical air motions.

:raw-latex:`\cite{fhfk14}` modeled vertical velocity as a white-noise
process with autocorrelation function:

.. math:: \overline{w(t)w(s)} = \sigma_w^2 \tau_{\rm d} \delta(t-s),

where :math:`\delta` is the Dirac distribution and the intensity of the
noise, :math:`\sigma_w^2`, will be called the standard derivation of the
vertical velocity fluctuations (due to the white nature of noise, a true
expectation value :math:`\overline{w^2}` is not defined) and
:math:`\tau_{\rm d}` a Lagrangian decorrelation time define here by the
relation used by :raw-latex:`\cite{rodean1997}`:

.. math::

   \tau_{\rm d} = \frac{2\sigma_w^2}{\varepsilon C_0},
   \label{eqn:taud}

where :math:`C_0` is a known constant.

Because it is linear in :math:`S_i`, Equation
`[eqn:squires_eqn] <#eqn:squires_eqn>`__ can be solved exactly, for any
given realisation of the noise term. By averaging the solutions over the
the noise and taking a steady-state limit (see
:raw-latex:`\cite{fhfk14}` for details) it can be shown that the
solution PDF is Gaussian with mean and variance given by:

.. math::

   \begin{aligned}
    \overline{S_i} &=& 
     S_{\rm E}\frac{ \left(\varepsilon/L^2\right)^{1/3}  }{ b_i B_0 {\cal M}_1  + \left(\varepsilon/L^2\right)^{1/3} }.
     \label{eqn:si_avg} \\
    \overline{S_i^2} &=&
     \frac{a^2_{\rm i} \sigma^2_w \tau_{\rm d}}{ 2\left(b_i B_0 {\cal M}_1  + \left(\varepsilon/L^2\right)^{1/3}\right)},
    \label{eqn:si_var}
   \end{aligned}

Equation `[eqn:si_avg] <#eqn:si_avg>`__ and
`[eqn:si_var] <#eqn:si_var>`__ completely specify the PDF,
:math:`F(S_i)`, of steady-state humidity variations for the subgrid
model. The liquid cloud fraction and liquid water mass mixing ratio are
given by

.. math::

   \begin{aligned}
    C_l^{sgt} &=& \int_{S_{i,wat}}^\infty d S_i F(S_i), \label{eqn:cloud_fraction} \\
    q_{cl}^{sgt} &=&  q_{sat\;ice}\int_{S_{i,wat}}^\infty d S_i (S_i -S_{i,wat}) F(S_i) \label{eqn:cloud_liquid},
   \end{aligned}

where :math:`S_{i,wat} = e_{sat\;wat}/e_{sat\;ice}-1` is the value of
ice supersaturation at water saturation. We use the superscription
‘:math:`sgt`’(=‘*s*\ ub\ *g*\ rid *t*\ urbulence’) to indicate that
:math:`C_l^{sgt}` and :math:`q_{cl}^{sgt}` are values of cloud fraction
and water content diagnosed from a parametrization of small-scale
turbulent processes.

.. _`sec:sgt_model_implement`:

Model implementation and closure relations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To implement the model of Section `3.8.2 <#sec:sgt_model_describe>`__ in
the Unified Model, closure relations are needed for the quantities
:math:`\sigma_w^2`, :math:`\varepsilon`, :math:`L`, :math:`\tau_{\rm d}`
and :math:`S_E`, subject to the constraining relationship given by Eq.
`[eqn:taud] <#eqn:taud>`__. In each model grid box, these parameters
specify the subgrid PDF, :math:`F(S_i)`, and from this the liquid cloud
fraction and water content produced by turbulence can be found using Eqs
`[eqn:cloud_fraction] <#eqn:cloud_fraction>`__ and
`[eqn:cloud_liquid] <#eqn:cloud_liquid>`__.

In addition we need to make some assumptions about how the diagnosed
values :math:`C_l^{sgt}` and :math:`q_{cl}^{sgt}` relate to the model
prognostic fields, :math:`C_l` and :math:`q_{cl}`. Two methods are
available for doing this. In the simplest case, the diagnosed values
:math:`C_l^{sgt}` and :math:`q_{cl}^{sgt}` are just treated as
increments to model prognostics (option one, in Sec.
`3.8.5 <#sec:sgt_increments>`__ below). A more complex option (see
option two, below) is to increment the model fields via the PC2 Erosion
functionality.

.. _`sec:sgt_closures`:

Closure relations
~~~~~~~~~~~~~~~~~

The vertical velocity variance, :math:`\sigma_w^2`, is available as a
diagnostic from the Boundary Layer scheme. Because the Boundary Layer
scheme is called after the Microphysics on each model timestep, the
diagnostic value is stored in a (non-advected) model prognostic field.
The scheme will operate only where there is diagnosed turbulence, i.e.,
non-zero :math:`\sigma_w^2`.

We take the mixing length scale, :math:`L`, to be proportional to the
vertical grid spacing in each grid box: :math:`L=\beta_{mix} \Delta z`,
where :math:`\Delta z` is calculated as the height different between the
:math:`\rho`-levels adjacent to the given :math:`\theta`-point. The
parameter, :math:`\beta_{mix}`, is an adjustable constant that the user
can define (see Section `3.8.6 <#sec:sgt_options>`__ below), however it
should be of order one.

To obtain :math:`\tau_{\rm d}` we impose an eddy size constraint:

.. math::

   \tau_{\rm d} = \frac{L}{\sigma_w} = \beta_{mix} \frac{\Delta z}{\sigma_w}
   \label{eqn:eddy_size}

Eq. `[eqn:taud] <#eqn:taud>`__ then determines the dissipation rate,
:math:`\varepsilon`, that is consistent with the other parameters. The
constant :math:`C_0=10` by default, but can be adjusted by the user.

The scheme is limited to act only in grid boxes where
:math:`\tau_{\rm d}` is less than a prescribed value,
:math:`\tau_{d}^{max}`. The default is
:math:`\tau_d^{max}=1200\;{\rm sec}`, which typically coincides with a
couple of model timesteps. The motivation for this is that a motion that
takes longer than a few timestep to decorrelate will be partially
resolved by the dynamics and therefore cannot be considered as ‘subgrid’
turbulence.

Finally, where :math:`T`, :math:`p` and :math:`q` appear in the
expressions for :math:`C_l^{sgt}` and :math:`q_{cl}^{sgt}`, these are
taken to be the grid box mean values. The first moment of the ice PSD,
:math:`{\cal M}_1`, is found from the parametrization, due to
:raw-latex:`\cite{fhbicc05}`, described in Section 4.1 of UMDP26.

.. _`sec:sgt_increments`:

Options for incrementing model prognostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the information in Section `3.8.4 <#sec:sgt_closures>`__ to obtain
closed expressions for the subgrid PDF of :math:`S_i`-fluctuations
allows :math:`C_l^{sgt}` and :math:`q_{cl}^{sgt}` to be calculated.
These will be non-zero only where there is turbulence as diagnosed by
the Boundary Layer scheme (and hence non-zero :math:`\sigma_w^2`). To
calculate :math:`C_l^{sgt}` and :math:`q_{cl}^{sgt}` the integrals in
Eqs `[eqn:cloud_fraction] <#eqn:cloud_fraction>`__ and
`[eqn:cloud_liquid] <#eqn:cloud_liquid>`__ are evaluated numerically
using discretisation based on user-specified number of bins.

Given :math:`C_l^{sgt}` and :math:`q_{cl}^{sgt}`, two options are
available for relating these to changes in the model prognostics:

Option one: direct increments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The values of :math:`C_l^{sgt}` and :math:`q_{cl}^{sgt}` can be added as
increments to the model prognostic fields, :math:`C_l` and
:math:`q_{cl}`. In this case

.. math::

   \begin{aligned}
    \left( \Delta C_l \right)_{sgt} &=&  C_l^{sgt} \\ 
    \left( \Delta q_{cl} \right)_{sgt} &=& q_{cl}^{sgt}, \\
    \left( \Delta q \right)_{sgt} &=& -\left( \Delta q_{cl} \right)_{sgt}, \\
    \left( \Delta T \right)_{sgt} &=& \frac{L_v}{c_p} \left( \Delta q_{cl} \right)_{sgt}, \\
    \left( \Delta C \right)_{sgt} &=&  C_l^{sgt} \\
   \end{aligned}

where the left hand sides denote the increments to :math:`C_l`,
:math:`q_{cl}`, :math:`T` and the total cloud fraction, :math:`C`, due
to the subgrid scheme. Some bounds-checking is then applied to ensure
that: (a) the resultant cloud fractions to not exceed one; (b) the
scheme does not condense out more liquid than there is available
moisture.

Option two: PC2 Erosion method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Option one gives a simple method for incrementing the model prognostics,
but it gives rise to a potential inconsistency with the PC2 cloud
scheme. This arises because the subgrid production scheme can elevate
cloud fraction to unity in grid boxes that are subsequently diagnosed by
PC2 Initiation to meet the criteria for clear-sky initiation. PC2 then
counteracts the scheme by removing some of the liquid cloud. To try to
mitigate against this issue, cloud fraction increments can be applied
using PC2 Erosion. In this case:

.. math::

   \begin{aligned}
    \left( \Delta q_{cl} \right)_{sgt} &=& q_{cl}^{sgt} - q_{cl}, \\
    \left( \Delta q \right)_{sgt} &=& -\left( \Delta q_{cl} \right)_{sgt}, \\
    \left( \Delta T \right)_{sgt} &=& \frac{L_v}{c_p} \left( \Delta q_{cl} \right)_{sgt}, \\
   \end{aligned}

where :math:`q_{cl}` is the liquid cloud amount prior to calling to the
turbulent production scheme. The cloud fraction increments are
calculated by calling PC2 Erosion with
:math:`\left( \Delta q_{cl} \right)_{sgt}` as input. See Section
`4.3 <#sec:turb>`__ for details on how the PC2 Erosion process works.
This method gives cloud fraction increments that are consistent with PC2
cloud scheme.

.. _`sec:sgt_options`:

Other user options
~~~~~~~~~~~~~~~~~~

The following variables and logical switches are optional inputs:

#. The logical ``l_dcfl_by_erosion`` provides a switch to apply cloud
   fraction increments using PC2 Erosion. Defaults to *FALSE*.

#. Setting the logical ``l_mixed_phase_t_limit`` to *TRUE* allows the
   user to use the variable ``mp_t_limit`` to define a temperature
   limit, :math:`T_{max}`, above which the scheme is not applied. The
   default is :math:`T_{max}=0^\circ\;{\rm C}`, so the scheme is only
   applied to cold clouds.

#. The input variable ``mp_tau_d_lim`` defines the upper limit,
   :math:`\tau_d^{max}`, on the value of :math:`\tau_d` above which the
   scheme is not applied. The default value is
   :math:`\tau_d^{max}=1200.0`, so the scheme is not applied in grid
   boxes where the decorrelation time scale exceeds :math:`1200`
   seconds.

#. ``nbins_mp`` is the number of bins used in the discretisation of the
   integrals in Eqs `[eqn:cloud_fraction] <#eqn:cloud_fraction>`__ and
   `[eqn:cloud_liquid] <#eqn:cloud_liquid>`__ for :math:`C_l^{sgt}` and
   :math:`q_{cl}^{sgt}`. The default value is :math:`100` bins.

#. ``mp_dz_scal`` is the scale parameter, :math:`\beta_{mix}`, in the
   definition of the mixing length, :math:`L=\beta_{mix}\Delta z`.

#. ``mp_czero`` defines the constant parameter :math:`C_0` (defaults to
   :math:`C_0=10`).

.. _`sec:app_um`:

Application to the Unified Model
================================

This section describes the way in which the physical concepts described
in the above section are applied to the sections of the Unified Model,
in order to build up the complete prognostic scheme. Description of the
actual subroutines themselves follow in section `5.2 <#sec:code>`__.
Note that the large-scale precipitation and convection schemes have
considerable documentation below, since these schemes have been heavily
modified for PC2. The other schemes use generic forcing scenarios, hence
their desciption here is much shorter. Remember, whenever a signficiant
:math:`\overline{T}` or :math:`\overline{q}` change occurs, PC2 must be
able to represent the corresponding condensation and changes in cloud
fractions.

.. _`sec:rad`:

Radiation
---------

The shortwave and longwave radiation schemes both alter the temperature
of the atmosphere, hence we need to calculate the corresponding
condensation and cloud fraction changes. For both shortwave and
longwave, we use the homogeneous forcing routines (section
`3.2 <#sec:homog>`__) for :math:`\overline{q_{cl}}` and :math:`C_l`,
(using eqn. `[eq:deltaqc_exp2] <#eq:deltaqc_exp2>`__ to calculate the
:math:`Q_c` forcing) and then the method in section `3.6 <#sec:ct>`__ to
calculate :math:`C_t` changes. There is no :math:`\overline{q_{cf}}`
change associated with this process since the deposition / sublimation
process is performed within the large-scale precipitation scheme (as it
also is in the absence of PC2).

It is reasonable to question whether homogeneous forcing is a reasonable
model to use when we know that a large proportion of the heating
associated with radiative transfer in the atmosphere comes from the
cloudy air and is not evenly spread across the gridbox. Possible
developments are discussed in section `7.4.7 <#sec:homog_improve>`__.

.. _`sec:precip`:

Large-scale precipitation
-------------------------

Precipitation processes have a large effect on cloud fractions. Here we
present the simple physical models that are applied to the transfer
terms included in the large-scale precipitation scheme. They are also
presented within the large-scale precipitation documentation ().

The basis of the physical model is that microphysical transfer processes
can be calculated separately in different partitions of the model cloud
(i.e. mixed phase cloud, liquid phase cloud, ice phase cloud or clear
sky). However, processes may change the size of these partitions. We
consider here separately each process that is modelled in the
large-scale precipitation scheme. The changes in
:math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}` and
:math:`\overline{q}` remain mathematically the same as in the non-PC2
version of the code (), we only need to introduce calculations for the
changes in cloud fractions. We will see that many of these changes can
be well modelled by assuming no change to the cloud fractions, and the
others by using simple assumptions.

Although the model may use two ice prognostic ice categories, only a
single ice cloud fraction is stored, the assumption being that the two
ice categories are completely overlapped with each other. Graupel is not
considered to contribute to the ice cloud fraction.

.. _`sec:lsp_fall`:

Fall of ice
~~~~~~~~~~~

The fall of ice is the process that contributes most to the growth of
ice cloud fraction in the model. The model results are therefore
sensitive to the formulation of this process. We will make the basic
assumption that a trail of falling ice does not reduce the horizontal
spread of ice cloud fraction at a particular level (hence
:math:`\overline{q_{cf}}` that leaves a gridbox does not reduce
:math:`C_f` in that gridbox). The in-cloud ice content simply reduces
due to the fall out of ice - it is the sublimation term (section
`4.2.4 <#sec:mp_depsub>`__) that erodes the fall streaks. However, ice
that falls into a clear layer from above may increase the ice cloud
fraction. We parametrize this by considering the horizontal overlap of
ice clouds between two model layers, and the fall speed of ice between
them. We will assume an overlap that is nearly, but not quite, maximum,
the difference being dependent upon the windshear and the time taken for
ice to fall between the levels.

.. math::

   O^{[k,k+1]} = \text{Max}( C_{i}^{[k+1]} - C_i^{[k]}  , 0)
   + w \frac{\Delta z^{[k]}}{v_i^{[k]}}
   \label{eq:overhang}

where :math:`O^{[k,k+1]}` is the amount of ice cloud ‘overhanging’ the
current (i.e. :math:`k`\ ’th) layer from the layer above, :math:`w` is a
parameter that is closely related to the windshear,
:math:`\Delta z^{[k]}` is the model layer thickness and
:math:`v_i^{[k]}` is the fallspeed of ice in the layer.
:math:`v_i^{[k]}` is calculated in the microphysics scheme and, if two
ice prognostics are used, is the mass-weighted average fall speed of the
two categories. The factor :math:`\frac{\Delta z}{v_i}` is simply the
time taken for the ice to fall through one model layer. Multiplying this
by the windshear would give an estimate to the amount of overlap between
a cloud source and its fall streak in the layer below (it is an
*estimate* since we assume that the cloud source is continuous and
unbroken). Although it is quite possible within PC2 to do this, to date
we have not programmed this link, and we use a constant, but tunable,
value of :math:`1.5 \times 10^{-4} s^{-1}` for :math:`w`.

The change in :math:`C_i` over the timestep is then given by the overlap
proportion multiplied by the how much (in the vertical dimension) of the
layer below can be filled by ice in the timestep:

.. math::

   \Delta C_i = \text{Max}(O^{[k,k+1]} , 1)  \text{Min} (v_i \frac{\Delta t}{\Delta z^{[k]}} , 1)
   \label{eq:lsp_fall}

where :math:`\Delta t` is the timestep. We now choose to assume a
minimum overlap between the liquid and the ice phases (as in section
`3.6 <#sec:ct>`__).

.. math::

   \Delta C_t = \text{Min} ( \Delta C_i , A_{clear} )
   \label{eq:lsp_fall_ct}

where :math:`A_{clear}` is the proportion of the gridbox that has
neither ice nor liquid cloud present.

**An inconsistency has been found in the way that the fall-of-ice term
is linked to the globally constant “wind-shear value” when calculting
the ice cloud fraction overhang. Consequently, the option not to use the
“wind shear value” when calculating the overhang is available in the
UMUI (from version 7.6 onwards).**

.. _`sec:lsp_homo`:

Homogeneous nucleation
~~~~~~~~~~~~~~~~~~~~~~

This will freeze all supercooled liquid water when a temperature
threshold is exceeded. Hence we turn all existing liquid and mixed phase
cloud to ice cloud. The cloud fraction changes are:

.. math::

   \begin{aligned}
   C_l \leftarrow 0  \nonumber \\ 
   C_i \leftarrow C_t \nonumber \\
   \Delta C_t = 0.
   \label{eq:lsp_homo}
   \end{aligned}

Heterogeneous nucleation
~~~~~~~~~~~~~~~~~~~~~~~~

This process will freeze a small amount of supercooled liquid water,
regardless of the previous presence of ice cloud. This will mean that
previously existing ‘liquid-only’ cloud is converted to mixed phase
cloud. These give the following changes:

.. math::

   \begin{aligned}
   \Delta C_l = 0  \nonumber \\
   C_i \leftarrow C_t \nonumber \\
   \Delta C_t = 0.
   \label{eq:lsp_het}
   \end{aligned}

.. _`sec:mp_depsub`:

Deposition and sublimation
~~~~~~~~~~~~~~~~~~~~~~~~~~

This term exerts one of the most important influences on the ice cloud
in the whole model (this applies to the control as well as for PC2).
Contained in the formulation is a subgrid-scale assumption that causes
equivalent effects to that for a moisture PDF under the ‘:math:`s`’
framework (section `2.2 <#sec:s_dist>`__). However, since
:math:`{q_{cf}}` changes slowly in response to local changes in
:math:`q` and :math:`T`, we cannot base the :math:`q_{cf}` response on
the same instantaneous condensation framework. It would be useful to
investigate in the future whether the two descriptions of the moisture
variability could be brought together. Because of its importance, we
describe the method below, although we note it is also described in .

We can calculate the local rate of change of :math:`q_{cf}`, given local
:math:`T` and :math:`q` etc. using the standard microphysical growth
equations (see ). However, it is critical to know the way in which the
moisture is correlated with the ice in the gridbox. We will assume there
exists a distribution of vapour in the gridbox. We know that the regions
where liquid cloud exists must be saturated with respect to liquid
water, hence we need only consider the part of the gridbox that does not
have liquid water present. The average value, :math:`q_a`, of :math:`q`
within the liquid-free part of the gridbox is thus

.. math::

   q_a = \frac{  \overline{q} - C_l q_{sat \, liq}(\overline{T}) } {1 - C_l}
   \label{eq:qa}

where we have assumed that the fluctuation of :math:`q_{sat~liq}` across
the gridbox due to temperature fluctuations is not significant compared
to the fluctuation of :math:`q` described below. We then parametrize a
width, :math:`b_i`, to the :math:`q` (not :math:`s`) fluctuations
*across the non-liquid cloud part of the gridbox*, based upon
:math:`RH_{crit}`. This is like that for the ‘:math:`s`’ distribution
width, :math:`b_s` but modified:

.. math::

   b_i = (1 - RH_{crit} ) q_{sat \, liq} ( 1 - \frac{1}{2}
   ~ \frac{\overline{q_{cf}}} {i q_{sat \, liq}(\overline{T})} ) .
   \label{eq:b_i}

where the factor :math:`( 1 - \frac{1}{2}
\frac{\overline{q_{cf}}} {i ~ q_{sat \, liq}(\overline{T})})` should be
limited to a minimum value of zero, but, for numerical reasons, is
limited to a minimum value of 0.001. We note that :math:`b_i` has a
similar form to :math:`b_s`, except the multiplier :math:`a_L` and the
factor in brackets. If we remember from (`[eq:s] <#eq:s>`__) that the
definition of ‘:math:`s`’ includes a factor :math:`a_L` we see that the
absence of the :math:`a_L` factor in (`[eq:b_i] <#eq:b_i>`__) is
consistent. The factor in brackets is a *parametrization* of the effect
that, when ice is present, deposition in the moistier parts and
sublimation in the drier parts of the gridbox must reduce the width of
the distribution of :math:`q` across the gridbox. It is a simple linear
function of :math:`\frac{\overline{q_{cf}}}{q_{sat~liq}(\overline{T})}`,
and is tunable using the factor :math:`i`, which takes the value of
0.04.

We note that this formulation isn’t totally consistent with the liquid
cloud formulation, which considers an underlying PDF across the whole
gridbox and does not have, in general, its width prescribed. Remember
that we do not calculate on-line the whole of the liquid - vapour PDF,
we only parametrize the single point :math:`G(-Qc)`, using equation
`[eqn22] <#eqn22>`__).

The width is then limited further to be no greater than
:math:`\overline{q}`, to make sure that there are no negative values of
:math:`q` predicted within the gridbox (possible at low temperatures
where :math:`q_{sat~liq}(\overline{T})` diverges from
:math:`q_{sat~ice}(\overline{T})`).

We then calculate the average value of :math:`q` in the ice-only and
clear-sky partitions of the gridbox. To do this, we make the further
assumption that the ice is correlated with the moistest part of the
distribution (an instantaneous condensation formulation would make the
same assumption). Some algebra retrieves the expressions:

.. math::

   \begin{aligned}
   q_{clear} = q_a - b_i A_{ice} ; \\
   q_{ice} = \frac {\overline{q} - C_l q_{sat~liq} - A_{clear} q_{clear} }
   {A_{ice}},
   \label{eq:q_clear_and_q_ice}
   \end{aligned}

where :math:`A_{ice}` is the proportion of the gridbox with ice cloud
but not liquid cloud and :math:`A_{clear}` is the proportion of the
gridbox without cloud. The numerical application will set
:math:`q_{clear}` to :math:`q_a` if :math:`A_{ice}` is zero. We now have
a representation of the :math:`q` values in each of the gridbox cloud
partitions, and can solve the microphysical transfer equation in each
partition.

The cloud fraction changes now need to be parametrized. We use the model
that deposition will *not* adjust the ice cloud *fraction* (increases
will be done within the fall-of-ice microphysics section). However,
deposition can decrease the liquid cloud fraction (locally, :math:`q`
can be reduced by deposition to below :math:`q_{sat~liq}`, hence this is
not inconsistent with the assumptions for the riming term below. This is
the principal sink of supercooled liquid cloud fraction in the model.
Sublimation will be allowed to decrease the ice cloud fraction (since
sublimation cannot act in liquid cloud, there is no impact on the liquid
cloud). To solve for these models, we will need to further split the
ice-only partition to give the proportion of that partition that is
above and below ice saturation. This gives, in general, an area of the
gridbox :math:`A_{ice1}` that contains ice and is above saturation where

.. math::

   A_{ice1} = \frac{1}{2} A_{ice} + \frac{1}{2} 
    \frac{ (q_{ice}-q_{sat~ice}(\overline{T})) }  {b_i},
   \label{eq:q_ice_above_sat}

having assumed that :math:`A_{ice1}` is between 0 and :math:`A_{ice}`.
If not, it is trivial to partition the gridbox, since the moisture in
the ice-only partition is either completely above or completely below
:math:`q_{sat~ice}(\overline{T})`. The corresponding area that contains
ice and is below saturation is given by
:math:`A_{ice2} = A_{ice} - A_{ice1}`. We can now parametrize the change
in cloud fractions. For deposition, we shall assume a uniform
distribution of local values of :math:`q_{cl}` about the local mean. If
we assume a uniform removal of local :math:`q_{cl}` then, with a little
algebra, we can obtain an expression for the change in :math:`C_l`:

.. math::

   \Delta C_l = C_l ( 1 - \frac {\Delta \overline{q_{cl}}} {\overline{q_{cl}}} )
   ^{\frac{1}{2}} - C_l 
   \label{eq:deltacfl_dep}

.

Since this occurs only in the mixed phase part of the gridbox, we can
say that :math:`\Delta C_t = 0`. We will also note that the change in
:math:`\overline{q_{cl}}` due to deposition is limited by the amount of
:math:`\overline{q_{cl}}` that is in the mixed phase partition in the
gridbox, hence (`[eq:deltacfl_dep] <#eq:deltacfl_dep>`__), although it
formally allows removal of :math:`C_l` from an ice-free partition, will
be unlikely to do so.

The sublimation forms the main method by which ice cloud is destroyed in
PC2, hence PC2 results are relatively sensitive to its formulation. Here
we make a similar assumption to that used for liquid in the deposition
term, except that we limit the changes only to the region of the gridbox
where ice is subliming.

.. math::

   \Delta C_i = A_{ice2} ( 1 + \frac{\Delta \overline{q_{cf}} }
   {  \overline{q_{cf}} (  \frac{A_{ice2}}{C_i}   )     }
    )^{\frac{1}{2}} - A_{ice2}
   \label{eq:deltacfi_sub}

.

The term :math:`\overline{q_{cf}} (  \frac{A_{ice2}}{C_i}   )` is the
amount of :math:`\overline{q_{cf}}` that is present in the subliming ice
region, hence its ratio with :math:`\Delta \overline{q_{cf}}` is the
fractional change in that region. The change in the total cloud fraction
must also be equal to the change above, since sublimation cannot occur
in the presence of liquid cloud:

.. math::

   \Delta C_t = \Delta C_i .
   \label{eq:deltacft_sub}

Riming
~~~~~~

This process acts only where mixed phase cloud occurs - although, in
theory, it could remove any supercooled liquid totally, the air would
remain saturated with respect to liquid water. Hence any subsequent
cooling would regenerate the same amount of liquid cloud. Hence we
choose to model this process as having *no effect* on the cloud
fractions.

Capture
~~~~~~~

This is the freezing of raindrops onto ice crystals by collision. This
does not alter the ice cloud *fraction* in the gridbox (although it does
alter :math:`\overline{q_{cf}}`, and it has no interaction with the
liquid cloud. Again, we therefore choose to model this process as having
*no effect* on the cloud fractions.

Evaporation of melting ice
~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we simply assume that ice cloud fraction is removed in proportion
to the ice content that is removed.

.. math::

   \Delta C_i = C_i \frac{ \Delta \overline{q_{cf}}}{\overline{q_{cf}}} .
   \label{eq:lsp_evapmeltsnow}

Because the evaporation cannot occur in the liquid part of the gridbox,
there is no change to :math:`C_t` (or to :math:`C_l`).

Melting
~~~~~~~

Again, the change in :math:`C_i` is calculated using the method in
(`[eq:lsp_evapmeltsnow] <#eq:lsp_evapmeltsnow>`__).

.. math::

   \Delta C_i = C_i \frac{ \Delta \overline{q_{cf}}}{\overline{q_{cf}}} .
   \label{eq:lsp_melt}

The change in :math:`C_t` is calculated assuming that there is no
correlation in the gridbox between where the ice melts and the liquid
cloud. Hence we must multiply (`[eq:lsp_melt] <#eq:lsp_melt>`__) by the
proportion of ice cloud fraction that exists without liquid cloud (i.e.
:math:`\frac{A_{ice}}{C_i}`).

.. math::

   \Delta C_t = C_i \frac{ \Delta \overline{q_{cf}}}{\overline{q_{cf}}}
   \frac{A_{ice}}{C_i} .
   \label{eq:lsp_melt2}

Evaporation of rain
~~~~~~~~~~~~~~~~~~~

Evaporation of rain will not, *on its own*, generate liquid cloud, since
a large-scale lifting process will be required in order to condense
water from the moistened air. We cannot, therefore, allow any change in
cloud fractions to occur as a result, subsequent changes are calculated
elsewhere in the model (e.g. by the lifting process, section
`4.8 <#sec:pres>`__).

Accretion
~~~~~~~~~

Accretion is the sweep-out of liquid water droplets by rain. We argue in
a similar way to the riming term, that this will not remove any liquid
cloud fraction, since a small amount of lifting will regenerate the same
amount of liquid cloud. Hence we choose to model this process as having
*no effect* on the cloud fractions. The arguments underlying the
formulation of the evaporation of rain and the accretion cloud fraction
changes may appear to be inconsistent in their limiting cases and the
subsequent response to lifting. However, when the limiting case is not
reached the formulations are both correct. For the moment, it is not
considered necessary to increase the complexity of the current, simple
representations.

Autoconversion
~~~~~~~~~~~~~~

As for accretion, the generation of rain directly from collision and
coalescence of liquid water droplets will not alter the cloud fractions.

Other microphysics terms
~~~~~~~~~~~~~~~~~~~~~~~~

There are already (i.e. also in the control) two numerical tidy-up terms
at the end of the microphysics section that remove small rain amounts
and provide an additional melting term for the snow. These do not change
the cloud fractions.

If there are small amounts of ice present at the end of the microphysics
then these are removed at the end of the microphysics timestep (also in
the control). PC2 responds by resetting the cloud fractions
appropriately, so :math:`C_t` is reset to :math:`C_l` etc.

.. _numerical-implementation-1:

Numerical implementation
~~~~~~~~~~~~~~~~~~~~~~~~

Note that after each process has been applied, we do *not* recalculate
the sizes of the ice-only, liquid-only and mixed phase partitions, but
use the values at the start of the microphysics (this includes the
values of :math:`C_i` used in the calculation of ‘in-cloud’ water
contents above. However, we do update the cloud fractions themselves
sequentially. We also recalculate after each process the overlaps
between the rain fraction (see ) and the cloud fractions.

There is also a final set of checks that :math:`C_l` and :math:`C_i` lie
between 0 and 1 and that :math:`C_t` is bounded between
:math:`\text{Max}(C_l, C_i)` (maximum overlap of liquid and ice) and
:math:`\text{Min}(C_l+C_i,1)` (minimum overlap of liquid and ice).

We should note in particular, that these parametrizations allow a
considerable reduction in :math:`\overline{q_{cl}}` without a
corresponding large reduction in :math:`C_l`. This is an underlying
feature of the PC2 scheme (discussed in :raw-latex:`\cite{wg03}`), and
necessarily implies the skewing of the underlying moisture PDF.
Subsequent parts of the model (e.g. the width narrowing, section
`3.3 <#sec:width>`__) will, of course, act on the modified fields to
adjust the cloud fractions further, but remember that these are separate
processes and modelled elsewhere in the timestep.

.. _`sec:turb`:

PC2 erosion
-----------

Original width-narrowing method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(selected by setting **i_pc2_erosion_method = 1** in the UM namelist).

In parallel with the homogeneous forcing part of the PC2 response to
convection, we introduce a new block of code that allows a background
change of the PDF width. At earlier versions of PC2 (PC2:65 and earlier)
this block was included as a separate section of code that was called in
as part of the atmphya parallel timestepping. This was later moved to
better numerically balance increments from the convection scheme with
the cloud fraction erosion term. From VN8.1 onwards, a further option
was introduced to implement the erosion prior to the microphysics
parametrization. This was primarily to allow PC2 to be run at
convection-resolving scales, at which the convection scheme is not
called and therefore the erosion is not called. We have empirically
selected a rate of change of width that depends upon the relative total
humidity of the grid box, such that there is more erosion in drier
gridboxes. This promotes more rapid erosion of shallow convective cloud,
which is the main effect that we seek to include, although the physical
implication that dry air is more turbulent than moist air does not match
the way the real atmosphere works, especially in the stratosphere. No
doubt the link can be improved upon with more research. The formulation
used is:

.. math::

   \frac{1}{b_s} \frac{\partial b_s}{\partial t} = \Upsilon exp ( - \frac{2.01 Q_c}{0.2 a_L q_{sat liq}(T_L)} )
   \label{eq:dbsbydtbs_turb}

where the 0.2 factor is chosen to be closely equivalent to
:math:`1 - RH_{crit}` and the value of 2.01 has been selected through
tuning. The code merges the two numerical values into a single quantity
(dbsdtbs-exp), equal to 10.05. We note that in PC2:64 (the library 6.4
code, a value of 0.62 is used rather than 2.01). As a guide to the
:math:`RH_T` dependence, note that when the value of :math:`RH_T` is
0.85, the value of :math:`\frac{1}{b_s} \frac{\partial b_s}{\partial t}`
is close to :math:`1 \times 10^{-4} s^{-1}`.

Note: the source-code for this erosion method (**pc2_hom_conv**,
**pc2_homog_plus_turb**, **pc2_delta_hom_turb**) includes an additional
term “dbsdtbs1” which scales with the rate of homogeneous forcing
:math:`\frac{\partial Q_c}{\partial t}`. However this term is always set
to zero on input to these routines so is never used.

The width-narrowing formulation of section `3.3 <#sec:width>`__ is used
to calculate increments in :math:`\overline{q_{cl}}` and :math:`C_l`.
Using the liquid - ice cloud overlap ideas of section `3.6 <#sec:ct>`__
then gives the associated :math:`C_t` change. This background narrowing
term, :math:`\Upsilon`, is originally based upon work by
:raw-latex:`\cite{sg03}`, although it is a parameter that has been
extensively tuned during PC2 development, a typical value would be
:math:`\Upsilon=-2.25 \times 10^{-5} s^{-1}`.

Numerical application of the original width-narrowing method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because of the strong link the mathematical expressions for width
narrowing (section `3.3 <#sec:width>`__) have with the expressions for
the homogeneous forcing (section `3.2 <#sec:homog>`__), we choose to
represent the timestepping of this process in exactly the same way as
for the homogeneous forcing (in fact, in the Unified Model code we use
the same subroutine, see section `5.2 <#sec:code>`__). As before, we use
a simple forward timestepping of :math:`C_l`, with :math:`Q_c` given by
(`[eq:qc_eq_qt-qs] <#eq:qc_eq_qt-qs>`__) and :math:`a_L` defined as
discussed in section `3.2.3 <#sec:homog_num_app>`__ and discretize eq
`[eq:dcdt_width] <#eq:dcdt_width>`__ as:

.. math::

   \Delta C_l^{[n+1]} = - G(-Q_c) Q_c \frac{1}{b_s}
   \frac{\partial b_s}{\partial t} \Delta t.
   \label{eq:dcl_turb_final}

Similarly to (`[eq:c_l^n+1] <#eq:c_l^n+1>`__), we then limit the cloud
fraction to 0 and 1 and then apply a mid-point value of :math:`C_l` to
calculate the change in :math:`\overline{q_{cl}}` (discretizing eq
`[eq:dqcldt_width] <#eq:dqcldt_width>`__):

.. math::

   \Delta q_{cl}^{[n+1]} = (q_{cl}^{[n]} - Q_c \frac{1}{2}(C_l^{[n]}+C_l^{[n+1]}))
   \frac{1}{b_s} \frac{\partial b_s}{\partial t} \Delta t.
   \label{eq:dqcl_turb_final}

In this case the value of :math:`\Delta q_{cl}` *is* limited to ensure
that no more :math:`\overline{q_{cl}}` is removed than the model has
available. This was chosen to ensure that the erosion process itself
contains this physical limit, not a numerical tidying-up process.

The option “l_fixbug_pc2_qcl_incr” ensures that qcl is set to zero if
the CFL has reached zero.

Cloud-surface-area hybrid erosion method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(selected by setting **i_pc2_erosion_method = 3** in the UM namelist).

:raw-latex:`\cite{morcrette_petch}` showed that changes to the erosion
parameter (:math:`\Upsilon` in Eqn.
`[eq:dbsbydtbs_turb] <#eq:dbsbydtbs_turb>`__) did not have as
significant an impact on the global work done by the erosion process as
might be expected. This was due to a feedback process whereby, reducing
the erosion parameter leads to more cloud water, more autoconversion of
cloud water to rain, more fall-out of rain and more drying of the layer,
hence increasing the
:math:`exp ( - \frac{2.01 Q_c}{0.2 a_L q_{sat liq}(T_L)} )` part of Eqn.
`[eq:dbsbydtbs_turb] <#eq:dbsbydtbs_turb>`__. Although the feedback is
physically plausible it crucially depends on the formulation of Eqn.
`[eq:dbsbydtbs_turb] <#eq:dbsbydtbs_turb>`__ and the dependence of the
rate of narrowing of the PDF on the moisture, a dependence that was
developed from a pragmatic rather than theoretical stand-point. The
option for an alternative way of calculating the erosion was introduced
at vn8.0

We use equation 30 from :raw-latex:`\cite{t93}` to specify the sink of
:math:`q_{cl}` due to erosion, i.e.

.. math::

   \frac{\partial q_{cl}}{\partial t}=-A K(q_{sat}-q_v)
   \label{eq:dqcldt_hybrid}

(note we have changed the sign as we have replaced the evaporation rate
:math:`E_2` in :raw-latex:`\cite{t93}` with
:math:`-\frac{\partial q_{cl}}{\partial t}` on the left-hand-side). In
the :raw-latex:`\cite{t93}` scheme, :math:`A` is set to the cloud
fraction (i.e. :math:`A=C_l`). Here we recall that "cloud erosion" is
meant to represent the evaporation of cloud water due to the mixing of
clear and cloudy air and that this can only happen on the edges of
cloud, where saturated air is exposed to sub-saturated air. If the cloud
fraction is small (e.g. 5\ :math:`\%`), then there are not many clouds,
so there is only a small surface area from which evaporation can occur.
Similarly if the cloud cover is very high (e.g. 95\ :math:`\%`) then
there is again not much surface area exposed to clear sky. A maximum in
exposed surface area is expected when the cloud cover is 50\ :math:`\%`.

By imagining that the grid-box is broken up into cubes whose horizontal
dimension equal the layer depth it is possible to work out what the
maximum lateral surface area would be, as a function of cloud fraction,
for different arrangements of cloudy cubes. The maximum lateral surface
area, is when the clear and cloudy cubes are arranged in a chess-board
pattern, and the minimum is when then are all grouped together into a
circular clump. Numerical tests using randomly distributed cloudy cube
shows that the variation in lateral surface area, :math:`S`, as a
function of cloud fraction can be expressed as:

.. math::

   S= - 2 C_l ^{2} + 2 C_l
   \label{eq:S_Cl}

The maximum normalised surface area of 0.5 occurs at a cloud fraction of
0.5. Using a cloud mask derived from satellite imagery shows that real
cloud fields do follow this kind of dependence, but that the peak
surface area is nearer to 0.35, meaning that real clouds are not as
randomly distributed as random ones and that there is some kind of
clumping together, which is what we might have expected. When it comes
to implementing such a scheme in the model, there will need to be a
tunable parameter to govern the rate of evaporation. This will not
affect the shape of the lateral surface area function. As a result the
details of whether the peak lateral surface area is 0.5 or 0.35 are
simply absorbed into the tunable parameter :math:`K`, which is supplied
from the UMUI (using the same text box as was used for supplying
:math:`\Upsilon`).

The exposed surface area associated with the tops and bottom of the
clouds is calculated assuming maximum overlap in adjacent layers and is
added to the lateral surface area to give a total surface area,

.. math::

   A=max(C_l(k)-C_l(k+1),0.0)+max(C_l(k)-C_l(k-1),0.0)+S
   \label{eq:A_top_and_bottom}

it is this value of :math:`A` which we use in Eqn.
`[eq:dqcldt_hybrid] <#eq:dqcldt_hybrid>`__.

Note that the contributions from the top and bottom interfaces of the
current model-level :math:`max(C_l(k)-C_l(k+1),0.0)` and
:math:`max(C_l(k)-C_l(k-1),0.0)` may optionally either be included or
excluded, depending on the UM namelist switch **i_pc2_erosion_method**.
Further note: at present these contributions are hardwired to be
excluded, as they prevented the erosion calculations from being
parallelised in the vertical direction using OpenMP, and no operational
model configurations were using them.

Having calculated a reduction in :math:`q_{cl}` using the
Tiedtke-surface-area method, we then work out the relative rate of
narrowing that would have given the same sink of :math:`q_{cl}`. This
value of :math:`\frac{1}{b_s} \frac{\partial b_s}{\partial t}` is then
used to calculate the change in :math:`C_l` using the same moisture PDF
assumptions as were used in the original PC2 erosion formulation. To
achieve this, we combine equations `[eq:dcdt_width] <#eq:dcdt_width>`__
and `[eq:dqcldt_width] <#eq:dqcldt_width>`__ from section
`3.3 <#sec:width>`__ to eliminate
:math:`\frac{1}{b_s} \frac{\partial b_s}{\partial t}` and write
:math:`\frac{\partial C_l}{\partial t}` as a function of
:math:`\frac{\partial \overline{q_{cl}}}{\partial t}`:

.. math::

   \frac{\partial C_l}{\partial t}
    = - \frac{ G(-Q_c) Q_c \frac{\partial \overline{q_{cl}}}{\partial t} }
             { (- C_l Q_c+\overline{q_{cl}}) }
   \label{eq:dcdt_hybrid}

Where the change in liquid water content
:math:`\frac{\partial \overline{q_{cl}}}{\partial t}` is given by eq
`[eq:dqcldt_hybrid] <#eq:dqcldt_hybrid>`__ above.

This combination of a Tiedkte sink term for :math:`q_{cl}`, a PC2 term
for :math:`C_l` and the introduction of some surface area dependence
leads to this formulation being referred to as a “hybrid”
cloud-surface-area erosion method.

.. _`sec:erosion_numerics`:

Numerical application of the hybrid erosion method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, we consider how to numerically discretise equations
`[eq:dqcldt_hybrid] <#eq:dqcldt_hybrid>`__ and
`[eq:dcdt_hybrid] <#eq:dcdt_hybrid>`__ to compute cloud increments due
to erosion. The simplest approach is an explicit forwards-in-time
discretisation:

.. math::

   \frac{ \Delta {q_{cl}}_{ero}}{\Delta t} = A(C_l^n) K(q_{sat}-q_v)
   \label{eq:hybrid_erosion_expl}

i.e. the increment is calculated by evaluating the term :math:`A` from
equations `[eq:S_Cl] <#eq:S_Cl>`__ and
`[eq:A_top_and_bottom] <#eq:A_top_and_bottom>`__ using the value of
cloud-fraction :math:`C_l` *before* erosion has been applied.

However, when the environment is significantly subsaturated (so that the
term :math:`(q_{sat}-q_v)` is large and negative), and long timesteps
:math:`\Delta t` are used (e.g. order 1000 s used in global climate
simulations), this discretization can suffer severe numerical overshoot.
i.e. the increment based on :math:`C_l^n` is large enough to reduce
:math:`q_{cl}` (and hence also :math:`C_l`) to less than zero within a
single timestep. If the continuous equation were solved analytically
this wouldn’t happen; as :math:`C_l` declines due to the erosion, so
will :math:`A(C_l)` and hence the erosion rate, so that :math:`q_{cl}`
and :math:`C_l` smoothly decline towards zero.

Three options are available in the code to address this problem,
selected by the UM namelist switch **i_pc2_erosion_numerics**, detailed
below. Single-Column Model tests indicate that the 2nd and 3rd options
yield much less timestep sensitivity for detrained cloud in shallow
cumulus regimes.

#. **Retain the explicit discretization, but limit the resulting erosion
   increments to ensure :math:`q_{cl}` and :math:`C_l` don’t go
   negative. (i_pc2_erosion_numerics=1)** Also, to ensure that some
   cloud remains at end-of-timestep where shallow cumulus is detraining
   into dry environments, the erosion calculation is fed copies of the
   fields with the current timestep’s convection increments subtracted
   off. This means any cloud detrained by convection during the current
   timestep cannot be eroded until the following timestep, and so is
   still present at end-of-timestep. As discussed in section
   `7.4.11 <#sec:timestepping>`__, this leads to a problematic timestep
   sensitivity, since the amount of cloud not subject to erosion is the
   convection increment, which scales with the timestep length.

   Having computed the erosion :math:`q_{cl}` increment using
   `[eq:hybrid_erosion_expl] <#eq:hybrid_erosion_expl>`__, the
   consistent :math:`C_l` increment is computed by discretising
   `[eq:dcdt_hybrid] <#eq:dcdt_hybrid>`__ as:

   .. math::

      \frac{\Delta {C_l}_{ero}}{\Delta t}
       = - \frac{ G(-Q_c)^n Q_c^n \frac{\Delta \overline{{q_{cl}}_{ero}}}{\Delta t} }
                { (- C_l^n Q_c^n + ( \overline{q_{cl}^n}
                                  + \frac{1}{2} \Delta \overline{{q_{cl}}_{ero}} ) ) }
      \label{eq:dcdt_hybrid_discr}

   i.e. all terms are treated explicitly (using the values before
   erosion), except for :math:`\overline{q_{cl}}` which takes the
   mid-point interpolated half-way between its values before and after
   erosion, to give some improvement in accuracy.

#. **Use an approximate implicit discretisation, which intrinsically
   yields a positive solution for :math:`q_{cl}` and :math:`C_l`.
   (i_pc2_erosion_numerics=2)** The copies of the fields passed to the
   erosion calculation are fully updated with the convection increments.
   We then write equation `[eq:dqcldt_hybrid] <#eq:dqcldt_hybrid>`__ in
   the form:

   .. math:: \frac{\partial q_{cl}}{\partial t} = q_{cl} f(q_{cl},C_l,(q_{sat}-q_v))

   (where the term
   :math:`f(q_{cl},C_l,(q_{sat}-q_v)) = \frac{A K(q_{sat}-q_v)}{q_{cl}}`
   will be treated explicitly, under the assumption that this ratio will
   evolve more slowly while erosion rapidly reduces both the numerator
   and the denominator).

   We then take a backwards-in-time implicit discretisation in terms of
   the leading factor :math:`q_{cl}`:

   .. math:: \frac{ q_{cl}^{n+1} - q_{cl}^{n}}{\Delta t} = q_{cl}^{n+1} f^n

   Now, the problem is somewhat complicated by the fact that in the
   code, erosion is calculated in parallel with the homogeneous forcing
   by convection, and we need to account for the homogeneous forcing
   increment :math:`\Delta q_{cl}^{hom}` in our implicit solution. We
   therefore write the above as:

   .. math::

      \Delta q_{cl}^{ero} = \Delta t
                         ( q_{cl}^n + \Delta q_{cl}^{hom} + \Delta q_{cl}^{ero} ) f^n

   Rearranging:

   .. math::

      \Delta q_{cl}^{ero} = \Delta t f^n q_{cl}^n \frac{ q_{cl}^n + \Delta q_{cl}^{hom} }
                                                 { q_{cl}^n - \Delta t f^n q_{cl}^n }

   Note that the term :math:`\Delta t f^n q_{cl}^n` is the erosion
   increment we would obtain from the purely explicit discretisation,
   :math:`\Delta q_{cl}^{ero\,expl}`. The implicit discretisation is
   implemented by first calculating :math:`\Delta q_{cl}^{ero\,expl}`
   using equation `[eq:hybrid_erosion_expl] <#eq:hybrid_erosion_expl>`__
   (as we do for **i_pc2_erosion_numerics=1**) but then rescaling it
   using the above expression, which becomes:

   .. math::

      \Delta q_{cl}^{ero} = \Delta q_{cl}^{ero\,expl} \frac{ q_{cl}^n + \Delta q_{cl}^{hom} }
                                                 { q_{cl}^n - \Delta q_{cl}^{ero\,expl} }
      \label{eq:hybrid_erosion_impl_qcl}

   Provided erosion is acting to reduce cloud-water
   (:math:`\Delta q_{cl}^{ero\,expl} < 0`), and homogeneous forcing by
   convection has not already completely removed the cloud
   (:math:`q_{cl}^n + \Delta q_{cl}^{hom} > 0`),
   `[eq:hybrid_erosion_impl_qcl] <#eq:hybrid_erosion_impl_qcl>`__ is
   guaranteed to yield a stable, positive solution for :math:`q_{cl}`.

   We also apply exactly the same argument to the equation for the
   cloud-fraction increment :math:`C_l`, and obtain:

   .. math::

      \Delta C_l^{ero} = \Delta C_l^{ero\,expl} \frac{ C_l^n + \Delta C_l^{hom} }
                                                 { C_l^n - \Delta C_l^{ero\,expl} }
      \label{eq:hybrid_erosion_impl_Cl}

   Where :math:`\Delta C_l^{ero\,expl}` is computed using eq
   `[eq:dcdt_hybrid_discr] <#eq:dcdt_hybrid_discr>`__, except that the
   term :math:`\frac{1}{2} \Delta \overline{{q_{cl}}_{ero}}` is omitted
   (interpolating to the mid-point value of :math:`\overline{q_{cl}}` in
   the denominator would be “double-counting” if we are already making
   an implicit correction to the full increment).

   In the case where the homogeneous forcing increments have already
   removed all of the cloud water content or fraction, erosion is not
   performed, and :math:`q_{cl}` and :math:`C_l` are both set to zero.
   In the case where erosion is actually acting to increase
   cloud-fraction, the code defaults to retaining the explicit
   discretisation solution :math:`\Delta q_{cl}^{ero\,expl}` and
   :math:`\Delta C_l^{ero\,expl}`. Otherwise, equations
   `[eq:hybrid_erosion_impl_qcl] <#eq:hybrid_erosion_impl_qcl>`__ and
   `[eq:hybrid_erosion_impl_Cl] <#eq:hybrid_erosion_impl_Cl>`__ are
   applied to yield the implicit solution.

#. **Use an analytic solution to the integration of the time-derivatives
   in (**\ `[eq:dqcldt_hybrid] <#eq:dqcldt_hybrid>`__\ **) and
   (**\ `[eq:dcdt_hybrid] <#eq:dcdt_hybrid>`__\ **) for greater
   accuracy. (i_pc2_erosion_numerics=3)**

   Two problems have been identified with the above implicit numerical
   method:

   - The implicit correction is applied completely independently to the
     increments for :math:`q_{cl}` and :math:`C_l`. So as with the
     explicit method, differing numerical error in the increments for
     the two variables can lead to them becoming inconsistent with
     eachother. It was found by experimentation that even with the
     implicit correction, it is possible for erosion to reduce
     :math:`C_l` by a bigger fraction than :math:`q_{cl}`, so that the
     in-cloud water content :math:`\frac{q_{cl}}{C_l}` is *increased*.
     Narrowing the PDF should only *decrease* the in-cloud
     water-content; occasional large increases due to numerical error
     can lead to spurious precipitation being produced by the
     microphysics scheme.

   - The implicit correction makes it impossible for erosion to reduce
     :math:`q_{cl}`, :math:`C_l` to zero. As we will show below, the
     analytic solution to the equations posed does in fact go to zero
     after a finite time under grid-mean subsaturation (although the
     erosion rate declines with :math:`C_l` as it approaches zero,
     :math:`C_l` approaches zero more slowly than :math:`q_{cl}`, so
     that both variables decrease following power-law curves not
     exponentials). When erosion (wrongly) can never entirely remove
     cloud, this allows tiny values of :math:`q_{cl}` and :math:`C_l` to
     spuriously spread across the domain via numerical diffusion from
     the model’s advection scheme.

   Under this option, we attempt to compute an analytic solution to the
   simultaneous differential equations
   `[eq:dqcldt_hybrid] <#eq:dqcldt_hybrid>`__ and
   `[eq:dcdt_hybrid] <#eq:dcdt_hybrid>`__ so that :math:`q_{cl}` and
   :math:`C_l` both decrease smoothly and consistently. The equations
   lead to somewhat different behaviour depending on whether the
   grid-mean state is subsaturated (:math:`Q_c < 0`), supersaturated
   (:math:`Q_c > 0`), or close to saturation (:math:`Q_c` near-zero). We
   can employ different approximations to integrate the equations in
   each case. In the code, we first test the value of :math:`Q_c` and
   compute the erosion increments as follows:

   #. **Grid-mean subsaturation (:math:`Q_c < 0`):**

      The relation between the erosion tendencies in
      liquid-cloud-fraction and liquid water content
      (`[eq:dcdt_hybrid] <#eq:dcdt_hybrid>`__) can be expressed in terms
      of *fractional* rates of change (dividing the top and bottom by
      :math:`-C_l Q_c`, and dividing both sides by :math:`C_l`):

      .. math::

         \frac{1}{C_l} \frac{\partial C_l}{\partial t}
          = \frac{ G(-Q_c) \frac{q_{cl}}{C_l^2} }{ 1 - \frac{q_{cl}}{C_l Q_c} } \;
            \frac{1}{q_{cl}} \frac{\partial q_{cl}}{\partial t}
         \label{eq:dcdt_hybrid_1}

      Under homogeneous forcing (section `3.2 <#sec:homog>`__), we
      defined the PDF height at the saturation boundary when near the
      cloudy end of the PDF as
      :math:`G(-Q_c) = \frac{n+1}{n+2} \frac{C_l^2}{q_{cl}}` (eq
      `[eqn20] <#eqn20>`__). In fact, :math:`G(-Q_c)` is set to some
      blend between this and the value near the clear end of the PDF (eq
      `[eqn21] <#eqn21>`__). But we will assume that when eroding cloud
      under grid-mean subsaturated conditions (:math:`Q_c < 0`),
      :math:`G(-Q_c)` follows this scaling with
      :math:`\frac{C_l^2}{q_{cl}}` even if its value differs somewhat
      from eq `[eqn20] <#eqn20>`__. Therefore the quantity
      :math:`c_1 = G(-Q_c) \frac{q_{cl}}{C_l^2}` remains constant during
      the erosion process, and eq
      `[eq:dcdt_hybrid_1] <#eq:dcdt_hybrid_1>`__ becomes:

      .. math::

         \frac{1}{C_l} \frac{\partial C_l}{\partial t}
          = \frac{ c_1 }{ 1 - \frac{q_{cl}}{C_l Q_c} } \;
            \frac{1}{q_{cl}} \frac{\partial q_{cl}}{\partial t}
         \label{eq:dcdt_hybrid_2}

      The term :math:`1 - \frac{q_{cl}}{C_l Q_c}` (which is :math:`> 1`
      since we are considering grid-mean subsaturation :math:`Q_c < 0`)
      usually remains close to 1 in practice, so we can assume its
      fractional variation over the timestep is small compared to the
      other terms, and treat it explicitly. We can therefore
      straightforwardly integrate eq
      `[eq:dcdt_hybrid_2] <#eq:dcdt_hybrid_2>`__ to obtain the scaling
      of :math:`C_l` with :math:`q_{cl}` as both are reduced by erosion:

      .. math::

         \frac{C_l}{{C_l}_0} = \left( \frac{q_{cl}}{{q_{cl}}_0} \right)^{b_1}
         \label{eq:cl_qcl_scaling}

      where :math:`{C_l}_0`, :math:`{q_{cl}}_0` are the values before
      erosion is applied, and the exponent is
      :math:`b_1 = \frac{ c_1 }{ 1 - \frac{q_{cl}}{C_l Q_c} }`. When
      :math:`G(-Q_c)` takes its value from the cloudy end of the PDF, we
      have :math:`c_1 = \frac{n+1}{n+2}`. Since the PDF power
      :math:`n > 0` and :math:`Q_c < 0` under the considered grid-mean
      subsaturation, we always have :math:`b_1 < 1`. This ensures that
      erosion reduces :math:`C_l` at a slower fractional rate than
      :math:`q_{cl}`, so that in-cloud water content
      :math:`\frac{q_{cl}}{C_l}` always decreases.

      Next we derive an integral solution for the decline of
      :math:`q_{cl}` with time. Ignoring the cloud surface-area
      contributions from the levels above and below (they are disabled
      in the code anyway), the erosion liquid water content tendency is
      obtained by combining `[eq:dqcldt_hybrid] <#eq:dqcldt_hybrid>`__
      and `[eq:S_Cl] <#eq:S_Cl>`__:

      .. math::

         \frac{\partial q_{cl}}{\partial t} = -K \, 2 C_l (1 - C_l) \, (q_{sat}(T)-q_v)
         \label{eq:dqcldt_hybrid_1}

      From eq `[SD2] <#SD2>`__, :math:`q_{sat}(T)-q_v = \frac{SD}{a_L}`,
      where :math:`SD` is the saturation defecit, and :math:`a_L` is the
      dimensionless factor defined in eq `[eq:a_L] <#eq:a_L>`__.
      Following the derivation in section
      `4.9.3 <#sec:smooth_initiation>`__ (eq
      `[eq:qc_plus_sd] <#eq:qc_plus_sd>`__), we can write this in terms
      of the liquid-water content: :math:`SD = q_{cl} - Q_c` (where
      :math:`Q_c` was defined in eq
      `[eq:qc_eq_qt-qs] <#eq:qc_eq_qt-qs>`__, and corresponds to the
      grid-mean supersaturation converted to an equivalent liquid water
      content). Substituting this into
      (`[eq:dqcldt_hybrid_1] <#eq:dqcldt_hybrid_1>`__) above, we obtain:

      .. math::

         \frac{\partial q_{cl}}{\partial t} = -\frac{K}{a_L} \, 2 C_l (1 - C_l) \,
                                             (q_{cl}-Q_c)
         \label{eq:dqcldt_hybrid_2}

      Substituting eq `[eq:cl_qcl_scaling] <#eq:cl_qcl_scaling>`__ for
      the leading factor of :math:`C_l` on the right-hand-side and
      rearranging:

      .. math::

         \left( \frac{q_{cl}}{{q_{cl}}_0} \right)^{-b_1} \frac{\partial q_{cl}}{\partial t}
            = -\frac{K}{a_L} \, 2 {C_l}_0 (1 - C_l) \, (q_{cl}-Q_c)

      In significantly subsaturated conditions the r.h.s. has only weak
      dependence on :math:`C_l` and :math:`q_{cl}` (:math:`C_l << 1`,
      :math:`q_{cl} << -Q_c`), so we can treat the whole r.h.s.
      explicitly (i.e. neglect its variation during each timestep), so
      that the above integrates to:

      .. math::

         \left[ \frac{{q_{cl}}_0}{1-b_1} \left( \frac{q_{cl}}{{q_{cl}}_0} \right)^{1-b_1}
         \right]_{{q_{cl}}_0}^{{q_{cl}}_{\Delta t}}
          = -\frac{K}{a_L} \, 2 {C_l}_0 (1 - C_l) \, (q_{cl}-Q_c) \Delta t

      Inserting the limits of the integral on the l.h.s. and
      rearranging, we obtain our analytical solution for :math:`q_{cl}`
      after time :math:`\Delta t`:

      .. math::

         {q_{cl}}_{\Delta t} = {q_{cl}}_0 \left( 1 - \frac{1-b_1}{{q_{cl}}_0}
             \frac{K}{a_L} \, 2 {C_l}_0 (1 - C_l) \, (q_{cl}-Q_c) \Delta t
                                     \right)^\frac{1}{1-b_1}
         \label{eq:qcl_int_hybrid}

      Note that :math:`q_{cl}` falls to zero after a finite time
      :math:`\frac{{q_{cl}}_0}{1-b_1} \frac{a_L}{K}
       \frac{1}{2 {C_l}_0 (1 - C_l) \, (q_{cl}-Q_c)}`. If the timestep
      :math:`\Delta t` is longer than this time, then erosion completely
      removes the cloud during the current timestep.

      We first set :math:`{q_{cl}}_0` and :math:`{C_l}_0` to the values
      already updated by homogeneous forcing, and then sequentially
      compute the updated :math:`q_{cl}` after erosion using
      (`[eq:qcl_int_hybrid] <#eq:qcl_int_hybrid>`__). Then we substitute
      this value into (`[eq:cl_qcl_scaling] <#eq:cl_qcl_scaling>`__) to
      compute the consistent updated value of :math:`C_l`. Finally, to
      improve accuracy, a small number of iterations are performed to
      find the solution with the explicitly-treated terms (the exponent
      :math:`b_1 = \frac{ c_1 }{ 1 - \frac{q_{cl}}{C_l Q_c} }` and the
      terms :math:`(1 - C_l)` and :math:`(q_{cl}-Q_c)` in eq
      `[eq:qcl_int_hybrid] <#eq:qcl_int_hybrid>`__) adjusted to values
      linearly-interpolated to half-way between the start and end of the
      erosion timestep.

   #. **Grid-mean supersaturation (:math:`Q_c > 0`):**

      In this case, erosion does not act to reduce :math:`C_l` and
      :math:`q_{cl}` towards zero. Instead, the narrow PDF limit it
      adjusts towards has no remaining subsaturated air, so that
      :math:`C_l = 1` and :math:`q_{cl} = Q_c`. In this case, we can
      repeat the above derivation, but considering the equations for the
      clear-fraction :math:`1-C_l` in place of :math:`C_l`, and the
      saturation defecit :math:`SD = q_{cl}-Q_c` in place of
      :math:`q_{cl}`. Assuming that :math:`G(-Qc)` follows the scaling
      for the clear end of the PDF (`[eqn21] <#eqn21>`__), this leads to
      a similar equation to
      (`[eq:cl_qcl_scaling] <#eq:cl_qcl_scaling>`__) but for the scaling
      as erosion reduces :math:`1-C_l` and :math:`SD` towards zero .:

      .. math::

         \frac{1-C_l}{1-{C_l}_0} = \left( \frac{SD}{SD_0} \right)^{b_2}
         \label{eq:ca_sd_scaling}

      with :math:`b_2 = \frac{ c_2 }{ 1 + \frac{SD}{(1-C_l) Q_c} }` and
      :math:`c_2 = G(-Q_c) \frac{SD}{(1-C_l)^2}` (note we must have
      :math:`0 < b_2 < 1`).

      And then the tendency equation for :math:`SD` is:

      .. math::

         \frac{\partial SD}{\partial t} = -\frac{K}{a_L} \, 2 (1 - C_l) C_l \, SD
         \label{eq:dsddt_hybrid_2}

      The one asymmetry between this and the :math:`q_{cl}` tendency
      equation (`[eq:dqcldt_hybrid_2] <#eq:dqcldt_hybrid_2>`__) is that
      for :math:`SD` the r.h.s. is directly proportional to the quantity
      in the time-derivative, whereas for :math:`q_{cl}` there is an
      additional :math:`Q_c` term which is constant during erosion.
      Substituting (`[eq:ca_sd_scaling] <#eq:ca_sd_scaling>`__) for the
      leading factor of :math:`(1 - C_l)` in
      (`[eq:dsddt_hybrid_2] <#eq:dsddt_hybrid_2>`__), integrating over
      time :math:`\Delta t` (neglecting the fractional variation of
      :math:`C_l` over the timestep) and rearranging, we obtain:

      .. math::

         {SD}_{\Delta t} = {SD}_0 \left( 1 + b_2
             \frac{K}{a_L} \, 2 (1-{C_l}_0) C_l \, \Delta t
                               \right)^{-\frac{1}{b_2}}
         \label{eq:sd_int_hybrid}

      Note that the additional power of :math:`SD` on the r.h.s. of
      (`[eq:dsddt_hybrid_2] <#eq:dsddt_hybrid_2>`__) leads to the
      integral solution having a negative exponent. This means that
      under grid-mean supersaturation, erosion makes :math:`SD` and
      :math:`1-C_l` approach but never quite reach zero, which is quite
      different behaviour to grid-mean subsaturation where
      :math:`q_{cl}` and :math:`C_l` go to zero over a finite time. This
      asymmetry is because the erosion rate is parameterised to be
      proportional to :math:`SD`, and this tends to zero as the PDF is
      narrowed under supersaturation, but remains finite positive under
      subsaturation.

      We first set :math:`{SD}_0 = {q_{cl}}_0 - Q_c` (where as above
      :math:`{q_{cl}}_0` is the value already updated by homogeneous
      forcing), then compute the value of :math:`SD` updated by erosion
      using (`[eq:sd_int_hybrid] <#eq:sd_int_hybrid>`__). Then we
      substitute this value into
      (`[eq:ca_sd_scaling] <#eq:ca_sd_scaling>`__) to compute the
      consistent updated value of :math:`1-C_l`. A small number of
      iterations are then performed to find the solution with the
      explicitly-treated terms (the exponent
      :math:`b_2 = \frac{ c_2 }{ 1 + \frac{SD}{(1-C_l) Q_c} }` and the
      term :math:`C_l` in eq `[eq:sd_int_hybrid] <#eq:sd_int_hybrid>`__)
      adjusted to values linearly-interpolated to half-way between the
      start and end of the erosion timestep. Then the final values of
      :math:`SD` and :math:`1-C_l` are used to increment
      :math:`q_{cl} = Q_c + SD` and :math:`C_l`, as prognosed by the
      rest of the model.

   #. **grid-mean saturation (:math:`Q_c` near-zero):**

      In this case, the PDF is centred on the saturation boundary, so
      that narrowing it does not change the cloud-fraction. In the limit
      :math:`Q_c = 0`, we have :math:`q_{cl} = SD`, and
      (`[eq:dqcldt_hybrid_2] <#eq:dqcldt_hybrid_2>`__) or
      (`[eq:dsddt_hybrid_2] <#eq:dsddt_hybrid_2>`__) becomes:

      .. math::

         \frac{1}{q_{cl}} \frac{\partial q_{cl}}{\partial t}
           = -\frac{K}{a_L} \, 2 C_l (1 - C_l)

      where everything on the r.h.s. is constant under erosion. This
      simply integrates to give exponential decline of :math:`q_{cl}`
      (and :math:`SD`) towards zero:

      .. math:: {q_{cl}}_{\Delta t} = {q_{cl}}_0 e^{ -\frac{K}{a_L} \, 2 C_l (1 - C_l) \Delta t }

Orographic and Gravity Wave Drag
--------------------------------

The Orographic and Gravity Wave Drag sections do not alter the
temperature or moisture content of the model gridboxes, hence PC2
assumes no change in the condensate and cloud fractions as a result of
these processes.

.. _`sec:advec`:

Advection
---------

The advection of :math:`\overline{q_{cl}}` and :math:`\overline{q_{cf}}`
are already performed separately by the semi-Lagrangian advection
scheme. Advection of the three cloud fractions :math:`C_l`, :math:`C_i`
and :math:`C_t` are all performed by PC2 in the same way.

Note that ascent or subsidence by advection entails a pressure change
following each parcel, which will cause an accompanying adiabatic
temperature change. These advective pressure and temperature changes
imply a homogeneous forcing, which yields a change in
:math:`\overline{q_{cl}}` and :math:`C_l` in addition to their transport
by the winds. This is described in section `4.8 <#sec:pres>`__.

If the UM namelist switch **l_pc2_sl_advection** is turned on, the PC2
homogeneous forcing response to advection is calculated straight after
the call to Semi-Lagrangian advection. Otherwise, the pressure change
from advection is combined with the Eulerian pressure change from the
dynamics Helmholtz solver, and the resulting homogeneous forcing of
liquid cloud is computed at the end of the timestep.

.. _`sec:bl`:

Boundary Layer
--------------

At a basic level, the boundary layer scheme works by mixing
:math:`\overline{q_T}` and :math:`\overline{T_L}`, and tracer mixing
:math:`\overline{q_{cf}}`. The condensation and :math:`C_l` changes are
represented using the homogeneous forcing representation. The forcing of
:math:`Q_c` can be written in :math:`\Delta \overline{q_T}` and
:math:`\Delta \overline{T_L}` terms using
(`[eq:deltaqc_exp] <#eq:deltaqc_exp>`__).

:math:`\overline{q_{cf}}` is already mixed using the tracer mixing
scheme. PC2 will calculate the corresponding :math:`C_i` change assuming
the inhomogeneous forcing scenario. Although this is not necessarily an
appropriate physical model to use, it is the only generic physical model
we have currently developed in order to convert increments in a
condensate to increments in a cloud fraction. We use a value of the
in-cloud water content :math:`q_c^S` based upon a linear combination of
the current in-cloud ice water content,
:math:`\frac{\overline{q_{cf}}}{C_i}`, and a fixed value.

.. math::

   q_C^S = C_i \frac{\overline{q_{cf}}}{C_i} + ( 1 - C_i) q_{cf0 \, BL}
   \label{eq:qcf_ci}

where :math:`q_{cf0 \, BL}` is a specified value of
:math:`1 \times 10^{-4} \, kg \, kg^{-1}`. We then use the inhomogeneous
forcing equation based upon (`[eq:dcdt_inhom2] <#eq:dcdt_inhom2>`__) but
for ice water content to write

.. math::

   \Delta C_i = \frac{(1 - C_i)}{q_C^S - \overline{q_{cf}}} Q4_i .
   \label{eq:deltaci_bl}

Since the physical model will have :math:`C_i` tend to 1 if the
denominator is small, we will, to avoid numerical problems, set
:math:`C_i` to 1 if
:math:`q_C^S - \overline{q_{cf}} < 1 \times 10^{-10} kg kg^{-1}`. Note
that we do not use the multiple phases injection source expressions
(section `3.5.1 <#sec:multiple>`__ and equation
`[eq:cff_ts] <#eq:cff_ts>`__). This is because the liquid water changes
are not associated with the plume model.

Equation `[eq:qcf_ci] <#eq:qcf_ci>`__ assumes that the change to the ice
water content has led to an increase in ice water content. However, if
the ince water content has reduced, the change to the ice cloud fraction
is not consistent. The option to "Use consistent formulation of ice
cloud fraction changes due to boundary-layer processes" ensure that if
the ice water content is reduced, the ice cloud fraction is reduced, in
such as way as to maintain the same in-cloud ice water content.

The :math:`C_t` changes are calculated using the minimum overlap method
of section `3.6 <#sec:ct>`__.

In *ni-imp-ctl* the control code inhibits the call to the diagnostic
cloud scheme if there is deep or shallow convection occurring and the
model level is less than *or equal to* the layer immediately above the
top of the boundary layer mixed layer (i.e. level ntml+1). This is in
order to ensure that there is no large-scale cloud present below the
base of the convective cloud, but additionally performs this calculation
at the level above, probably in order that latent heating from
large-scale condensation does not inhibit the convection. A similar
thing is performed for PC2, with any large-scale cloud being evaporated
if the same criteria are met, *except that it is not performed on the
level above the boundary layer mixed layer*. This choice (i.e. ntml) is
seen to give improved results in PC2, and is arguably a more physical
reasonable choice anyway than using ntml+1.

.. _`sec:convec`:

Convection
----------

This section concentrates specifically upon the PC2 interface to the
convection scheme. In the current formulation of the UM, only a
mass-flux convection scheme exists, and this is what is described below.
Work to interface PC2 to the developing turbulence based convection
scheme is commented upon in section `7.4.4 <#sec:tbcs>`__.

An alternative way of calculating cloud fraction increments is currently
under development and is described in section
`4.7.10 <#sec:conv-simpler>`__.

A traditional view of convective parametrization is a scheme that
transports vapour, :math:`q`, heat, :math:`\theta`, and momentum,
:math:`u` and :math:`v` winds within a single column. It does not
consider transport sideways to adjoining columns, and (at least in the
Gregory-Rowntree scheme used in the UM) is considered independent of any
resolved scale vertical air motions. This necessitates the view of
compensating subsidence within the column, whereas some conceptual
models of tropical convection would have the bulk of the ascent in the
convective cores and the associated descent thousands of miles away in
the downward branch of the Hadley circulation. The parametrization
schemes traditionally overlook the existence of condensate in the model
column. The non-PC2 version of the mass-flux convection scheme used in
the UM would have the same large-scale liquid and ice prognostics before
and after convection occurs (apart from a bolt-on evaporation below
convective cloud base), with no regard at all to what happens to it or
its effect on the rest of the convection. Within PC2 we have had to work
to more fully incorporate the condensate into the convection scheme.

Introduction to the convective mass flux scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Within the mass flux scheme the net change in :math:`\overline{q_{cl}}`
and :math:`C_l` etc. comes from two distinct sources. Firstly, the
condensate and cloud fraction injected from the plume (the :math:`Q4`
terms, section `3.5 <#sec:inhomog>`__); secondly, the condensation
response to the vapour and heat changes associated with the detrainment
and compensating subsidence. Strictly, we will see that the :math:`Q4`
terms also include the contribution to the condensate transport by the
compensating subsidence - this casts doubt on the validity of the
application of the injection forcing scenario to calculate the
equivalent cloud fraction change, since ideally the cloud fractions
ought to be transported by the compensating subsidence in a similar way
to the condensate transport (which is documented below).

We therefore split the convective contribution in
(`[eq:dqcldt_and_dcdt] <#eq:dqcldt_and_dcdt>`__) into two parts:

.. math::

   \frac{\partial \overline{q_{cl}}}{\partial t} |_{convection} =
   Q4_l + Q_{environment}
   \label{eq:inhomg_plus_homog}

where :math:`Q_{environment}` is the condensation associated with
changes in the vapour and temperature from the detrainment and
compensating subsidence. Similar splits are made for the cloud
variables, where the injection forcing, section `3.5 <#sec:inhomog>`__,
is used to calculate the first term from :math:`Q4_l`. Section
`4.7.3 <#subsect:q4calculation>`__ looks at the issue of the calculation
of :math:`Q4_l` etc., and section `4.7.4 <#sec:conv_homog>`__ looks at
the calculation of :math:`Q_{environment}`, and its associated cloud
fraction change. We first look at the basic transport equations in a
mass flux convection scheme.

.. _`subsect:basmaseqs`:

Basic Equations for a Convective Mass Flux Scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first consider a generic mass-flux scheme before its application to
PC2. As discussed by Grant and Stirling (personal communication), the
equations for convective tendencies are most simply applied to a
variable, :math:`{\chi}`, that is conserved under moist adiabatic
processes (e.g. total water content). In this case,

.. math::

   {\ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv}} = 
   - \frac{1}{\overline{\rho}} \, \ensuremath{\frac{\partial \, \overline{\rho w^{'} \ensuremath{{\chi}_{\rm{ }}^{\rm{E'}}}}}{\partial \, z}} 
   \label{eq:chibasic}

To parametrize `[eq:chibasic] <#eq:chibasic>`__, the current UM
convection scheme takes a mass flux approximation

.. math::

   \left({\overline{\rho w^{'} \ensuremath{{\chi}_{\rm{ }}^{\rm{E'}}}}} \right)_{\rm{conv}} = M^{\rm{P}} \, 
   \left({ \ensuremath{{\chi}_{\rm{ }}^{\rm{P}}} - \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}} } \right)
   \label{eq:massflux}

which can be differentiated to give

.. math::

   - \frac{1}{\overline{\rho}} \, \ensuremath{\frac{\partial \, \overline{\rho w^{'} \ensuremath{{\chi}_{\rm{ }}^{\rm{E'}}}}}{\partial \, z}} =
   \ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{P}}} \, M^{\rm{P}}}{\partial \, p}} - 
   \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}} \,  \ensuremath{\frac{\partial \, M^{\rm{P}}}{\partial \, p}} -   
   M^{\rm{P}} \, \ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}}}{\partial \, p}}  
   \label{eq:eddyflux}

The bulk cloud model plume equations for mass and :math:`{\chi}` are:

.. math::

   \begin{aligned}
   - \ensuremath{\frac{\partial \, M^{\rm{P}}}{\partial \, p}} & = & 
   \left({ \varepsilon \, M^{\rm{P}} - \mu \, M^{\rm{P}} - \delta \, M^{\rm{P}} } \right)
   \label{eq:dbydpmassflux}  \\
   - \ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{P}}} \, M^{\rm{P}}}{\partial \, p}} & = & \left({ 
   \varepsilon \, M^{\rm{P}} \, \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}}
   - \mu \, M^{\rm{P}} \, \ensuremath{{\chi}_{\rm{ }}^{\rm{R}}} - \delta \, M^{\rm{P}} \, \ensuremath{{\chi}_{\rm{ }}^{\rm{P}}} 
   } \right)\label{eq:dbydpmfchi}  
   \end{aligned}

Equations `[eq:eddyflux] <#eq:eddyflux>`__,
`[eq:dbydpmassflux] <#eq:dbydpmassflux>`__ and
`[eq:dbydpmfchi] <#eq:dbydpmfchi>`__ can then be substituted into
`[eq:chibasic] <#eq:chibasic>`__ to give:

.. math::

   {\ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv}} = 
   - M^{\rm{P}} \, \ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}}}{\partial \, p}} 
   + \mu \, M^{\rm{P}} \, \left({ \ensuremath{{\chi}_{\rm{ }}^{\rm{R}}} - \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}} } \right)
   + \delta \, M^{\rm{P}} \, \left({ \ensuremath{{\chi}_{\rm{ }}^{\rm{P}}} - \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}} } \right)
   \label{eq:chimassflux}

while :math:`{\chi}_{\rm{}}^{\rm{P}}` is obtained from the vertical
gradient derived by combining `[eq:dbydpmassflux] <#eq:dbydpmassflux>`__
and `[eq:dbydpmfchi] <#eq:dbydpmfchi>`__ :

.. math::

   M^{\rm{P}} \, \ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{P}}}}{\partial \, p}} =  
   \varepsilon \, M^{\rm{P}} \, \left({ \ensuremath{{\chi}_{\rm{ }}^{\rm{P}}} - \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}} } \right)- 
   \mu         \, M^{\rm{P}} \, \left({ \ensuremath{{\chi}_{\rm{ }}^{\rm{P}}} - \ensuremath{{\chi}_{\rm{ }}^{\rm{R}}} } \right)
   \label{eq:gradchipar}

Within the model, eqn `[eq:chimassflux] <#eq:chimassflux>`__ would take
a discretized form which actually depends upon whether the model level,
k, is above or at the lowest cloud level (k = cb). Note that the formal
cloud base lies at the half-level below, i.e. on the layer boundary
which is also the top of the turbulent mixed boundary layer. A simple
discretized form of `[eq:chimassflux] <#eq:chimassflux>`__, setting
:math:`{ \mu = 0 }`, is:

.. math::

   \begin{aligned}
   {\ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv, \, k}} & = & m_{\rm{k+1/2}} \, 
   \frac{ \left({\ensuremath{{\chi}_{\rm{k+1}}^{\rm{E}}} - \ensuremath{{\chi}_{\rm{k}}^{\rm{E}}}} \right)}
   {{\Delta z}_{\rm{k \, \rightarrow \, k+1}}} 
   + {\delta}_{\rm{k}} \, m_{\rm{k}} \, \left({ \ensuremath{{\chi}_{\rm{k}}^{\rm{P}}} - \ensuremath{{\chi}_{\rm{k}}^{\rm{E}}} } \right)
   \qquad \ldots \; \mbox{for k $>$ cb} \label{eq:chidisck} \\
   {\ensuremath{\frac{\partial \, \ensuremath{{\chi}_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv, \, cb}} & = & m_{\rm{cb+1/2}} \, 
   \frac{ \left({\ensuremath{{\chi}_{\rm{cb+1}}^{\rm{E}}} - \ensuremath{{\chi}_{\rm{cb}}^{\rm{E}}}} \right)}
   {{\Delta z}_{\rm{cb \, \rightarrow \, cb+1}}} 
   - m_{\rm{cb}} \, 
   \left({ \ensuremath{{\chi}_{\rm{i,cb}}^{\rm{P}}} - \ensuremath{{\chi}_{\rm{cb}}^{\rm{E}}} } \right)\label{eq:chidisccb}
   \end{aligned}

where the initial parcel value :math:`{\chi}_{\rm{i,cb}}^{\rm{P}}` may
be chosen to produce a fixed increment or place a closure condition on
the cloud base flux. In fact, the convection equations (see ) differ
from `[eq:chidisck] <#eq:chidisck>`__ and
`[eq:chidisccb] <#eq:chidisccb>`__ because a different discretization is
used, but the principle is unaltered.

The model convection variables are NOT conserved under moist adiabatic
processes because precipitation processes deplete the column moisture
and condensation processes affect the temperature, specific humidity and
cloud condensate variables. Surprisingly, however, the form of
eqn `[eq:chimassflux] <#eq:chimassflux>`__ is retained even though the
basic equation `[eq:chibasic] <#eq:chibasic>`__ acquires additional
terms for temperature and specific humidity:

.. math::

   \begin{aligned}
   {\ensuremath{\frac{\partial \, \ensuremath{T_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv}} = Q1 & \equiv &
   \left({ \frac{L}{c_{P}} } \right)\, {\overline{Q}}_{\rm{par}}
   - \frac{1}{\overline{\rho}} \, \ensuremath{\frac{\partial \, \overline{\rho w^{'} \ensuremath{T_{\rm{ }}^{\rm{E'}}}}}{\partial \, z}} 
   \label{eq:defineq1} \\ 
   {\ensuremath{\frac{\partial \, \ensuremath{q_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv}} = Q2 & \equiv & - {\overline{Q}}_{\rm{par}}
   - \frac{1}{\overline{\rho}} \, \ensuremath{\frac{\partial \, \overline{\rho w^{'} \ensuremath{q_{\rm{ }}^{\rm{E'}}}}}{\partial \, z}} 
   \label{eq:defineq2} 
   \end{aligned}

where :math:`{\overline{Q}}_{\rm{par}}` is the rate of condensation
which occurs in the ascending plumes.

The reason that `[eq:defineq1] <#eq:defineq1>`__ and
`[eq:defineq2] <#eq:defineq2>`__ retain this form is due to cancellation
from the bulk cloud terms equivalent to
`[eq:dbydpmfchi] <#eq:dbydpmfchi>`__ which are modified in the same way
as `[eq:defineq1] <#eq:defineq1>`__ and
`[eq:defineq2] <#eq:defineq2>`__. The change is seen in the vertical
gradient equations based upon `[eq:gradchipar] <#eq:gradchipar>`__

.. math::

   \begin{aligned}
   M^{\rm{P}} \, \ensuremath{\frac{\partial \, \ensuremath{T_{\rm{ }}^{\rm{P}}}}{\partial \, p}} & = & 
   \varepsilon \, M^{\rm{P}} \, \left({ \ensuremath{T_{\rm{ }}^{\rm{P}}} - \ensuremath{T_{\rm{ }}^{\rm{E}}} } \right)- 
   \mu         \, M^{\rm{P}} \, \left({ \ensuremath{T_{\rm{ }}^{\rm{P}}} - \ensuremath{T_{\rm{ }}^{\rm{R}}} } \right)-
   \left({ \frac{L}{c_{P}} } \right)\, {\overline{Q}}_{\rm{par}}  \label{eq:gradtpar} \\
   M^{\rm{P}} \, \ensuremath{\frac{\partial \, \ensuremath{q_{\rm{ }}^{\rm{P}}}}{\partial \, p}} & = & 
   \varepsilon \, M^{\rm{P}} \, \left({ \ensuremath{q_{\rm{ }}^{\rm{P}}} - \ensuremath{q_{\rm{ }}^{\rm{E}}} } \right)- 
   \mu         \, M^{\rm{P}} \, \left({ \ensuremath{q_{\rm{ }}^{\rm{P}}} - \ensuremath{q_{\rm{ }}^{\rm{R}}} } \right)+
   {\overline{Q}}_{\rm{par}}  \label{eq:gradqpar}  \\
   M^{\rm{P}} \, \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{ }}^{\rm{P}}}}{\partial \, p}} & = & 
   \varepsilon \, M^{\rm{P}} \, \left({ \ensuremath{l_{\rm{ }}^{\rm{P}}} - \ensuremath{l_{\rm{ }}^{\rm{E}}} } \right)
   - {\overline{Q}}_{\rm{par}} + PPN  \label{eq:gradlpar} 
   \end{aligned}

The final calculation of rates in the current condensation scheme (,
section 10) assumes a further condensation term,
:math:`{\overline{Q}}_{\rm{reset}}`, which acts to make the net rate of
change of condensate equal zero, and a final assumption is made that the
environment values of condensate remain zero (and also that
:math:`l_{\rm{ }}^{\rm{R}}` = :math:`l_{\rm{ }}^{\rm{P}}`). The result
is basic equations

.. math::

   \begin{aligned}
   {\ensuremath{\frac{\partial \, \ensuremath{T_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv}} & = & Q1 -
   \left({ \frac{L}{c_{P}} } \right)\, {\overline{Q}}_{\rm{reset}}
   \label{eq:basictold} \\ 
   {\ensuremath{\frac{\partial \, \ensuremath{q_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv}} & = & Q2 + {\overline{Q}}_{\rm{reset}}
   \label{eq:basicqold} \\
   0 \equiv {\ensuremath{\frac{\partial \, \ensuremath{l_{\rm{ }}^{\rm{E}}}}{\partial \, t}}}_{\rm{conv}} & = & {\overline{Q}}_{\rm{par}} - 
   {\overline{Q}}_{\rm{reset}} - PPN
   - \frac{1}{\overline{\rho}} \, \ensuremath{\frac{\partial \, \overline{\rho w^{'} \ensuremath{l_{\rm{ }}^{\rm{E'}}}}}{\partial \, z}} \nonumber \\
   & = & 
   \mu \, M^{\rm{P}} \, \ensuremath{l_{\rm{ }}^{\rm{P}}} + \delta \, M^{\rm{P}} \, \ensuremath{l_{\rm{ }}^{\rm{P}}} -
   {\overline{Q}}_{\rm{reset}} 
   \label{eq:basiclold} 
   \end{aligned}

By analogy with equations `[eq:defineq1] <#eq:defineq1>`__ and
`[eq:defineq2] <#eq:defineq2>`__, we can define a :math:`Q4` from
`[eq:basiclold] <#eq:basiclold>`__ and state that for the control
convection scheme :math:`Q4 = 0`. The PC2 scheme requires a reassessment
of these assumptions because we wish to allow non-zero environment
condensate values and to allow them to change.

.. _`subsect:q4calculation`:

Calculation of Grid-Box Averaged Condensate Rate (Q4)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The PC2 condensation scheme allows convection to feed cloud condensate
(ice or liquid) directly into the large scale and to update the cloud
amount accordingly.

Define

.. math::

   \begin{aligned}
   \left({ \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{l}}^{\rm{ }}}}{\partial \, t}} } \right)_{\rm{conv}} = Q4_{\rm{l}} & \equiv &  
   {\overline{Q}}_{\rm{l, par}} - {\overline{Q}}_{\rm{l, reset}} - RAIN - 
   \frac{1}{\overline{\rho}} \, \ensuremath{\frac{\partial \, \overline{\rho w^{'} \ensuremath{l_{\rm{l}}^{\rm{'}}}}}{\partial \, z}}
   \label{eq:defineq4l} \\
   \left({ \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{f}}^{\rm{ }}}}{\partial \, t}} } \right)_{\rm{conv}} = Q4_{\rm{f}} & \equiv &  
   {\overline{Q}}_{\rm{f, par}} - {\overline{Q}}_{\rm{f, reset}} - SNOW - 
   \frac{1}{\overline{\rho}} \, \ensuremath{\frac{\partial \, \overline{\rho w^{'} \ensuremath{l_{\rm{f}}^{\rm{'}}}}}{\partial \, z}}
   \label{eq:defineq4f}
   \end{aligned}

where the PC2 assumption thus far has been that
:math:`{\overline{Q}}_{\rm{l, reset}} = 0
= {\overline{Q}}_{\rm{f, reset}}`.

- The current convection scheme assumes that parcel condensate is single
  phase (ie. either all liquid or all frozen) and this is seriously
  hard-wired into the code. Thus we can treat the precipitation and
  parcel condensation processes in :math:`Q4_{\rm{l}}` and
  :math:`Q4_{\rm{f}}` separately without worrying about cross-transfer
  between the two because at most only one set will ever be active in a
  given grid box at one time. However, even for the inactive (zero
  parcel condensate) phase, convection mixes environmental air into the
  parcel and can therefore maintain a non-zero :math:`Q4`. Enablement of
  multiple phase condensate in the current scheme is a task requiring
  great caution as the formulations are extremely sensitive to errors in
  assignment of condensate phase.

Based on `[eq:gradlpar] <#eq:gradlpar>`__, the vertical dependence of
condensate is calculated as

.. math::

   \begin{aligned}
   \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{l}}^{\rm{P}}}}{\partial \, p}} & = & \varepsilon \, 
   \left({ \ensuremath{l_{\rm{l}}^{\rm{P}}} - \ensuremath{l_{\rm{l}}^{\rm{E}}} } \right)- 
   \frac{{\overline{Q}}_{\rm{l, par}}}{M^{\rm{P}}} - 
   \frac{RAIN}{M^{\rm{P}}}  \label{eq:vertparl} \\
   \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{f}}^{\rm{P}}}}{\partial \, p}} & = & \varepsilon \,  
   \left({ \ensuremath{l_{\rm{f}}^{\rm{P}}} - \ensuremath{l_{\rm{f}}^{\rm{E}}} } \right)- 
   \frac{{\overline{Q}}_{\rm{f, par}}}{M^{\rm{P}}} - 
   \frac{SNOW}{M^{\rm{P}}}\label{eq:vertparf} 
   \end{aligned}

Following , equations `[eq:dbydpmassflux] <#eq:dbydpmassflux>`__,
`[eq:vertparl] <#eq:vertparl>`__ and `[eq:vertparf] <#eq:vertparf>`__
are discretized:

.. math::

   \begin{aligned}
   M_{\rm{k} + 1} & = & M_{\rm{k}} \, 
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)\,
   \left({ 1 - \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)\,
   EPSS_{\rm{k}}  \label{eq:discdmfbydp} \\
   \ensuremath{l_{\rm{l \, k + 1}}^{\rm{P}}} & = & \left({ 
   \ensuremath{l_{\rm{l \, k}}^{\rm{P}}} +  
   \varepsilon_{\rm{k} + 1/4} \, \Delta p_{\rm{k} + 1/4} \, \ensuremath{l_{\rm{l \, k}}^{\rm{E}}} +
   \varepsilon_{\rm{k} + 3/4} \, \Delta p_{\rm{k} + 3/4} \, 
   \left[{1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4}} \right]\,
   \ensuremath{l_{\rm{l \, k + 1}}^{\rm{E}}}
   } \right)\, / \, \left({EPSS_{\rm{k}}} \right)\nonumber \\
   { } & { } & + \left({ {\overline{Q}}_{\rm{l} \, \rm{k} + 1} \, / \, M_{\rm{k} + 1}} \right)
   - \left({ RAIN_{\rm{k} + 1} \, / \, M_{\rm{k} + 1} } \right)
   \label{eq:discvparl} \\
   \ensuremath{l_{\rm{f \, k + 1}}^{\rm{P}}} & = & \left({ 
   \ensuremath{l_{\rm{f \, k}}^{\rm{P}}} + 
   \varepsilon_{\rm{k} + 1/4} \, \Delta p_{\rm{k} + 1/4} \, \ensuremath{l_{\rm{f \, k}}^{\rm{E}}} +
   \varepsilon_{\rm{k} + 3/4} \, \Delta p_{\rm{k} + 3/4} \, 
   \left[{1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4}} \right]\,
   \ensuremath{l_{\rm{f \, k + 1}}^{\rm{E}}}
   } \right)\, / \, \left({EPSS_{\rm{k}}} \right)\nonumber \\
   { } & { } & + \left({ {\overline{Q}}_{\rm{f} \, \rm{k} + 1} \, / \, M_{\rm{k} + 1}} \right)
   - \left({ SNOW_{\rm{k} + 1} \, / \, M_{\rm{k} + 1} } \right)
   \label{eq:discvparf}
   \end{aligned}

where :math:`EPSS_{\rm{k}} = 
\left({1 + \varepsilon_{\rm{k} + 3 / 4} \, \Delta p_{\rm{k} + 3 / 4}} \right)\, 
\left({1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4}} \right)`.

The condensation and precipitation terms in equations
`[eq:discdmfbydp] <#eq:discdmfbydp>`__,
`[eq:discvparl] <#eq:discvparl>`__ and
`[eq:discvparf] <#eq:discvparf>`__ make the equations implicit. They are
therefore solved by starting with an ascent in which condensation and
precipitation terms are suppressed:

.. math::

   \begin{aligned}
   \ensuremath{l_{\rm{l \, k + 1}}^{\rm{P}}} & = & \frac{\left({ 
   \ensuremath{l_{\rm{l \, k}}^{\rm{P}}} +  
   \varepsilon_{\rm{k} + 1/4} \, \Delta p_{\rm{k} + 1/4} \, \ensuremath{l_{\rm{l \, k}}^{\rm{E}}} +
   \varepsilon_{\rm{k} + 3/4} \, \Delta p_{\rm{k} + 3/4} \, 
   \left[{1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4}} \right]\,
   \ensuremath{l_{\rm{l \, k + 1}}^{\rm{E}}}
   } \right)}{EPSS_{\rm{k}}}  \label{eq:discvparldry} \\
   \ensuremath{l_{\rm{f \, k + 1}}^{\rm{P}}} & = & \frac{\left({ 
   \ensuremath{l_{\rm{f \, k}}^{\rm{P}}} + 
   \varepsilon_{\rm{k} + 1/4} \, \Delta p_{\rm{k} + 1/4} \, \ensuremath{l_{\rm{f \, k}}^{\rm{E}}} +
   \varepsilon_{\rm{k} + 3/4} \, \Delta p_{\rm{k} + 3/4} \, 
   \left[{1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4}} \right]\,
   \ensuremath{l_{\rm{f \, k + 1}}^{\rm{E}}}
   } \right)}{EPSS_{\rm{k}}}  \label{eq:discvparfdry}
   \end{aligned}

At the base of the convective plume (ie. the level immediately above
cloud base), :math:`l_{\rm{l \, k}}^{\rm{P}}` is initialized to
:math:`l_{\rm{l \, i}}^{\rm{P}}` and :math:`l_{\rm{f \, k}}^{\rm{P}}` to
:math:`l_{\rm{f \, i}}^{\rm{P}}`, where the initial values are chosen
such that the modified form of `[eq:chidisccb] <#eq:chidisccb>`__
produces zero fluxes at cloud base:

.. math::

   \begin{aligned}
   Q4_{\rm{l}}(cb) = 0 & = & M_{\rm{cb+1/2}}^{\rm{P}} \, 
   \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{l}}^{\rm{E}}}}{\partial \, p}}  - M_{\rm{cb}}^{\rm{P}}\, 
   \left({ \ensuremath{l_{\rm{l}}^{\rm{P \, i}}} - \ensuremath{l_{\rm{l}}^{\rm{E}}}(\rm{cb}) } \right)\label{eq:q4lcbi} \\ 
   Q4_{\rm{f}}(cb) = 0 & = & M_{\rm{cb+1/2}}^{\rm{P}} \, 
   \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{f}}^{\rm{E}}}}{\partial \, p}}  - M_{\rm{cb}}^{\rm{P}}\, 
   \left({ \ensuremath{l_{\rm{f}}^{\rm{P \, i}}} - \ensuremath{l_{\rm{f}}^{\rm{E}}}(\rm{cb}) } \right)\label{eq:q4fcbi} 
   \end{aligned}

As the convection scheme makes the single phase assumption for parcel
condensate, it may be necessary to melt or freeze entrained condensate
at this point and adjust the temperature accordingly.

.. math::

   \begin{aligned}
   \theta_{\rm{k + 1}}^{\rm{P}} = \theta_{\rm{k + 1}}^{\rm{P}} - 
   \left(\frac{L_{\rm{F}}}{C_{p} \, \Pi_{\rm{k + 1}}} \right)\, \ensuremath{l_{\rm{f \, k + 1}}^{\rm{P}}} 
   & \; \ldots \; & \mbox{ if \ensuremath{l_{\rm{f \, k + 1}}^{\rm{P}}} is melted }  
   \label{eqn:meltlf} \\
   \theta_{\rm{k + 1}}^{\rm{P}} = \theta_{\rm{k + 1}}^{\rm{P}} + 
   \left(\frac{L_{\rm{F}}}{C_{p} \, \Pi_{\rm{k + 1}}} \right)\, \ensuremath{l_{\rm{l \, k + 1}}^{\rm{P}}} 
   & \; \ldots \; & \mbox{ if \ensuremath{l_{\rm{l \, k + 1}}^{\rm{P}}} is frozen }  
   \label{eqn:freezell} 
   \end{aligned}

Once a final value for the condensation term
:math:`{\overline{Q}}_{\rm{x} \, \rm{k} + 1} \, / \, M_{\rm{k} + 1}` has
been calculated from the parcel specific humidity equations, it can then
be added to the parcel condensate to give a final pre-precipitation
value.

- In practice, the rates :math:`{\overline{Q}}_{\rm{x} \, \rm{k} + 1}`
  and :math:`PPN` are not calculated explicitly in the code. Instead,
  their effect is applied directly as increments to the temperature and
  moisture fields.

The precipitation calculation is unaltered.

.. math::

   P_{\rm{k} + 1} = \left({ \ensuremath{l_{\rm{k + 1}}^{\rm{P}}} - \ensuremath{l_{\rm{MIN}}^{\rm{P}}} } \right)\, 
   M_{\rm{k} + 1} \, / \, g
   \label{eq:precip}

where :math:`l_{\rm{k + 1}}^{\rm{P}}` =
:math:`l_{\rm{l \, k + 1}}^{\rm{P}}` +
:math:`l_{\rm{f \, k + 1}}^{\rm{P}}`.

- Actually, given that the precipitation calculation appears to be based
  upon the hydrostatic equation, it is debatable whether it is even
  suitable for use with the New Dynamics model and I guess therefore
  that this needs revisiting at some point.

This reduces the parcel condensate to :

.. math::

   \begin{aligned}
   \ensuremath{l_{\rm{l \, k + 1}}^{\rm{P}}} & = & \left({ 
   \frac{\ensuremath{l_{\rm{l \, k + 1}}^{\rm{P}}}}{\ensuremath{l_{\rm{k + 1}}^{\rm{P}}}}
   } \right)\, \ensuremath{l_{\rm{MIN}}^{\rm{P}}}  \label{eq:vparlfinal} \\
   \ensuremath{l_{\rm{f \, k + 1}}^{\rm{P}}} & = & \left({
   \frac{\ensuremath{l_{\rm{f \, k + 1}}^{\rm{P}}}}{\ensuremath{l_{\rm{k + 1}}^{\rm{P}}}}
   } \right)\, \ensuremath{l_{\rm{MIN}}^{\rm{P}}}  \label{eq:vparffinal} 
   \end{aligned}

The final parcel condensate values are then used in the rate calculation
based upon eqn `[eq:basiclold] <#eq:basiclold>`__:

.. math::

   \begin{aligned}
   Q4_{\rm{l}}(k) & = &  M_{\rm{k+1/2}}^{\rm{P}} \, \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{l}}^{\rm{E}}}}{\partial \, p}}   + 
   \left({ {\mu}_{\rm{k}} \, M_{\rm{k}}^{\rm{P}} + 
   {\delta}_{\rm{k}} \, M_{\rm{k}}^{\rm{P}} } \right)\, 
   \left({ \ensuremath{l_{\rm{l}}^{\rm{P}}}(\rm{k}) - \ensuremath{l_{\rm{l}}^{\rm{E}}}(\rm{k}) } \right)-
   {\overline{Q}}_{\rm{l, reset}}  \label{eq:q4lmassf} \\
   Q4_{\rm{f}}(k) & = &  M_{\rm{k+1/2}}^{\rm{P}} \, \ensuremath{\frac{\partial \, \ensuremath{l_{\rm{f}}^{\rm{E}}}}{\partial \, p}}  +
   \left({ {\mu}_{\rm{k}} \, M_{\rm{k}}^{\rm{P}} + 
   {\delta}_{\rm{k}} \, M_{\rm{k}}^{\rm{P}} } \right)\, 
   \left({ \ensuremath{l_{\rm{f}}^{\rm{P}}}(\rm{k}) - \ensuremath{l_{\rm{f}}^{\rm{E}}}(\rm{k}) } \right)- 
   {\overline{Q}}_{\rm{f, reset}}  \label{eq:q4fmassf}
   \end{aligned}

Note that, as a side-effect, the environment equations for potential
temperature and specific humidity are also altered because the
condensate is no longer re-evaporated at the end
(:math:`{\overline{Q}}_{\rm{l, reset}} = 0
= {\overline{Q}}_{\rm{f, reset}}`):

.. math::

   \begin{aligned}
   \frac{\Delta \, \theta_{\rm{k}}^{\rm{E}}}{\Delta \, t} = 
   \left(\frac{ M_{\rm{k}} }{ \Delta \, p_{\rm{k}} } \right)
   \left[{
   \left({ 1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4} } \right)
   \left({ 1 - \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \theta_{\rm{k + 1}}^{\rm{E}} - \theta_{\rm{k}}^{\rm{E}} } \right)
   } \right . & + &  \nonumber \\
   \left({ \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \theta_{\rm{k}}^{\rm{R}} - \theta_{\rm{k}}^{\rm{E}} } \right)
   & + &  \nonumber \\
   \left . {
   \left({ \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \theta_{\rm{k}}^{\rm{P}} - \theta_{\rm{k}}^{\rm{E}} } \right)
   } \right]& { } &  \label{eq:enviroth}
   \end{aligned}

and

.. math::

   \begin{aligned}
   \frac{\Delta \, q_{\rm{k}}^{\rm{E}}}{\Delta \, t} = 
   \left(\frac{ M_{\rm{k}} }{ \Delta \, p_{\rm{k}} } \right)
   \left[{
   \left({ 1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4} } \right)
   \left({ 1 - \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ q_{\rm{k + 1}}^{\rm{E}} - q_{\rm{k}}^{\rm{E}} } \right)
   } \right . & + &  \nonumber \\
   \left({ \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ q_{\rm{k}}^{\rm{R}} - q_{\rm{k}}^{\rm{E}} } \right)
   & + &  \nonumber \\
   \left . {
   \left({ \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ q_{\rm{k}}^{\rm{P}} - q_{\rm{k}}^{\rm{E}} } \right)
   } \right]& { } &  \label{eq:enviroq}
   \end{aligned}

Similarly, eqns `[eq:q4lmassf] <#eq:q4lmassf>`__ and
`[eq:q4fmassf] <#eq:q4fmassf>`__ have a discretized form as follows:

.. math::

   \begin{aligned}
   \frac{\Delta \, \ensuremath{l_{\rm{l \, k}}^{\rm{E}}}}{\Delta \, t} = 
   \left(\frac{ M_{\rm{k}} }{ \Delta \, p_{\rm{k}} } \right)
   \left[{
   \left({ 1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4} } \right)
   \left({ 1 - \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \ensuremath{l_{\rm{l \, k + 1}}^{\rm{E}}} - \ensuremath{l_{\rm{l \, k}}^{\rm{E}}} } \right)
   } \right . & + &  \nonumber \\
   \left({ \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \ensuremath{l_{\rm{l \, k}}^{\rm{P}}} - \ensuremath{l_{\rm{l \, k}}^{\rm{E}}} } \right)
   & + &  \nonumber \\
   \left . {
   \left({ \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \ensuremath{l_{\rm{l \, k}}^{\rm{P}}} - \ensuremath{l_{\rm{l \, k}}^{\rm{E}}} } \right)
   } \right]& { } &  \label{eq:enviroll}
   \end{aligned}

and

.. math::

   \begin{aligned}
   \frac{\Delta \, \ensuremath{l_{\rm{f \, k}}^{\rm{E}}}}{\Delta \, t} = 
   \left(\frac{ M_{\rm{k}} }{ \Delta \, p_{\rm{k}} } \right)
   \left[{
   \left({ 1 + \varepsilon_{\rm{k} + 1 / 4} \, \Delta p_{\rm{k} + 1 / 4} } \right)
   \left({ 1 - \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \ensuremath{l_{\rm{f \, k + 1}}^{\rm{E}}} - \ensuremath{l_{\rm{f \, k}}^{\rm{E}}} } \right)
   } \right . & + &  \nonumber \\
   \left({ \delta_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ 1 - \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \ensuremath{l_{\rm{f \, k}}^{\rm{P}}} - \ensuremath{l_{\rm{f \, k}}^{\rm{E}}} } \right)
   & + &  \nonumber \\
   \left . {
   \left({ \mu_{\rm{k}} \, \Delta p_{\rm{k} + 1 / 2} } \right)
   \left({ \ensuremath{l_{\rm{f \, k }}^{\rm{P}}} - \ensuremath{l_{\rm{f \, k}}^{\rm{E}}} } \right)
   } \right]& { } &  \label{eq:envirolf}
   \end{aligned}

.. _`sec:conv_homog`:

Background condensation
~~~~~~~~~~~~~~~~~~~~~~~

The modification to the convective plume will result in the transport,
detrainment and entrainment of condensate, in addition to the transport
of vapour and heat. Although condensation processes within the plume are
treated, it does not treat condensation in the environment, which is
forced by the compensating subsidence. We wish to relate the
environmental increments of vapour and temperature to a forcing that can
be applied in the environment. Because we know that any detrained air
associated with detrained liquid water from the plume must be saturated
with respect to liquid water, we are able to translate the environmental
changes into forcings.

Here we will consider that the vapour change in the gridbox is as a
result of *saturated with respect to liquid water* air being injected
from the plume and background air being displaced. We do not consider
whether the background air is at saturation yet, for we wish to derive
the expression for the required condensation if this is not the case. We
consider only liquid water clouds, ice clouds have no background
condensation applied as we do not make the instantaneous condensation
assumption.

Hence we can write

.. math::

   \Delta \overline{q} = \Delta C_S ( q_{sat liq}(\overline{T_{s}}) - \overline{q} )
   + (1 - \Delta C_S) \Delta \overline{q_{background}}

where :math:`\Delta C_S` is the volume of plume air that is detrained
into the gridbox, as discussed by :raw-latex:`\cite{bwg03}`. :math:`T_s`
is the temperature of the air injected into the gridbox by the plume.
The first term is simply the difference between the value of :math:`q`
in the plume and what was previously in the gridbox, and the second term
is the effect of a background change of :math:`q` that will be applied
across the part of the gridbox that is not associated with the injected
air. We write this as:

.. math::

   (1 - \Delta C_S) \Delta \overline{q_{background}} = \Delta \overline{q} - 
   \Delta C_S ( q_{sat liq}(\overline{T_{s}}) - \overline{q} ) .
   \label{eqn:1mcs}

Now we recognise that

.. math::

   \Delta \overline{q} = Q2~ \Delta t
   \label{eqn:Q2}

where :math:`Q2` is the rate of moistening of the whole gridbox due to
convection. Remember that, at this stage, we haven’t done any
condensation outside of the plume. Hence to calculate the condensation
we should apply the background change in :math:`\overline{q}` as a
uniform forcing for the background air. Hence
(`[eqn:1mcs] <#eqn:1mcs>`__) becomes, using (`[eqn:Q2] <#eqn:Q2>`__),

.. math::

   (1 - \Delta C_S) A_q |_{background} \Delta t =  Q2 ~ \Delta t - \Delta C_S
   ( q_{sat liq}(T_{s}) - \overline{q} ) .
   \label{eqn:Aq}

where :math:`A_q |_{background}` is the currently unknown background
forcing of :math:`q` (see :raw-latex:`\cite{gwb02}`) and
:math:`\Delta t` is the timestep. We can do the same analysis for the
temperature change, and obtain

.. math::

   (1 - \Delta C_S) A_T |_{background} \Delta t  =  Q1~ \Delta t -
   \Delta C_S (T_s - \overline{T} )
   \label{eqn:AT}

where Q1 is the rate of warming in the gridbox due to convection and
:math:`A_T |_{background}` is the currently unknown background forcing
of temperature.

The full change of liquid water content in the gridbox is that injected,
:math:`Q4~\Delta t`, plus the amount of condensation in the background
from the uniform forcings (see :raw-latex:`\cite{wg03}`). Note that the
uniform forcings are only applied across a proportion
:math:`1 - \Delta C_S` of the gridbox. Hence these two terms give, using
the homogeneous forcing equations (`[dqcldt] <#dqcldt>`__) and
(`[eq:deltaqc_exp2] <#eq:deltaqc_exp2>`__),

.. math::

   \Delta \overline{q_{cl}} |_{convection} = Q4 \Delta t
   + (1 - \Delta C_S) a_L C_l (A_q |_{background} \Delta t
   - \alpha A_T |_{background} \Delta t ).  
   \label{eqn:qclconv}

Using (`[eqn:Aq] <#eqn:Aq>`__) and (`[eqn:AT] <#eqn:AT>`__) to expand
the forcing terms in (`[eqn:qclconv] <#eqn:qclconv>`__) gives

.. math::

   \Delta \overline{q_{cl}} |_{convection} = Q4 \Delta t
   + \Delta t a_L C_l ( Q2 - \alpha Q1) - \Delta C_S a_L C_l 
   (q_{sat liq}(T_s)-\overline{q} - \alpha (T_s - \overline{T})) .

We now note that

.. math:: q_{sat} (T_s) - q_{sat liq} (\overline{T}) = \alpha (T_s - \overline{T} )

and hence the final result

.. math::

   \Delta \overline{q_{cl}} |_{convection} = Q4 \Delta t
   + \Delta t ~ a_L C_l (  ( Q2 - \alpha Q1) - \Delta C_S 
   (q_{sat liq}(\overline{T}) - \overline{q} ) ) .
   \label{eqn:dqcl}

There is thus an extra term,
:math:`-\Delta C_S  (q_{sat}(\overline{T})-\overline{q} )`, which needs
to be included in addition to the standard application of the
homogeneous forcing of :math:`Q1` and :math:`Q2` (this is represented by
the second term of the expression). This has arisen from the requirement
that the vapour injected by the plume is saturated. We need simply to
retrieve the value of :math:`\Delta C_S` to complete the
parametrization. This can be straightforwardly obtained from
(`[eq:dcldt_inhom] <#eq:dcldt_inhom>`__), which links the net change of
liquid cloudy volume due to the injection, :math:`\Delta C_{injection}`,
with :math:`\Delta C_S`.

.. math:: \Delta C_{injection} = (g_l - C_l) \Delta C_S

where :math:`g_l` is 1 if the injected cloud is of liquid phase and 0 if
it is of ice phase. We already know :math:`\Delta C_{injection}` from
the injection forcing arguments (`[eq:dctdt_xl] <#eq:dctdt_xl>`__) above
that link it to :math:`Q4`. We therefore complete the parametrization by
calculating :math:`\Delta C_S` based on whether :math:`\Delta C` is
positive or negative. If :math:`\Delta C` is positive, we assume that
the plume must be of liquid phase and hence

.. math::

   \Delta C_S = \frac{\Delta C_{injection}} {1 -C_l} .
   \label{eqn:cs1}

If :math:`\Delta C_l` is negative, we assume that the plume must be of
ice phase and hence

.. math::

   \Delta C_S = - \frac{\Delta C_{injection}} {C_l} .
   \label{eqn:cs2}

Here we have still assumed that the vapour content in the detrained
plume is equal to :math:`q_{sat liq}`. A better assumption may be to
replace the :math:`q_{sat liq}` term in (`[eqn:dqcl] <#eqn:dqcl>`__)
with a :math:`q_{sat}` expression that depends on the volume fraction of
detrained condensate that is liquid phase, :math:`g_l`.

If :math:`\Delta C_{injection}` is zero, we assume that
:math:`\Delta C_S` is 0 also. Equations
(`[eqn:dqcl] <#eqn:dqcl>`__),(`[eqn:cs1] <#eqn:cs1>`__), and
(`[eqn:cs2] <#eqn:cs2>`__) form the parametrization for
:math:`\Delta \overline{q_{cl}}|_{convection}`. The representation of
:math:`\Delta C_{convection}` is similar in form to
:math:`\Delta \overline{q_{cl}}|_{convection}`:

.. math::

   \Delta C |_{convection} = \Delta C_{injection}
   + \Delta t ~ a_L G(-Q_c) (  ( Q2 - \alpha Q1) - \Delta C_S  (q_{sat}(\overline{T}) 
   - \overline{q} ) ) .

where the specification of :math:`G(-Q_c)` follows
(`[eqn22] <#eqn22>`__). Note that the code includes the numerical limit
restriction that :math:`\Delta C_S` is between 0 and 1.

Thus we are able to parametrize the net condensation and cloud changes
associated with the :math:`Q1` and :math:`Q2` terms in a physically more
consistent way than using simple homogeneous application of these terms.

As an aside, we note that in the :raw-latex:`\cite{t93}` scheme the
condensation and cloud fraction change associated with the compensating
subsidence is taken out of the convection term by adding the vertical
motion associated with the compensating subsidence to the large-scale
vertical velocity before the :raw-latex:`\cite{t93}` equivalent of the
homogeneous forcing term is applied. By doing so it ensures that any
balance between these two terms (as the tropical circulation is commonly
analysed to show) is removed before the net effect is calculated,
leading to more accurate numerical behaviour.

Homogeneous forcing of the environment by convective-subsidence pressure change
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To this end, the code includes an option to perform the homogeneous
forcing of liquid cloud by convection using the “pressure forcing” from
the convective subsidence, consistent with the pressure forcing by
large-scale advection (see sections `4.5 <#sec:advec>`__ and
`4.8 <#sec:pres>`__). This approach replaces the above method of
homogeneous forcing by convection if the UM namelist switch
**l_pc2_homog_conv_pressure** is turned on. By applying the same
homogeneous forcing method for advection and convectively-forced
subsidence, we should get the correct zero net change in liquid cloud in
the common situation where the large-scale ascent and convective
subsidence are in balance (implying no net vertical displacement of
environment parcels).

Under this option, the increments to :math:`\overline{q_{cl}}` and
:math:`C_l` produced by the convection scheme are assumed to already
include the effects of entrainment, detrainment (i.e. injection) and
compensating subsidence (i.e. vertical advection) as expressed by
equation `[eq:chimassflux] <#eq:chimassflux>`__, but exclude the effects
of homogeneous forcing of clouds in the enviroment. Note that taking
equation `[eq:chimassflux] <#eq:chimassflux>`__ with :math:`\chi` set to
water vapour :math:`q`, detrainment of saturated air into a subsaturated
environment will imply a positive tendency of :math:`\overline{q}`, but
this is *not* a homogeneous forcing, since the increase in
:math:`\overline{q}` is entirely due to injecting new parcels of
saturated air without altering the existing environment parcels. Setting
:math:`\chi` to be :math:`\overline{q_{cl}}` or :math:`C_l` in equation
`[eq:chimassflux] <#eq:chimassflux>`__, there is a simply-calculated
source of cloud water and fraction wherever the detrained air is cloudy
(:math:`C_l=1` in the detrained parcel), and we assume these terms have
been calculated this way inside the convection scheme.

Since entrainment and detrainment do not constitute a homgeneous forcing
and are already accounted for in the convection scheme, the only
component of the convective forcing of liquid cloud that needs to be
done by the PC2 call after convection is the homogeneous forcing by the
subsidence term. This is in essence a vertical advection (environmental
forced descent by updrafts, or forced ascent by downdrafts). The
homogeneous forcing can be calculated from the expected pressure change
(and accompaying adiabatic temperature change) following the environment
as it is vertically displaced. Conveniently, the UM already holds the
convective mass-flux in units of Pa s\ :math:`^{-1}`, so it already
expresses the pressure vertical velocity forced by subsidence in the
environment:

.. math::

   \Delta p^E = \Delta t \left( M_{up} - M_{dwn} \right)
   \label{eq:delta_p_conv}

where :math:`M_{up}` is the updraft mass-flux, :math:`M_{dwn}` is the
downdraft mass-flux, and :math:`\Delta t` is the model timestep length.
The adiabatic temperature change following an environment parcel
subsided from pressure :math:`p - \Delta p^E` to :math:`p` is then given
by:

.. math::

   \Delta T^E = \theta^E \left( \left(\frac{p}{p_{ref}}\right)^\kappa
                              - \left(\frac{p - \Delta p^E}{p_{ref}}\right)^\kappa
                         \right)
   \label{eq:delta_t_conv}

where :math:`\theta^E` is the environment potential temperature,
:math:`p_{ref}` is the reference pressure used to define potential
temperature, and :math:`\kappa = \frac{R_d}{c_p}` is the ratio of the
gas constant for dry air over its heat capacity at constant pressure.
`[eq:delta_p_conv] <#eq:delta_p_conv>`__ and
`[eq:delta_t_conv] <#eq:delta_t_conv>`__ are passed into the PC2
homogeneous forcing routine after convection as the forcings to be
applied (with the forcings to all other variables set to zero).

Convective cloud amount
~~~~~~~~~~~~~~~~~~~~~~~

It is a debatable point whether the convective cloud fraction should be
set to zero. Although this was one of the original key concepts of PC2,
the cloud that is detrained from the convection scheme is into the
*environment*, and does not represent the tower cloud. However, it
should be able to represent recently detrained cloudy air in a more
accurate way than by simply appealing to a diagnostic large-scale cloud
scheme. There are similar issues associated with the cloud fraction
predicted from the Tiedtke scheme. Probably the most consistent
interpretation is the inclusion of a tower cloud fraction within PC2,
but not an anvil cloud. However, we need to consider carefully any
double counting (or non-counting) implications. In the PC2:64
formulation, we can represent the large optical depths associated with
new anvils, although we also tend to overestimate the optical depth of
shallow convective clouds. Hence we choose to apply neither a diagnostic
anvil or tower cloud, so similar to Tiedtke, and let the large-scale
cloud fraction represent the convection completely.

Strictly speaking these choices are independent of the PC2 scheme, being
simply choices that are available as part of the existing convection
scheme, but they are clearly directly related to the rest of the cloud
scheme formulation.

CAPE scaling
~~~~~~~~~~~~

The CAPE scaling option in the mass-flux convection scheme scales its
increments by the calculated values of
:math:`\frac{1}{CAPE} \frac{dCAPE}{dt}`. This applies also to all the
PC2 calculated condensate and cloud fraction increments. Additionally,
in order to achieve reasonable mass flux profiles, it has proved
necessary to adjust the calculation of :math:`\frac{dCAPE}{dt}` to use
increments of :math:`\Delta \theta` (potential temperature) and
:math:`\Delta q` calculated using a non-PC2 calculation of these terms.
Hence we consider any detrained condensate to have been evaporated when
we calculate :math:`\frac{dCAPE}{dt}`.

Convective precipitation
~~~~~~~~~~~~~~~~~~~~~~~~

The amount of condensate detrained from convective plumes, and hence the
amount of moisture in the upper levels of the atmosphere, is very
dependent upon the amount of convective precipitation that is allowed to
fall from the column. The standard parametrization of this is that any
condensate greater than a specified value (dependent on :math:`T`) is
precipitated, leaving the rest to be detrained.

PC2 incorporates a tuning to this function of temperature by applying
the additional restriction that the limit may not fall to less than
:math:`2 \times 10^{-4}~kg~kg^{-1}`. This implies a difference at
temperatures less than around :math:`-42 ^{\circ} C`, with the tuning
allowing less precipitation and greater detrainment. This change is
necessary in order to produce thick enough anvil clouds.

.. _`sec:plume_phase`:

Phase of condensate
~~~~~~~~~~~~~~~~~~~

The phase of the convective condensate *carried in the plume* is
determined by a single phase change temperature TICE, with condensate
entirely in the ice phase at colder temperatures and condensate entirely
in the liquid phase at warmer temperatures. For PC2:66, this temperature
is -10 :math:`^{\circ}` C.

.. math::

   \delta_{xl} = \left\{  \begin{array}{ll}
                     1,                             &  T_{plume} \ge -10 ^{\circ} C \\
                     0,                             &  T_{plume} <  -10 ^{\circ} C
                 \end{array} \right.

.. math::

   \delta_{xi} = \left\{  \begin{array}{ll}
                     0,                             &  T_{plume} \ge -10 ^{\circ} C \\
                     1,                             &  T_{plume} <  -10 ^{\circ} C
                 \end{array} \right.

.. _`sec:conv-simpler`:

Tidier way of coupling convection and PC2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This area is still under development. But in brief, work is udner way to
ensure that the convective plume smoothly transitions from detraining
liquid to detraining ice, rather than using the abrupt change implied by
the current formulation of the convection scheme. Additionally, rather
than using inhomogeneous increments to condensate (combining detrainment
and subsidence advection) to calculate cloud fraction increments, an
alternative is to use the detrainment of condensate to simply grow cloud
fraction to ensure a specified in-cloud liquid water content. The cloud
fraction are then advected downwards byt he subsidence advection. The
increments to cloud fraction from detrainment and subsidence are then
combined.

Prognostic dust approach
~~~~~~~~~~~~~~~~~~~~~~~~

A prognostic dust approach is implemented in the micro-physics scheme
under large-sale-precipitation where by the heterogeneous nucleation
temperature can be defined to vary three dimensionally globally as an
arc-tangent function of the mineral dust distribution in the model
(documented in ). By default, both liquid and ice are detrained
simultaneously at the same height, and the fraction of condensate that
is ice linearly ramps as a function of temperature. i.e. condensate is
assumed to be all-liquid when T is greater than one tuneable threshold;
all-ice when T is less than another tuneable threshold, and vary
linearly in-between (the threshold values are given by starticeTkelvin
and alliceTdegC in the UM cloud-scheme namelist. The new heterogeneous
nucleation temperatures calculated in the large-scale-precipitation are
passed to the convection scheme and are used as the above detrainment
temperature thresholds by maintaining a similar linear ramp. For e.g.,
condensate is assumed to be all-liquid for T :math:`\geq` :math:`tnuc_n`
and all-ice for T :math:`\leq` :math:`tnuc_n` - 10.0

.. _`sec:conv_input_profs`:

Condensation adjustment in the profiles input to the convection scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The convection scheme itself is highly sensitive to the input
environment temperature and moisture profiles *before* the convection
increments (or PC2 response) are calculated. In particular, the parcel
buoyancy (and hence the CAPE and mass-flux scaling) maybe radically
different depending on whether a “large-scale” condensation /
evaporation adjustment is performed before the convection call.

Where there is large-scale ascent, the profiles after Semi-Lagrangian
advection may have become supersaturated and unrealistically unstable,
until the expected condensation adjustment is performed. If the
convection scheme “sees” these unrealistic intermediate profiles, it is
likely to predict an excessive, unrealistic mass-flux.

To address this problem, there are two namelist switches that enable
additional condensation adjustments from PC2 before the convection call:

- **l_pc2_sl_advection**: performs homogeneous forcing response to
  Semi-Lagrangian advection immediately after the advection calculation,
  instead of at the end of the timestep (see section
  `4.8 <#sec:pres>`__).

- **l_cloud_call_b4_conv**: performs an additional call to PC2
  initiation (and PC2 checks) before the convection scheme (see section
  `4.9 <#sec:init2>`__). This should catch any instances where
  large-scale ascent or other processes have brought the profiles after
  advection to near or beyond saturation, in grid-points where there was
  no liquid cloud already present (and so no homogeneous forcing
  response).

.. _`sec:pres`:

Response to pressure changes
----------------------------

A pressure change following the parcel during the timestep will result
in an adiabatic temperature change which will force condensation, hence
we must include this temperature change forcing within PC2. The majority
of this pressure change comes from vertical advection (although not
all). Remember that the advection (section `4.5 <#sec:advec>`__), on its
own, does not cause condensation, it merely moves the existing cloud
field.

Using the semi-Lagrangian advection in the same way as is performed for
:math:`\overline{q_{cl}}` etc., the PC2 scheme will obtain the value of
the model prognostic *Exner*, (:math:`\prod`) on the departure points
(:math:`\prod_{dep}`). *Exner* is defined as

.. math::

   \label{eq:exner}
   \prod = \frac{T}{\theta} = \left( \frac{p}{p_{ref}} \right)^{\kappa}

where :math:`\theta` is the potential temperature, :math:`p_{ref}` is a
reference pressure set to 1000 hPa, and :math:`\kappa =
\frac{c_p - c_v}{c_p}` , where :math:`c_v` is the heat capacity of dry
air at constant volume. The *Exner* quantity is kept as a prognostic
variable in the model (this is unchanged from the control model), and
the value of :math:`\prod` on the departure points represents the
initial value in the timestep, since there is no update to :math:`\prod`
until the end of the timestep. After the second physics updates have
been performed (*atmos-physics2*), the model (including the control)
recalculates the value of *Exner* (:math:`\prod^{[n+1]}`). From
:math:`\prod_{dep}` and :math:`\prod^{[n+1]}` we can calculate, using
the definition (`[eq:exner] <#eq:exner>`__), the values of departure
pressure and temperature:

.. math:: \overline{p}_{dep} = p_{ref} {\prod_{dep}}^{\frac{1}{\kappa}}

.. math:: \overline{T}_{dep} = \theta \prod_{dep}

Hence we obtain the net forcing values

.. math::

   \Delta \overline{T} = \overline{T}^{[n+1]} - \overline{T}_{dep}
   \label{eq:deltatsl}

and

.. math::

   \Delta \overline{p} = \overline{p}^{[n+1]} - \overline{p}_{dep} .
   \label{eq:deltapsl}

where :math:`\overline{T}^{[n+1]}` and :math:`\overline{p}^{[n+1]}` are
the temperature and pressure at the arrival point, after the dynamics
call. (`[eq:deltatsl] <#eq:deltatsl>`__) and
(`[eq:deltapsl] <#eq:deltapsl>`__) are passed to the homogeneous forcing
routine in order to calculate the condensation and cloud fraction
changes associated with the pressure change.

We include this forcing towards the end of the timestep. There are two
reasons for this: firstly, values of :math:`\prod^{[n+1]}` are not
calculated by the control model until after the physics is complete;
secondly, it makes sense to locate this process in the timestep in a
similar location to where the large-scale cloud scheme is included in
the control (i.e. after the implicit part of the boundary layer has
finished).

However, there is a counter argument that says we should include this
process immediately after the dynamics, since we can then apply a
forcing on an initial state that has not already been modified by the
dynamics, boundary layer and convection schemes. This improves the
numerics of the problem, since the homogeneous forcing is designed to
take time level n values as inputs.

These issue are optionally addressed by turning on the UM namelist
switch **l_pc2_sl_advection**. Under this switch, the PC2 homogeneous
forcing response to pressure change is split:

#. Forcing by the *Lagrangian* component of pressure change, performed
   immediately after the Semil-Lagrangian advection scheme (before the
   call to atmos_physics2). This calculates the pressure change from the
   departure point value of *Exner* described above, to the
   start-of-timestep value of *Exner* at the arrival point.

#. Forcing by the *Eulerian* component of pressure change, performed at
   the end of the timestep (after the dynamics Helmholtz solver). This
   calculates the pressure change from the start-of-timestep *Exner* at
   the arrival point, to the end-of-timestep *Exner*.

Having to calculate the pressure forcing twice obviously adds some
computational cost, but has several advantages:

- As noted above, the PC2 homogeneous forcing calls can now take as
  input the temperature and water-vapour content *before* the pressure
  change has been applied, as intended. This should improve the
  numerical accuracy.

- Most of the condensation or evaporation from the dynamics comes from
  the *Lagrangian* component of the pressure change, which has now moved
  from the end of the timestep to before the dynamics Helmholtz solver.
  This means that any latent heating from condensation forced by ascent
  is now accounted for by the solver within the same timestep. This
  improves the numerical accuracy of the dynamics-physics coupling.

- If the condensation forced by resolved ascent is only added on at the
  end of the timestep, the profiles passed into atmos_physics2 can
  contain out-of-balance thermodynamic states (e.g. if the profile has
  been lifted by advection, it maybe supersaturated / unrealistically
  unstable before the resulting condensation is added on). This may
  adversely affect the convection scheme, which must act upon the
  profiles passed into atmos_physics2.

The splitting of the pressure forcing call under the
**l_pc2_sl_advection** switch was originally implemented to make the
profiles passed to convection more realistic.

.. _`sec:init2`:

Initiation
----------

As discussed in section `3.4 <#sec:init>`__, there are occasions when
:math:`\overline{q_{cl}}` and :math:`C_l` need to be initiated from 0 or
1. The application of the initiation is given in section
`3.4 <#sec:init>`__. The initiation forms a new, separate block of PC2
code to perform this calculation, and is located immediately following
the pressure change response (section `4.8 <#sec:pres>`__). Also, if the
UM namelist switch **l_cloud_call_b4_conv** is set to true, an
additional call to PC2 initiation is performed before the convection
scheme, to ensure that the condensation response to advection and other
forcings earlier in the timestep has been accounted for in the profiles
passed to the convection scheme, even if there was no cloud already
present for homogeneous forcing to act upon. (see section
`4.7.12 <#sec:conv_input_profs>`__).

There are currently 3 options for the conditions under-which initiation
may occur. For all of these options, if using the bimodal cloud scheme
to do initiation within PC2, then the tests on :math:`RH_T` relative to
:math:`RH_{crit}` are replaced by equivalent tests for whether the
saturation boundary lies within the bounds of the bimodal scheme’s
assumed PDF, as described in section `3.4.3 <#sec:bimodal_init>`__.

“Original” initiation logic
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This option is selected by setting the UM namelist switch
**i_pc2_init_logic = 1** (Original)

The initiation will be called if the liquid cloud fraction is either 0
or 1 and appropriate :math:`RH` criteria hold, along with other
restrictions. :math:`C_l` is initiated away from 0 if

-  :math:`RH_T > RH_{crit} + RH_{crit \, tol}` **and**

-  Cumulus convection has *not* been diagnosed from the boundary-layer
  in the current column **and**

-  The current level is not below the surface mixed-layer LCL **and**

-  :math:`C_l = 0` **and**

-  :math:`RH_T^{[n+1]} > RH_T^{[n]}` ,

where :math:`RH_{crit \, tol}` is a specified tolerance parameter, of
value 0.01, and :math:`RH_T` is defined in (`[eq:rht] <#eq:rht>`__).
:math:`RH_T^{[n]}` is the start of timestep value of :math:`RH_T` (i.e.
at time level n) and :math:`RH_T^{[n+1]}` is the value when initiation
is called. Additionally, there is another possibility for the last of
the relations. This second option also allows initiation when the water
is supercooled:

-  :math:`C_l < 0.05` *and* :math:`\overline{T} < 0 ^{\circ} C` .

Equivalently, :math:`C_l` is initiated away from 1 if

-  :math:`RH_T < 2 - RH_{crit} - RH_{crit \, tol}` **and**

-  :math:`C_l = 1` **and**

-  :math:`RH_T^{[n+1]} < RH_T^{[n]}` .

“Simplified” initiation logic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This option is selected by setting the UM namelist switch
**i_pc2_init_logic = 2** (Simplified)

Under this option, the conditions for initiation are:

Either:

- :math:`RH_T > RH_{crit} + RH_{crit \, tol}` **and**

- :math:`C_l < C_{tol}` **and**

- The current level is not below the surface mixed-layer LCL **and**

- :math:`RH_T^{[n+1]} > RH_T^{[n]}`

Or:

- :math:`RH_T < 2 - RH_{crit} - RH_{crit \, tol}` **and**

- :math:`C_l > 1 - C_{tol}` **and**

- :math:`RH_T^{[n+1]} < RH_T^{[n]}`

where :math:`C_{tol}` can be set via the UM namelist; its original
standard value is 0.005. Note this threshold is also used to remove
small cloud-fractions after initiation; see section
`4.9.4 <#sec:checks2>`__.

This is very similar to the “Original” initiation logic described above,
but with the following differences:

- The condition that the boundary-layer hasn’t diagnosed cumulus
  convection in the column is removed. Note that this condition
  spuriously suppressed initiation in the free troposphere *above* any
  cumulus cloud produced by the convection scheme.

- :math:`C_l` only needs to be within a numerical tolerance
  :math:`C_{tol}` from 0 or 1, rather than having to be *exactly* 0 or
  1.

- The different threshold when initiating super-cooled cloud is removed.

.. _`sec:smooth_initiation`:

“Smooth” initiation logic
~~~~~~~~~~~~~~~~~~~~~~~~~

This option is selected by setting the UM namelist switch
**i_pc2_init_logic = 3** (Smooth)

There is a fundamental numerical problem with the above options, in that
the initiation process is not permitted to have any effect at all unless
:math:`C_l` goes to (near) 0 or 1, but can predict values of :math:`C_l`
very different to 0 or 1 when it does activate. This leads to unphysical
sudden noisy jumps in :math:`C_l` and :math:`q_{cl}` when initiation
occurs. For example, if erosion causes :math:`C_l` to steadily decline,
it will continue to decline (even when the grid-mean :math:`RH_T`
exceeds :math:`RH_{crit}`) until it reaches the threshold (0 or
:math:`C_{tol}`). At this point, initiation suddenly increases
:math:`C_l` and :math:`q_{cl}` to the values predicted by the diagnostic
cloud scheme. Erosion may then gradually remove them again, and the
cycle repeats. There is no physical reason for this internal mode of
variability in the scheme.

Another problem arises if we consider the sensitivity to model
resolution. Suppose we have many adjacent small grid-boxes with similar
:math:`RH_T`, a few containing cloud, the rest containing no cloud. If
the whole region cools to the point where :math:`RH_T > RH_{crit}`, then
new cloud will initiate in the cloud-free grid-boxes, but not in the
cloudy grid-boxes. Now suppose we run a coarse-grained version of the
same simulation; the many small grid-boxes are replaced by a single
grid-box containing the average :math:`C_l` over the small grid-boxes.
Since we now have just one grid-box already containing partial
cloud-cover, initiation of new cloud can no longer occur anywhere.

To address these problems, there is an option to use a much simpler /
numerically better-posed initiation method; always allow the diagnostic
cloud scheme to be called (provided it is expected to predict nonzero
cloud water, i.e. :math:`RH_T > RH_{crit}` in the case of the Smith
scheme). The :math:`q_{cl}` predicted by the diagnostic cloud scheme is
then taken as a minimum limit applied to the prognostic :math:`q_{cl}`.
This amounts to taking the diagnostic cloud scheme’s assumed PDF as a
minimum allowed width to the actual prognostic moisture PDF. The
prognostic :math:`C_l` and :math:`q_{cl}` are incremented as follows:

- If :math:`{q_{cl}}_{diag} > q_{cl}`:

  :math:`\Delta q_{cl} = {q_{cl}}_{diag} - q_{cl}
  \quad \refstepcounter{equation}(\theequation)\label{eq:dqcl_init}`

  - If :math:`Q_C < 0`:

    :math:`\Delta C_{l} = \frac{\Delta q_{cl}}{{q_{cl}}_{diag}}
                    \left( {C_{l}}_{diag} - C_{l} \right)
      \quad \refstepcounter{equation}(\theequation)\label{eq:dcl_init1}`

  - If :math:`Q_C > 0`:

    :math:`\Delta C_{l} = \frac{\Delta SD}{{SD}_{diag}}
                    \left( {C_{l}}_{diag} - C_{l} \right)
      \quad \refstepcounter{equation}(\theequation)\label{eq:dcl_init2}`

- Otherwise:

  :math:`\Delta q_{cl} = 0`

  :math:`\Delta C_{l} = 0`

where the subscript :math:`_{diag}` denotes the liquid cloud water
content and fraction predicted by the diagnostic cloud scheme (either
Smith or Bimodal).

Equation `[eq:dcl_init1] <#eq:dcl_init1>`__ simply sets the
cloud-fraction to a weighted mean of the pre-existing and
diagnostic-scheme cloud-fractions, in proportion to the fraction of the
water content that was created by initiation versus that which was
already there. If the pre-existing :math:`q_{cl}` is zero,
`[eq:dqcl_init] <#eq:dqcl_init>`__ and
`[eq:dcl_init1] <#eq:dcl_init1>`__ simply set :math:`q_{cl}` and
:math:`C_l` to their new diagnosed values, as in the previous options.
Crucially, in the limit that the pre-existing :math:`q_{cl}` approaches
:math:`{q_{cl}}_{diag}`, the increments to :math:`q_{cl}` and
:math:`C_l` smoothly go to zero. This is important to make the
initiation process numerically well-posed, so that it yields a smooth,
continuous solution.

Note that when we are initiating from :math:`C_l = 1` instead of
:math:`C_l = 0`, we expect the pre-existing :math:`q_{cl}` to be nonzero
even when there is no pre-existing sub-grid PDF width. In this case, the
completely uninitiated state will have zero saturation deficit
:math:`SD`, rather than zero :math:`q_{cl}`. Therefore, in this case the
increment to :math:`C_l` is calculated based on the fractional increase
in :math:`SD` from initiation (equation
`[eq:dcl_init2] <#eq:dcl_init2>`__), instead of the fractional increase
in :math:`q_{cl}`.

Whether to increment :math:`C_l` based on the increase in :math:`q_{cl}`
or :math:`SD` is determined based on the sign of :math:`Q_C`, which is
defined as in equation `[eq:qc_eq_qt-qs] <#eq:qc_eq_qt-qs>`__
(reproduced here for clarity):

.. math:: Q_c = a_L \left( \overline{q_T} - q_{sat}(\overline{T_L}) \right)

The saturation deficit :math:`SD` is defined by equation
`[SD2] <#SD2>`__:

.. math:: SD = a_L \left( q_{sat}(\overline{T}) - \overline{q} \right)

Under the reasonable approximation that :math:`q_{sat}` varies linearly
between :math:`\overline{T}` and :math:`\overline{T_L}`, so that the
values of :math:`\alpha` and :math:`a_L` are the same in both of these
equations, and:

.. math:: q_{sat}(\overline{T_L}) = q_{sat}(\overline{T}) - \alpha \frac{L}{c_p} q_{cl}

we obtain:

.. math::

   q_{cl} = Q_c + SD
   \label{eq:qc_plus_sd}

It can be seen that when :math:`Q_C > 0` (total-water super-saturation),
it represents the value :math:`q_{cl}` would have if the whole grid-box
were saturated (:math:`SD = 0`, :math:`C_l = 1`). Note that
:math:`q_{cl}` cannot fall below :math:`Q_C`, since :math:`SD` cannot be
negative. Since :math:`Q_c` is invariant under condensation /
evaporation, we must have :math:`\Delta SD = \Delta q_{cl}` (hence the
implementation of `[eq:dcl_init2] <#eq:dcl_init2>`__ in the code simply
uses :math:`q_{cl} - Q_c` in place of :math:`SD`, and
:math:`\Delta q_{cl}` in place of :math:`\Delta SD`).

.. _`sec:checks2`:

Additional checks after PC2 initiation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The initiation is followed immediately by a section of resetting code.
For numerical reasons, it is possible to obtain very low, but non zero,
values of :math:`C_l` (and equivalently values very close to, but not
equal to, 1). The code will reset these clouds to either a fraction of 0
or 1, as appropriate. We choose to apply these terms here and not in the
Bounds Checking part of the code (section `4.10 <#sec:checks>`__)
because these are not required to obtain consistency between fields, but
are ‘tidying up’ pieces of code, although they may reasonably also be
applied in the Bounds Checking. Care needs to be taken when choosing the
thresholds, since we do not wish to reset small values that are
genuinely created by a physics scheme in the model.

We first calculate :math:`RH_T` using (`[eq:rht] <#eq:rht>`__) and
compare this to the critical relative humidity, :math:`RH_{crit}`. The
liquid cloud fraction will be reset to 1 if:

-  :math:`RH_T > 2 - RH_{crit}` and :math:`C_l \ge C_{high}`

-  or :math:`C_l \ge C_{high 2}`

where :math:`C_{high}` and :math:`C_{high 2}` are defined in
`[eq:chigh-chigh2] <#eq:chigh-chigh2>`__. The evaporation is done by
calculating :math:`SD` using (`[SD2] <#SD2>`__) with
(`[eq:a_L] <#eq:a_L>`__) and (`[eq:alpha_exp] <#eq:alpha_exp>`__) and
evaporating the equivalent amount of liquid into the gridbox to take it
to saturation, according to (`[eq:qsdcheck1] <#eq:qsdcheck1>`__) below.

Similarly, the equivalent check for low values of :math:`RH_T` is
performed. The liquid cloud fraction will be reset to 0 if:

-  :math:`RH_T < RH_{crit}` and :math:`C_l \le C_{low}`

-  or :math:`C_l \le C_{low 2}` .

The remaining :math:`\overline{q_{cl}}` is evaporated into the gridbox
using (`[eq:qclcheck] <#eq:qclcheck>`__) below.

The thresholds :math:`C_{high}`, :math:`C_{high 2}`, :math:`C_{low}` and
:math:`C_{low 2}` are set using the parameters :math:`C_{tol}` and
:math:`C_{tol 2}`, according to:

.. math::

   \begin{aligned}
   C_{high} = 1 - C_{tol},  \nonumber \\
   C_{high 2} = 1 - C_{tol 2},  \nonumber \\
   C_{low} = C_{tol},  \nonumber \\
   C_{low 2} =  C_{tol 2},
   \label{eq:chigh-chigh2}
   \end{aligned}

where the parameters :math:`C_{tol}` and :math:`C_{tol 2}` can be set
via the UM namelist variables **cloud_pc2_tol** and **cloud_pc2_tol_2**.
The original standard values of these parameters are
:math:`C_{tol} = 0.005` and a lower value :math:`C_{tol 2} = 0.001`.

Investigations in SCM runs using the comorph convection scheme (which
behaves more smoothly and so typically gives smaller increments to
:math:`C_l` over a single timestep than other schemes which exhibit
intermittent behaviour) suggested these thresholds are too high to avoid
spuriously resetting physical values of :math:`C_l` to zero. Detrainment
from sparse shallow cumulus, or advection of cloud into a neighbouring
grid-box under light winds, commonly give increments which increase
:math:`C_l` from zero to a value less than :math:`0.005` in one timestep
(but would eventually increase :math:`C_l` to a significant value over
subsequent timesteps if the checks did not keep resetting :math:`C_l` to
zero).

Note that if these checks are relaxed by lowering the thresholds
:math:`C_{tol}` and :math:`C_{tol 2}` to near-zero, similar checks are
still performed independently by the bounds checking described in
section `4.10 <#sec:checks>`__, but with a much lower threshold of
:math:`C_{tol 3} = 1 \times 10^{-12}`.

.. _`sec:checks`:

Bounds checking
---------------

Ideally, model prognostics would never become inconsistent with one
another. However, even although the mathematical solution of the
governing equations may be well behaved, due to numerical inaccuracies
values may become inconsistent. For the cloud and condensate quantities,
there are a number of consistencies that must apply. The bounds checking
forms a subroutine that will, if necessary, adjust :math:`\overline{q}`,
:math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}`, :math:`C_l`,
:math:`C_i`, :math:`C_t` and, for latent heating, :math:`\overline{T}`,
to ensure consistency between these values.

The bounds checking is performed three times during the timestep.
Firstly, after the parallel part of the physics (*atmos-physics1*) is
complete; secondly, before the initiation (section `3.4 <#sec:init>`__)
is called; thirdly, after the initiation is called.

Firstly, if :math:`C_l > 1 - C_{tol 3}` then :math:`C_l` is set to 1.
Accordingly, :math:`C_t` is set to 1 as well. :math:`C_{tol 3}` is a
tiny numerical tolerance set to :math:`1 \times 10^{-12}`, a value
intended to be in the realm of floating point rounding error rather than
anything that represents a physical solution.

.. _section-1:

The second check is to reset :math:`\overline{C_l}` to zero. This may be
performed for two reasons. Firstly, if the amount of
:math:`\overline{q_{cl}}` is very small
(:math:`\overline{q_{cl}} < q_{c0}`, where
:math:`q_{c0} = 1 \times 10^{-10} kg kg^{-1}`), so we avoid carrying
negligible, but non-zero values of :math:`\overline{q_{cl}}` and
:math:`C_l`. Secondly, if :math:`C_l < C_{tol 3}` then we reasonably
reset :math:`C_l` to zero. :math:`C_t` gets reset, as it must if there
is no liquid cloud, to be equal to :math:`C_i`.

.. _`sec:pc2_checks_sd`:

The next check complements the first but updates the moisture fields. We
firstly calculate :math:`SD` using (`[SD2] <#SD2>`__) and
(`[eq:alpha_exp] <#eq:alpha_exp>`__). We then check whether
:math:`SD < 0`. This check catches instances where we have grid-mean
supersaturation, which ought to be impossible (under the instantaneous
condensation assumption made by PC2, condensation should occur to
instantly adjust any supersaturated regions of the gridbox to
saturation, so we *must always* have :math:`SD \ge 0`. When this
condition is violated, we condense water vapour to adjust to grid-mean
saturation. :math:`-SD` corresponds to the amount of vapour that must be
condensed to achieve this, so we have:

.. math::

   \begin{aligned}
   \overline{q} \leftarrow \overline{q} + SD  \nonumber \\
   \overline{q_{cl}} \leftarrow \overline{q_{cl}} - SD \nonumber \\
   \overline{T} \leftarrow \overline{T} - \frac{L_c}{c_p} SD
   \label{eq:qsdcheck1}
   \end{aligned}

The original version of this check on :math:`SD` (which may increase
:math:`\overline{q_{cl}}`), made no accompanying changes to liquid cloud
fraction. However, increases in :math:`\overline{q_{cl}}` without any
increase in :math:`C_l` can lead to spurious high in-cloud condensate
which is then converted to rain by the microphysics at the next
time-step. There are currently 4 options for how to treat :math:`C_l`
when increasing :math:`\overline{q_{cl}}` under this saturation
adjustment, selected by the UM large-scale cloud namelist switch
**i_pc2_checks_cld_frac_method**:

- **i_pc2_checks_cld_frac_method = 0** - Original method; :math:`C_l` is
  left unaltered.

- **i_pc2_checks_cld_frac_method = 1** - Set :math:`C_l` and :math:`C_t`
  to 1.

- **i_pc2_checks_cld_frac_method = 2** - If :math:`\overline{q_{cl}}`
  and :math:`C_l` were already nonzero before the adjustment, increase
  :math:`C_l` at the same fractional rate as :math:`\overline{q_{cl}}`,
  so that the in-cloud water content
  :math:`\frac{\overline{q_{cl}}}{C_l}` is conserved. Otherwise,
  increase :math:`C_l` so-as to yield a prescribed in-cloud water
  content set to 0.5 g kg\ :math:`^{-1}`. :math:`C_t` is then increased
  by the same amount as :math:`C_l`, to maintain consistency.

- **i_pc2_checks_cld_frac_method = 3** - This is the same as option 2
  above, except in the case where :math:`\overline{q_{cl}}` or
  :math:`C_l` was zero before the adjustment. In this case, :math:`C_l`
  is set based on an empirical power-law function of
  :math:`\overline{q_{cl}}`.

.. _section-2:

Next we check whether :math:`SD > 0`, *and* :math:`C_l = 1` (the first
of our checks has ensured that :math:`C_l` is no greater than 1). This
check catches instances where we have total cloud-cover in a
subsaturated grid-box, which ought to be impossible (if the whole
grid-box is full of liquid cloud, then it must be at grid-mean
saturation, i.e. :math:`SD = 0`). When this happens, we adjust
:math:`\overline{q}` and :math:`\overline{q_{cl}}` to take :math:`SD` to
zero, *provided* that :math:`\overline{q_{cl}} > SD`. Remember that
:math:`SD` corresponds to the amount of vapour that must be *evaporated*
into the gridbox to give saturation, so we simply make exactly the same
adjustments as we do for removing supersaturated states above
(`[eq:qsdcheck1] <#eq:qsdcheck1>`__), except that here :math:`SD` is
positive rather than negative.

Our proviso that :math:`\overline{q_{cl}} > SD` ensures that we do not
make :math:`\overline{q_{cl}}` negative by this adjustment. If
:math:`\overline{q_{cl}} < SD`, then we cannot bring the gridbox to
saturation, but it is still wrong to allow :math:`C_l = 1` in a
subsaturated gridbox! This was identified as a bug in the
bounds-checking code, which sometimes caused instances of
:math:`C_l = 1` to spuriously persist in dry environments. This
behaviour is currently controlled by a temporary logical in the
**temp_fixes** namelist:

- If **l_pc2_checks_sdfix** is set to false, the code simply does
  nothing when it finds instances of :math:`C_l = 1`, :math:`SD > 0` and
  :math:`SD > \overline{q_{cl}}`, allowing such artefacts to persist.

- If **l_pc2_checks_sdfix** is set to true, in these instances we simply
  evaporate all the remaining liquid water, and reset :math:`C_l` to
  zero:

  .. math::

     \begin{aligned}
     \overline{q} \leftarrow \overline{q} + \overline{q_{cl}}  \nonumber \\
     \overline{T} \leftarrow \overline{T} - \frac{L_c}{c_p} \overline{q_{cl}}
     \nonumber \\
     \overline{q_{cl}} \leftarrow 0 \nonumber \\
     C_l \leftarrow 0\nonumber \\
     C_t \leftarrow C_i
     \label{eq:qsdcheck2}
     \end{aligned}

.. _section-3:

The next check is similar to above but for the :math:`C_l = 0`
situation.

If :math:`\overline{q_{cl}} < q_{c0}` or :math:`C_l = 0` then we
evaporate the small amount of :math:`\overline{q_{cl}}` that remains in
the gridbox:

.. math::

   \begin{aligned}
   \overline{q} \leftarrow \overline{q} + \overline{q_{cl}} \nonumber \\
   \overline{q_{cl}} \leftarrow 0 \nonumber \\
   \overline{T} \leftarrow \overline{T} - \frac{L_c}{c_p} \overline{q_{cl}}
   \label{eq:qclcheck}
   \end{aligned}

.. _section-4:

Next, if :math:`C_i > 1` then :math:`C_i` is set to 1. Accordingly,
:math:`C_t` is set to 1 as well.

.. _section-5:

The following check is on the ice water content,
:math:`\overline{q_{cf}}`, and the ice fraction :math:`C_i`. If
:math:`\overline{q_{cf}} < q_{c0}` we simply condense some vapour to
remove the negative quantity.

.. _section-6:

However, instead of removing small amounts of :math:`\overline{q_{cf}}`
when :math:`C_i = 0` but :math:`\overline{q_{cf}} > 0`, we choose
instead to create some :math:`C_i` to keep consistency. This is to allow
small, but significant, amounts of :math:`\overline{q_{cf}}` created by
the microphysics scheme to be maintained.

.. math::

   C_i \leftarrow \frac { \overline{q_{cf}} }{q_{cf0}}
   \label{eq:cf_reset}

where the ‘in-cloud’ ice content
:math:`q_{cf0} = 1 \times 10^{-4} kg kg^{-1}`.

.. _section-7:

The next two checks are on the total cloud fraction, :math:`C_t`, to
ensure that it takes on a value that is physically possible, given the
values of :math:`C_l` and :math:`C_i`. We have, firstly, the maximum
overlap situation and then the minimum overlap situation.

.. math::

   \begin{aligned}
   C_t \leftarrow \text{Max}( C_t, C_i, C_l ) \nonumber \\
   C_t \leftarrow \text{Min}( C_t , C_l + C_i, 1)
   \label{eq:ctchecks}
   \end{aligned}

.. _section-8:

Finally, there is a homogeneous nucleation term applied, similar to that
in the large-scale precipitation (section `4.2.2 <#sec:lsp_homo>`__).
This is a fast microphysics process, and must act to ensure that no
liquid cloud created by the initiation is allowed to persist in this
phase if the temperature is cold enough. Hence, if
:math:`\overline{T} < T_{homo}` then

.. math::

   \begin{aligned}
   \overline{q_{cf}} \leftarrow \overline{q_{cf}} + \overline{q_{cl}} \nonumber \\
   \overline{q_{cl}} \leftarrow 0 \nonumber \\
   \overline{T} \leftarrow \overline{T} + \frac{L_f}{c_p} \overline{q_{cl}} \nonumber \\
   C_i \leftarrow C_t \nonumber \\
   C_l \leftarrow 0.
   \label{eq:homochecks}
   \end{aligned}

.. _`sec:qpos`:

Qpos checks
~~~~~~~~~~~

The implementation of the PC2 code includes an additional bounds check
after the *atmos-physics-2* part of the model timestep has been
completed. This check is necessary to trap a rare failure, and uses the
*Qpos* subroutines to check that :math:`\overline{q_{cl}}` is greater or
equal to 0.

During trialling prior to operational implementation, it was found that
relying on Q-Pos to deal with negative condensate values was very
expensive, as the Q-Pos routine does a lot of communications between
different processors. It may be preferable to deal with the cause of
negative condensate amounts at their source. The option to “Ensure
consistent sinks of qcl and CFL” prevents the QCL increment from trying
to remove too much liquid condensate and hence reduces the models
reliance on Q-Pos to deal with the inconsistencies.

.. _`sec:da`:

Data Assimilation
-----------------

The data assimilation section in the model will output assimilation
increments that represent changes to :math:`\overline{q}` and
:math:`\overline{T}` which *include* the condensation contributions. We
hence need to calculate equivalent increments to
:math:`\overline{q_{cl}}`, :math:`C_l` and :math:`C_t`. We assume that
the assimilation has not calculated these using a different method. We
consider the homogeneous framework and assume that there is a forcing
value of :math:`Q_c` that exists that will produce the known increment
to :math:`\overline{q}` and :math:`\overline{T}`.

Discritising (`[dqcldt] <#dqcldt>`__) we have, using
(`[eq:deltaqc_exp] <#eq:deltaqc_exp>`__) and expanding
:math:`\Delta T_L` in terms of :math:`\Delta T` and
:math:`\Delta q_{cl}`,

.. math::

   \Delta \overline{q_{cl}} = C_l ( a_L ( \Delta \overline{q} -
   \alpha \Delta \overline{T} - \beta \Delta \overline{p}) + \Delta \overline{q_{cl}} ).
   \label{eq:da1}

Remember that :math:`Q_c` (and hence :math:`\Delta Q_c`) is independent
of condensation. Rearranging, we obtain

.. math::

   \Delta \overline{q_{cl}} = \frac{1}{1 - C_l} C_l 
   a_L ( \Delta \overline{q} - \alpha \Delta \overline{T} - \beta \Delta \overline{p})
   \label{eq:da2}

and hence an expression for the condensate increment,
:math:`\Delta \overline{q_{cl}}`, that accompanies the known increments
to :math:`\overline{q}` and :math:`\overline{T}`. The similar analysis,
from (`[dcdt] <#dcdt>`__) and (`[eq:da1] <#eq:da1>`__) gives

.. math::

   \Delta C_l = \frac{1}{1 - C_l} G(-Q_c) 
    a_L ( \Delta \overline{q} - \alpha \Delta \overline{T}
   - \beta \Delta \overline{p} ) .
   \label{eq:da3}

Hence the equation set is equivalent to the use of the homogeneous
forcing set, except for the multiplier :math:`\frac{1}{1 - C_l}`.
Although this is a clean solution, we need to be very careful with the
ill-conditioning of this solution near :math:`C_l = 1`.

In practice, the ill-conditioning of (`[eq:da2] <#eq:da2>`__) and
(`[eq:da3] <#eq:da3>`__) becomes too numerically awkward for us to apply
the full solution based on homogeneous forcing, although, for
completeness, we outline it in Appendix `6 <#sec:appendix-da>`__. Hence
we have chosen to apply a much simpler model. Here we use simply the
data assimilation increments :math:`\Delta \overline{q}` and
:math:`\Delta \overline{T}` within the standard homogeneous forcing
(section `3.2 <#sec:homog>`__), even though we are fully aware that this
is inconsistent (because :math:`\Delta \overline{q}` and :math:`\Delta
\overline{T}` are not forcings, but are forcings plus the condensation.
This allows us an *estimate* of :math:`\Delta \overline{q_{cl}}` and
:math:`\Delta{C_l}`, via the homogeneous forcing routine (and
:math:`\Delta C_t` via the standard updating described in section
`3.6 <#sec:ct>`__). These are the quantities applied as the equivalent
data assimilation increments for :math:`\Delta \overline{q_{cl}}`,
:math:`\Delta{C_l}` and :math:`\Delta C_t`. The increments
:math:`\Delta \overline{q}` and :math:`\Delta \overline{T}` remain those
that the data assimilation scheme itself calculated.

Appendix `6 <#sec:appendix-da>`__ gives, for completeness, the
alternative numerical technique for the solution of
(`[eq:da2] <#eq:da2>`__) and (`[eq:da3] <#eq:da3>`__). However, we
stress that this technique is not used within the current PC2
formulation.

.. _`sec:um`:

Implementation in the Unified Model
===================================

This section considers the implementation of PC2 within the Unified
Model code and provides a brief guide to its use.

In general, we have written PC2 so that the timestepping of the cloud
fraction variables within the *atm_step_4a* subroutine is treated as
much as possible in a similar way to the condensate variables. Hence,
wherever the condensed water variables :math:`q_{cl}` and :math:`q_{cf}`
are updated, the cloud fractions need to be updated consistently.

.. _`sec:acf`:

Area cloud fraction
-------------------

Two area cloud fraction parametrizations are available for use with PC2.

The area cloud fraction of Cusack (documented in ) has been adapted by
:raw-latex:`\cite{boutle_morcrette10}` so it can be used with PC2 (and
is available from the UMUI as the “Cusack” option from version 7.6
onwards). This method aims to reproduce some of the detail of the
thermodynamic profile lost due to the coarseness of the grid. The
interpolation/extrapolation technique is used prior to PC2 initiation
(which is then called with three times as many levels) and it is used,
along with the homogeneous forcing idea at the start of the timestep to
allow more cloud to be seen by radiation.

The diagnostic area cloud fraction of :raw-latex:`\cite{bhi05}` has also
been implemented in the model (available from the UMUI at version 6.4
onwards), and this is used in PC2:64. This method diagnoses the area
cloud fraction given the volume cloud fraction, taking into account the
size of the grid box. The setting of the area cloud fraction is
performed at the end of the timestep.

.. _`sec:code`:

Code Structure
--------------

A detailed description of the UM’s timestep structure, showing where in
the model all the PC2 cloud scheme subroutine calls are made, is given
in the subsections below.

Note that there are three different subroutines that all do the PC2
homogeneous forcing, with slightly different details:

- ***pc2_delta_hom_turb*** outputs increments due to the condensation or
  evaporation, but doesn’t update the fields themselves.

- ***pc2_homog_plus_turb*** just updates the fields that are passed in,
  instead of outputting separate increment arrays.

- ***pc2_hom_conv*** outputs increments but includes additional
  calculations for various cloud erosion formulations.

Note that code exists in the first two of these routines to do erosion,
but they can only do it via an input fixed rate of narrowing of the
moisture PDF (which is currently set to zero in all instances). PC2
development has settled on a more complicated treatment of erosion,
which has only been implemented in *pc2_hom_conv*. This can either be
called after the convection scheme (within *pc2_from_conv_ctl*), or
before the microphysics scheme (within *pc2_turbulence_ctl*).

Note there is also an optional call to *pc2_turbulence_ctl* after the
microphysics scheme, which is used only to estimate the cloud fraction
change consistent with the turbulent production of liquid cloud (see
section `3.8 <#sec:turb_qcl_scheme>`__).

Most PC2 code is protected by IF tests on the namelist input *i_cld_vn*
= 2 (PC2 in the GUI). However, within the convection scheme, the code is
controlled by logicals *l_calc_dxek* (which is just set to true if using
PC2, and set false otherwise), and *l_q_interact*, which controls
whether to allow the interactive detrainment and entrainment of
condensate.

There is also a switch (currently hardwired to .false. in the code)
called *l_pc2_reset*. Turning this on (not recommended!) does 2 things:

- Convective entrainment and detrainment of condensate is disabled, by
  setting *l_q_interact* to false.

- The prognostic cloud variables are overwritten by a call to the
  diagnostic cloud scheme at the end of the timestep, in subroutine
  *qt_bal_cld*. NOTE: this functionality will no longer work, because
  inside *qt_bal_cld* the call to the diagnostic cloud scheme is now
  protected by IF tests on using either the Smith or bimodal cloud
  schemes. If using PC2, no cloud scheme is called here, and required
  output variables are just left unset!

The location of the various cloud scheme routine calls within the UM is
summarised in the list below.

Subroutines only called for the Smith scheme are highlighted in blue,
those only called for PC2 are in green, and those only called for the
bimodal scheme are in purple.

Main Tree from atm_step_4a
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container:: itemize

   | **atm_step_4a**
   | \* (performs one timestep of the Unified Model...)

   .. container:: itemize

      .. container:: tcolorbox

         | **atm_step_alloc_4a**
         | \* (does miscellaneous initialisations in atm_step)

         - | pc2_rhtl
           | \* (calculate start-of-timestep Relative Humidity, used by
             PC2 initiation)

      .. container:: tcolorbox

         | **atmos_physics1**
         | \* (calls explicit “slow” physics routines...)

         .. container:: itemize

            .. container:: tcolorbox

               | **microphys_ctl**
               | \* (interface to microphysics scheme)

               - | pc2_turbulence_ctl
                 | \* (Perform optional erosion of liquid-cloud; done
                   here if NOT doing erosion after convection, e.g. if
                   no convection scheme is used).

                 - | pc2_hom_conv
                   | \* (called here just to do erosion)

               - | ls_cld
                 | \* (Smith scheme without area cloud fraction
                   calculation, to set initial cloud fields passed into
                   microphysics)

               - | **ls_ppn**
                 | \* (microphysics scheme)

               - | **mphys_turb_gen_mixed_phase**
                 | \* (turbulent production of liquid cloud)

               - | pc2_turbulence_ctl
                 | \* (optionally use the PC2 pdf-width-change code to
                   calculate the cloud-fraction change from the above
                   turbulent production of liquid cloud)

                 - | pc2_hom_conv
                   | \* (called here just to calculate the cloud
                     fraction increment consistent with the turbulent
                     qcl increment)

            .. container:: tcolorbox

               | **rad_ctl**
               | \* (interface to radiation scheme)

               - | **sw_rad**
                 | \* (short-wave radiation scheme)

               - | pc2_homog_plus_turb
                 | \* (PC2 homogeneous forcing of liquid-cloud by SW
                   radiation heating)

               - | **lw_rad**
                 | \* (long-wave radiation scheme)

               - | pc2_homog_plus_turb
                 | \* (PC2 homogeneous forcing of liquid-cloud by LW
                   radiation tendency)

            .. container:: tcolorbox

               **atmos_physics1_alloc_pc2** (wrapper for PC2
               self-consistency checks at end of atmos_physics1)

               - Add increments from microphysics + radiation onto
                 start-of-timestep fields to form updated fields.

               - | pc2_checks
                 | \* (self-consistency checks on cloud fractions and
                   water contents)

               - Convert corrected updated fields back to increments.

      Begin loop over solver outer cycles

      .. container:: itemize

         .. container:: tcolorbox

            | **atm_step_phys_reset**
            | \* (for PC2, on subsequent solver outer cycles, reset
              cloud-fractions to saved values after atmos_physics1)

         .. container:: tcolorbox

            | **eg_sl_moisture**
            | \* (large-scale advection of cloud water contents and
              fractions)

         .. container:: tcolorbox

            | pc2_pressure_forcing_only
            | \* (Optionally calculate homogeneous forcing of liquid
              cloud by the pressure change along the trajectory from
              departure point to arrival point).

            - | pc2_homog_plus_turb
              | \* (generic homogeneous forcing routine used here).

         .. container:: tcolorbox

            | **atmos_physics2**
            | \* (calls “fast” physics routines...)

            .. container:: itemize

               .. container:: tcolorbox

                  | **ni_bl_ctl**
                  | \* (interface to explicit boundary-layer and surface
                    scheme calls, including calculation of TKE and
                    TKE-based :math:`RH_{crit}`)

               .. container:: tcolorbox

                  | bm_calc_tau
                  | \* (calculates turbulence properties used in the
                    bimodal cloud scheme, based on the boundary-layer
                    scheme TKE and mixing-length)

               .. container:: tcolorbox

                  | **cloud_call_b4_conv**
                  | \* (routine for optional cloud-scheme calls before
                    convection)

                  - | ls_arcld
                    | \* (Smith scheme with area cloud fraction; see
                      `5.2.2 <#subsubsec:smith_acf>`__ for a drill-down
                      inside this routine)

                  - | bm_ctl
                    | \* (bimodal scheme)

                  - Set area cloud fraction equal to bulk cloud fraction

                  - | pc2_initiation_ctl
                    | \* (interface to PC2 initiation and
                      consistency-checks; see
                      `5.2.3 <#subsubsec:pc2_initiation>`__ for a
                      drill-down inside this routine)

               .. container:: tcolorbox

                  | **ni_conv_ctl** or **other_conv_ctl**
                  | \* (interface routines to various convection
                    schemes...)

                  - | **glue_conv_5a/6a**
                    | \* (calls deep, shallow and mid-level convection
                      schemes)

                    - | **deep/shallow/mid_conv**
                      | \* (convection scheme main routines)

                      - **convec2** (completes lifting of the convective
                        parcel by one model-level)

                        - **parcel** (calculates new parcel properties
                          at next level)

                        - **environ** (calculates grid-mean increments
                          to primary fields; includes PC2 partitioning
                          of detrained condensate mass between liquid
                          and ice phases)

                        - pc2_environ (calculates increments to PC2
                          cloud fractions due to convective detrainment
                          and subsidence)

                  - | pc2_from_conv_ctl
                    | \* (PC2 calculations after convection)

                    - | pc2_hom_conv
                      | \* (homogeneous forcing by convection, and
                        erosion of liquid-cloud)

               .. container:: tcolorbox

                  | **ni_imp_ctl**
                  | \* (interface to boundary-layer implicit solver)

                  - | **imp_solver**
                    | \* (implicitly solves vertical diffusion to find
                      :math:`T_l` and :math:`q_T` updated by turbulent
                      fluxes).

                  - | pc2_bl_inhom_ice
                    | \* (inhomogeneous forcing of ice-cloud)

                  - | pc2_delta_hom_turb
                    | \* (homogeneous forcing of liquid cloud by the
                      turbulent fluxes)

                  - | pc2_bl_forced_cu
                    | \* (adds diagnosed “forced cumulus” cloud fraction
                      and water content onto the PC2 prognostics)

                  - Calculate area cloud fraction:

                    | ls_acf_brooks
                    | \* (for the Brooks epirical method)

                    | pc2_hom_arcld
                    | \* (for the Cusack vertical interpolation method)

                    - | pc2_homog_plus_turb
                      | \* (generic homogeneous forcing routine used to
                        interpolate)

                  - | ls_arcld
                    | \* (interface to diagnostic Smith scheme and area
                      cloud fraction; see
                      `5.2.2 <#subsubsec:smith_acf>`__ for a drill-down
                      inside this routine)

                  - | bm_ctl
                    | \* (bimodal cloud scheme)

                  - Set area cloud fraction equal to bulk cloud fraction

                  - | **diagnostics_bl**
                    | \* (outputs boundary-layer diagnostics to STASH)

                    - | **ls_cld**
                      | \* (Smith scheme used here to calculate various
                        diagnostics of near-surface temperature and
                        humidity, by extrapolating pressure, :math:`T_l`
                        and :math:`q_t` down to the desired height and
                        then re-diagnosing :math:`q_{cl}`.

         .. container:: tcolorbox

            | **atm_step_ac_assim**
            | \* (interface to Data Assimilation analysis increments...)

            - **ac_ctl** (control routine for Data Assimilation analysis
              increments...)

              - **ac** (main analysis increment routine)

              - | pc2_assim
                | \* (PC2 reponse to the analysis increments; see
                  `5.2.4 <#subsubsec:pc2_assim>`__ for a drill-down
                  inside this routine)

              - ls_acf_brooks (calculate area cloud fraction using
                Brooks empirical method if active)

              - ls_arcld (call diagnostic Smith scheme with area cloud
                fraction again to account for the analysis increments;
                see `5.2.2 <#subsubsec:smith_acf>`__ for a drill-down
                inside this routine)

         .. container:: tcolorbox

            | **eg_sl_helmholtz**
            | \* (dynamics pressure solver; updates pressure, and the
              winds used to perform advection on the next solver outer
              cycle)

      End loop over solver outer cycles

      .. container:: tcolorbox

         | pc2_pressure_forcing
         | \* (interface to miscellaneous PC2 calculations at
           end-of-timestep)

         - | pc2_homog_plus_turb
           | \* (homogeneous forcing of liquid-cloud by the dynamics
             pressure change; optionally either uses total pressure
             change including the Lagrangian component following the
             winds, or only the Eulerian component from the dynamics
             solver)

         - | pc2_initiation_ctl
           | \* (interface to PC2 initiation and consistency-checks; see
             `5.2.3 <#subsubsec:pc2_initiation>`__ for a drill-down
             inside this routine)

      .. container:: tcolorbox

         | **qt_bal_cld**
         | \* (calculates end-of-timestep cloud state consistent with
           final pressure...)

         - | ls_arcld
           | \* (interface to diagnostic Smith scheme and area cloud
             fraction; see `5.2.2 <#subsubsec:smith_acf>`__ for a
             drill-down inside this routine)

         - | bm_ctl
           | \* (bimodal cloud scheme)

         - Set area cloud fraction equal to bulk cloud fraction

      .. container:: tcolorbox

         | **iau**
         | \* (incremental analysis update; part of data assimilation)

         - | pc2_assim
           | \* (PC2 reponse to the analysis increments; see
             `5.2.4 <#subsubsec:pc2_assim>`__ for a drill-down inside
             this routine)

         - | initial_pc2_check
           | \* (wrapper for optional self-consistency checks on
             prognostic cloud variables if not doing PC2 response to
             analysis increments)

           - | pc2_checks
             | \* (self-consistency checks on cloud fractions and water
               contents)

Drill-downs within some routines in the call tree are listed separately
below, to avoid duplication (since these routines are called in multiple
different places in the tree)...

.. _`subsubsec:smith_acf`:

Smith scheme with area cloud fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container:: itemize

   .. container:: tcolorbox

      | ls_arcld
      | \* (interface to diagnostic Smith scheme and area cloud
        fraction)

      - If no area cloud fraction scheme:

        | ls_cld
        | \* (just directly call Smith scheme)

        Set area cloud fraction equal to bulk cloud fraction.

      - If using Cusack vertical interpolation method:

        Interpolate fields onto finer vertical grid

        | ls_cld
        | \* (call Smith scheme using higher vertical resolution fields)

        Coarse-grain cloud fields back to model grid, but set area cloud
        fraction to max of bulk cloud fraction over corresponding
        fine-grid levels.

      - If using Brooks empirical area cloud fraction method:

        | ls_cld
        | \* (just directly call Smith scheme)

        | ls_acf_brooks
        | \* (estimate area cloud fraction)

.. _`subsubsec:pc2_initiation`:

PC2 initiation
~~~~~~~~~~~~~~

.. container:: itemize

   .. container:: tcolorbox

      | pc2_initiation_ctl
      | \* (interface to PC2 initiation and consistency-checks)

      - | pc2_checks
        | \* (self-consistency checks on cloud fractions and water
          contents)

      - PC2 initiation of liquid-cloud:

        | pc2_bm_initiate
        | \* (using the bimodal cloud scheme)

        | pc2_arcld
        | \* (using the Smith scheme with the Cusack vertical
          interpolation method)

        - | pc2_initiate
          | \* (initiation using the Smith scheme, called here on a
            finer vertical grid as per the Cusack method)

        | pc2_initiate
        | \* (using the Smith scheme with no area cloud representation)

      - | pc2_checks2
        | \* (further self-consistency checks on cloud-fractions)

      - | pc2_checks
        | \* (repeat the first lot of self-consistency checks again,
          just in case we broke something in the mean-time!)

      - | pc2_hom_arcld
        | \* (finds area cloud fraction using a version of the Cusack
          method, where the cloud fraction on the finer vertical grid is
          estimated by applying homogeneous forcing relative to the
          original grid fields)

        - | pc2_homog_plus_turb
          | \* (generic homogeneous forcing routine used to interpolate)

.. _`subsubsec:pc2_assim`:

PC2 Data Assimilation
~~~~~~~~~~~~~~~~~~~~~

.. container:: itemize

   .. container:: tcolorbox

      | pc2_assim
      | \* (PC2 reponse to the analysis increments)

      - | pc2_homog_plus_turb
        | \* (generic PC2 homogeneous forcing routine used here for
          liquid-cloud)

      - Estimate change in ice-cloud fraction from the assimilation
        increment to ice-cloud mass.

      - | pc2_total_cf
        | \* (update bulk cloud fraction due to change in ice cloud
          fraction)

      - | pc2_checks
        | \* (self-consistency checks on prognostic cloud fractions and
          water contents)

.. _`sec:diags`:

Diagnostics
-----------

Nearly all diagnostics retain their meaning when PC2 is run. However,
there are a few that are subtly modified.

The convective diagnostics that use the convective cloud base and top
calculations remain the same if PC2 is used with a zeroed convective
cloud fraction. These values are not reset by the convection scheme,
since the model is still predicting convection between the diagnosed
levels.

The visibility diagnostics need modifying if the convective cloud
fraction is switched off, since they use the convective cloud fraction
within their calculation. Here we use a value of 0.2 for the convective
cloud amount if there is convective precipitation but the
two-dimensional convective cloud amount is zero. This will be the case
if the PC2 scheme has zeroed the convective cloud amount.

There are a number of increment diagnostics that are required to fully
diagnose the moisture cycle within PC2. Since most physics sections can
cause condensation, condensate and cloud fraction increment diagnostics
have been written for each of these sections.

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{C_t}` and
  :math:`\overline{C_l}` increments from SW radiation,
  :math:`\overline{T}` increment from SW Radiation without including the
  condensation: **Section 1** .

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{C_t}` and
  :math:`\overline{C_l}` increments from LW radiation,
  :math:`\overline{T}` increment from LW Radiation without including the
  condensation: **Section 2** .

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}`,
  :math:`\overline{C_t}`, :math:`\overline{C_l}` and
  :math:`\overline{C_f}` increments from Boundary Layer: **Section 3** .

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}`,
  :math:`\overline{C_t}`, :math:`\overline{C_l}` and
  :math:`\overline{C_f}` increments from Large-scale precipitation:
  **Section 4** .

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}`,
  :math:`\overline{C_t}`, :math:`\overline{C_l}` and
  :math:`\overline{C_f}` increments from Convection,
  :math:`\overline{q}`, :math:`\overline{q_{cl}}`,
  :math:`\overline{q_{cf}}`, :math:`\overline{C_t}`,
  :math:`\overline{C_l}` and :math:`\overline{C_f}` increments from the
  inhomogeneous part of the Convection scheme only: **Section 5** .

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}`,
  :math:`\overline{C_t}`, :math:`\overline{C_l}` and
  :math:`\overline{C_f}` increments from the Advection: **Section 12** .

However, there are a number of parts of PC2 that do not fit into a
pre-existing section of code, and hence the associated increment
diagnostics are not easily placed within the UM framework. These
increments were available using a modification set or branch and a
user-STASHmaster file up to version 7.5. From version 7.6 these
diagnostics are available as standard.

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{C_t}`, and
  :math:`\overline{C_l}` increments from the PC2 erosion section:
  **Section 4** or **Section 5** depending on where the erosion is
  called.

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}`,
  :math:`\overline{C_t}`, :math:`\overline{C_l}` and
  :math:`\overline{C_f}` increments from the Bounds Checking after
  atmphya: **Section 4** .

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}`,
  :math:`\overline{C_t}`, :math:`\overline{C_l}` and
  :math:`\overline{C_f}` increments from the Initiation and Bounds
  checking at the end of the timestep: **Section 16** .

-  :math:`\overline{T}`, :math:`\overline{q}`,
  :math:`\overline{q_{cl}}`, :math:`\overline{q_{cf}}`,
  :math:`\overline{C_t}`, :math:`\overline{C_l}` and
  :math:`\overline{C_f}` increments from the Pressure Forcing section:
  **Section 16** .

Single Column Model
-------------------

The updating in the single column model follows the same timestepping as
that in the full model, but the changes to atm-step are mirrored within
scm_main. The method used is to store the driving SCM forcing increments
of vapour, liquid and temperature across the forcing subroutine. The
forcing of pressure is set to zero. After atmos_physics2 has been
called, a PC2 section of code calls the homogeneous forcing subroutine
with these increments. This therefore treats the response of PC2 to the
prescribed dynamical forcing in the SCM as homogeneous. Following this
calculation, the initiation scheme is called, as usual. Finally, the
area cloud fraction is set to the bulk cloud fraction and
:math:`\overline{\Theta}` (potential temperature) is made consistent
with :math:`\overline{T}` (dry-bulb temperature) which was changed by
the condensation in the PC2 response to the homogeneous forcing. The
rest of the SCM uses the same PC2 code as the full model.

Note that the change to PC2 homogeneous forcing from advection under the
UM namelist switch **l_pc2_sl_advection** (see section
`4.8 <#sec:pres>`__) is also mirrored in the Single-Column Model. If
this switch is turned on, the PC2 homogeneous forcing call using the SCM
forcing increments is moved straight after the call to the forcing
routine, so that the condensation adjustment is performed before the
call to atmos_physics2. If **l_pc2_sl_advection** is turned on, the PC2
response to SCM forcings is also improved as follows...

The SCM forcings may comprise one or both of the following:

- (a) Prescribed tendencies or relaxation applied to T,q.

- (b) Interactive vertical advection applied to T,q.

For the latter, we can calculate the pressure change experienced by
vertically-advected parcels, and so calculate the PC2 homogeneous
forcing response in the same way as we do for Semi-Lagrangian advection
in the full model (see section `4.8 <#sec:pres>`__). For the former, we
don’t know if the prescribed T,q tendencies are due to advection,
radiation, or some other process, so we calculate the PC2 homogeneous
forcing response as if the tendencies are applied "in-situ".

To split the PC2 homogeneous response into these 2 components, the SCM
forcing routine outputs:

- (a) The forcing increments to T,q excluding the contribution from
  interactive vertical advection.

- (b) The value of exner pressure at departure points, consistent with
  the vertical advection.

The PC2 homogeneous forcing responses to these 2 forcing components are
then calculated by 2 separate PC2 calls in scm_main.

Limited Area Boundary Conditions
--------------------------------

Cloud fractions on the limited area boundaries are fully updateable.
Writing of cloud fraction Limited Area Boundary Conditions (LBCs) will
be automatic if PC2 is selected. A PC2 LAM may be run from an LBC file
with or without cloud fraction LBCs (this is specified by the logical
l-pc2-lbc, which is set in the UMUI). If there are no cloud fraction
lbcs then around the edge of the domain the checking and initiation
routines will be applying significant increments to the cloud and
condensate fields near the boundaries, but this does not have an adverse
effect well away from the boundaries. If there are no cloud fraction
LBCs the cloud fraction fields themselves are not forced to zero around
the edge of the domain but are allowed to freely find their own value. A
PC2 run that outputs lbcs will, by default, always output cloud
fractions as part of the LBCs file.

Parameter values
----------------

Table `1 <#tab:pc2_names>`__ summarizes the values of parameters used in
the PC2 scheme and their location within various comdecks. Those
parameters marked as ‘Num’ are those that are not part of the
mathematical equation set that is being solved, but are required in
order to achieve a stable, realistic, numerical solution. These include,
for instance, thresholds for resetting cloud fractions back to 0 or 1.
Those marked ’Phy’ are physical quantities that form an integral part of
the equation set that we wish to solve. Those marked ’Clo’ form part of
a closure needed to form the equation set, but are less readily related
to physical quantities. Variables marked ’Diag’ form a part of the
diagnostic output routines.

.. container:: center

   .. container::
      :name: tab:pc2_names

      .. table:: PC2 parameter values and locations

         +----------+----------+----------+----------+----------+----------+
         | Symbol   | Code     | Des      | Value    | Location | Notes    |
         |          | variable | cription |          |          | and ref. |
         +==========+==========+==========+==========+==========+==========+
         | -        | init-it  | Number   | 10       | p        | Num:     |
         |          | erations | of       |          | c2-const | `3.4     |
         |          |          | it       |          |          | .2 <#sec |
         |          |          | erations |          |          | :numapp_ |
         |          |          | in       |          |          | init>`__ |
         |          |          | in       |          |          |          |
         |          |          | itiation |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :math:`  | cloud    | Bounds   | 0.005    | UM       | Num:     |
         | C_{tol}` | -pc2-tol | checking |          | namelist | `4.9     |
         |          |          | :ma      |          |          |  <#sec:i |
         |          |          | th:`C_l` |          |          | nit2>`__ |
         |          |          | t        |          |          |          |
         |          |          | hreshold |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :        | cloud-p  | Bounds   | 0.001    | UM       | Num:     |
         | math:`C_ | c2-tol-2 | checking |          | namelist | `4.9     |
         | {tol 2}` |          | :ma      |          |          |  <#sec:i |
         |          |          | th:`C_l` |          |          | nit2>`__ |
         |          |          | t        |          |          |          |
         |          |          | hreshold |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :math:`R | rh       | :        | 0.01     | p        | Num:     |
         | H_{tol}` | crit-tol | math:`RH |          | c2-const | `4.9     |
         |          |          | _{crit}` |          |          |  <#sec:i |
         |          |          | t        |          |          | nit2>`__ |
         |          |          | olerance |          |          |          |
         |          |          | in       |          |          |          |
         |          |          | in       |          |          |          |
         |          |          | itiation |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :math    | ls-bl0   | Fixed    | :ma      | imp-ctl  | Clo:     |
         | :`q_{cf0 |          | value of | th:`1.0  |          | `        |
         |  \, BL}` |          | BL       | \times 1 |          | 4.6 <#se |
         |          |          | in-plume | 0^{-4} \ |          | c:bl>`__ |
         |          |          | :        | , kg \,  |          |          |
         |          |          | math:`\o | kg^{-1}` |          |          |
         |          |          | verline{ |          |          |          |
         |          |          | q_{cf}}` |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :math:`  | one-     | Fixed    | :ma      | pc2-chck | Num:     |
         | q_{cf0}` | over-qcf | in-cloud | th:`1.0  |          | `4.10    |
         |          |          | :        | \times 1 |          | <#sec:ch |
         |          |          | math:`\o | 0^{-4} \ |          | ecks>`__ |
         |          |          | verline{ | , kg \,  |          |          |
         |          |          | q_{cf}}` | kg^{-1}` |          |          |
         |          |          | if       |          |          |          |
         |          |          | :math:`  |          |          |          |
         |          |          | C_f`\ =0 |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :        | pdf-mer  | Merging  | 0.5      | p        | Clo:     |
         | math:`m` | ge-power | power    |          | c2-const | `3.2     |
         |          |          | for      |          |          |  <#sec:h |
         |          |          | :math:`  |          |          | omog>`__ |
         |          |          | G(-Q_c)` |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :        | p        | Shape    | 0.0      | p        | Phy:     |
         | math:`n` | df-power | p        |          | c2-const | `3.2     |
         |          |          | arameter |          |          |  <#sec:h |
         |          |          | for      |          |          | omog>`__ |
         |          |          | :math:`  |          |          |          |
         |          |          | G(-Q_c)` |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :        | w        | Wind     | :mat     | p        | Phy:     |
         | math:`w` | ind-shea | shear in | h:`1.5 \ | c2-const | `        |
         |          | r-factor | fallout  | times 10 |          | 4.2.1 <# |
         |          |          | of ice   | ^{-4} \, |          | sec:lsp_ |
         |          |          | term     |  s^{-1}` |          | fall>`__ |
         +----------+----------+----------+----------+----------+----------+
         | :        | i        | Scaling  | 0.04     | p        | Phy:     |
         | math:`i` | ce-width | factor   |          | c2-const | `4       |
         |          |          | for      |          |          | .2.4 <#s |
         |          |          | r        |          |          | ec:mp_de |
         |          |          | eduction |          |          | psub>`__ |
         |          |          | in       |          |          |          |
         |          |          | :ma      |          |          |          |
         |          |          | th:`b_i` |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :        | dbsdtb   | Rate of  | :math:   | UM       | Phy:     |
         | math:`a` | s-turb-0 | r        | `-2.25 \ | namelist | `3.3     |
         |          |          | eduction | times 10 |          |  <#sec:w |
         |          |          | of PDF   | ^{-5} \, |          | idth>`__ |
         |          |          | width    |  s^{-1}` |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :        | dbsdtb   | Rate of  | 0        | p        | Phy:     |
         | math:`b` | s-turb-1 | r        |          | c2-const | `3.3     |
         |          |          | eduction |          |          |  <#sec:w |
         |          |          | of PDF   |          |          | idth>`__ |
         |          |          | width    |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         |          | dbsd     | Redn of  | 0        | p        | Phy:     |
         |          | tbs-conv | PDF      |          | c2-const | `3.3     |
         |          |          | width in |          |          |  <#sec:w |
         |          |          | co       |          |          | idth>`__ |
         |          |          | nvection |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         |          | dbs      | V        | 10.05    | p        | Phy:     |
         |          | dtbs-exp | ariation |          | c2-const | `3.3     |
         |          |          | of       |          |          |  <#sec:w |
         |          |          | erosion  |          |          | idth>`__ |
         |          |          | on RH    |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :        | RHCRIT   | Critical |          | UM       | Phy:     |
         | math:`RH |          | RH for   |          | namelist | `3.4     |
         | _{crit}` |          | cloud    |          |          |  <#sec:i |
         |          |          | f        |          |          | nit>`__, |
         |          |          | ormation |          |          | `4       |
         |          |          |          |          |          | .2.4 <#s |
         |          |          |          |          |          | ec:mp_de |
         |          |          |          |          |          | psub>`__ |
         +----------+----------+----------+----------+----------+----------+
         | :math:   | condensa | Minimum  | :m       | pc2-chck | Num:     |
         | `q_{c0}` | te-limit | allowed  | ath:`1 \ |          | `4.10    |
         |          |          | co       | times 10 |          | <#sec:ch |
         |          |          | ndensate | ^{-10} \ |          | ecks>`__ |
         |          |          |          | , kg \,  |          |          |
         |          |          |          | kg^{-1}` |          |          |
         +----------+----------+----------+----------+----------+----------+
         | :math:`q | ls0      | Lower    | :        | enviro?a | Num:     |
         | _c^{S0}` |          | limit of | math:`5  |          | `3.5.    |
         |          |          | plume    | \times 1 |          | 2 <#sec: |
         |          |          | co       | 0^{-5} \ |          | multi_nu |
         |          |          | ndensate | , kg \,  |          | mapp>`__ |
         |          |          |          | kg^{-1}` |          |          |
         +----------+----------+----------+----------+----------+----------+
         |          | *Har     | Conv     | 0.2      | imp-ctl2 | Diag:    |
         |          | d-wired* | cloud    |          |          | `5.3     |
         |          |          | fraction |          |          |  <#sec:d |
         |          |          | for      |          |          | iags>`__ |
         |          |          | vi       |          |          |          |
         |          |          | sibility |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         |          | *Har     | Limit on | 0.001    | lspice3d | Num:     |
         |          | d-wired* | width of |          |          | `4       |
         |          |          | ice      |          |          | .2.4 <#s |
         |          |          | dist     |          |          | ec:mp_de |
         |          |          | ribution |          |          | psub>`__ |
         +----------+----------+----------+----------+----------+----------+
         |          | *Har     | :ma      | 0.05     | pc2-init | Num:     |
         |          | d-wired* | th:`C_l` |          |          | `4.9     |
         |          |          | limit    |          |          |  <#sec:i |
         |          |          | for init |          |          | nit2>`__ |
         |          |          | if       |          |          |          |
         |          |          | :math:`T |          |          |          |
         |          |          |  < 0 ^{\ |          |          |          |
         |          |          | circ} C` |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         |          | *Har     | T        | :mat     | imp-ctl  | Num:     |
         |          | d-wired* | olerance | h:`1.0 \ |          | `        |
         |          |          | on calc. | times 10 |          | 4.6 <#se |
         |          |          | of       | ^{-10} \ |          | c:bl>`__ |
         |          |          | :math    | , kg \,  |          |          |
         |          |          | :`q_C^s` | kg^{-1}` |          |          |
         |          |          | in BL    |          |          |          |
         +----------+----------+----------+----------+----------+----------+

PC2 also recommends some tunings of the existing convection scheme
parameters. These cannot be placed in the library code, since they would
interact with non-PC2 simulations, hence would need to be specified with
modification sets. We have included those parameters that have been
investigated throughout testing, although only two are different between
PC2:64 and a non-PC2 run.

.. container:: center

   .. container::
      :name: tab:pc2_conv_names

      .. table:: PC2 parameter values and locations relating to the
      convection. \*These values are those used in HadGAM

         +----------+----------+----------+----------+----------+----------+
         | Code     | Des      | Value in | Value in | Location | Notes    |
         | variable | cription | PC2      | Control  |          | and      |
         |          |          |          |          |          | r        |
         |          |          |          |          |          | eference |
         +==========+==========+==========+==========+==========+==========+
         | TICE     | Tem      | :math:   | :math:`0 | tice.cdk | Phy:     |
         |          | perature | `-10 ^{\ |  ^{\circ | or UMUI  | `4.7     |
         |          | at which | circ} C` | } C`\ \* |          | <#sec:co |
         |          | plume    |          |          |          | nvec>`__ |
         |          | freezes  |          |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | QSTICE   | App      | :m       | :m       | qs       | Phy:     |
         |          | roximate | ath:`3.5 | ath:`3.5 | tice.cdk | `4.7     |
         |          | qs       |  \times  |  \times  | or UMUI  | <#sec:co |
         |          | at(TICE) | 10^{-3}` | 10^{-3}` |          | nvec>`__ |
         +----------+----------+----------+----------+----------+----------+
         | *Har     | Limit on | 0.5      | :math:   | cloudw   | Phy:     |
         | d-wired* | conv.    | :        | `0.5 \,  |          | `4.7     |
         |          | cond.    | math:`q_ | q_{sat}` |          | <#sec:co |
         |          | after    | {sat}, 2 |          |          | nvec>`__ |
         |          | precip   |  \times  |          |          |          |
         |          |          | 10^{-4}` |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | Anvil    | Shape    | 0        | 0.3\*    | UMUI     | Phy:     |
         | factor   | p        |          |          |          | `4.7     |
         |          | arameter |          |          |          | <#sec:co |
         |          | for      |          |          |          | nvec>`__ |
         |          | conv.    |          |          |          |          |
         |          | cloud    |          |          |          |          |
         |          | anvil    |          |          |          |          |
         +----------+----------+----------+----------+----------+----------+
         | Tower    | Shape    | 0        | 0.25\*   | UMUI     | Phy:     |
         | factor   | p        |          |          |          | `4.7     |
         |          | arameter |          |          |          | <#sec:co |
         |          | for      |          |          |          | nvec>`__ |
         |          | conv.    |          |          |          |          |
         |          | cloud    |          |          |          |          |
         |          | tower    |          |          |          |          |
         +----------+----------+----------+----------+----------+----------+

How to run the PC2 scheme
-------------------------

Running PC2 is straightforward, but you should seek advice as to
modification sets that you need to include to ensure you are running the
most up-to-date version of PC2. The following is a brief checklist of
the options in the UMUI which need to be selected in order to run PC2.
No hand-edits are required.

- In the LS cloud panel (atmos-science-section-LScloud) push the button
  marked ’use the PC2 cloud scheme’.

- If you wish to use PC2 in the diagnostic only mode, also push ’run the
  PC2 scheme in diagnostic only mode’. If you wish to run PC2 fully then
  do not push this button

- In the large-scale precipitation section
  (atmos-science-section-LSprecip) select the 3D large-scale
  precipitation scheme.

- The specification of the LA boundary conditions can be set in the
  atmos-InFiles-OtherAncil-LBC panel.

- You will need to select modsets to include update the library code to
  the PC2 version described here. Seek advice on these.

- You may wish to adjust the convective anvil parameters in
  atmos-science-section-convec. Again, seek advice.

More information
----------------

Information on results of the scheme and how to run the PC2 code at
various model versions is available on the PC2 web site.

.. _`sec:appendix-da`:

Appendix: Alternative PC2 - Data Assimilation formulations
==========================================================

In this alternative method to section `4.11 <#sec:da>`__ we will assume
that there exists a homogeneous forcing, :math:`\Delta Q_c`, that gives
changes, net of condensation, of :math:`\Delta\overline{q}` and
:math:`\Delta\overline{T}`. If we can recover what :math:`\Delta Q_c` is
then we can use this to calculate the liquid, :math:`\overline{q_{cl}}`,
and liquid cloud fraction, :math:`C_l`, increments.

As in section `4.11 <#sec:da>`__, we start by discretising
(`[dqcldt] <#dqcldt>`__) to give

.. math::

   \Delta \overline{q_{cl}} = C_l \Delta Q_c
   \label{eq:dqcldt_discrete}

and hence, using the discrete form of :math:`\Delta Q_c` from
(`[eq:deltaqc_exp2] <#eq:deltaqc_exp2>`__) gives

.. math::

   \Delta \overline{q_{cl}} = C_l ( a_L ( \Delta \overline{q} - \alpha \Delta 
   \overline{T} ) + \Delta \overline{q_{cl}} ) ,

which rearranges to

.. math::

   \Delta \overline{q_{cl}} = \frac{1}{1-C_l} C_l a_L ( \Delta \overline{q}
   - \alpha \Delta \overline{T} - \beta \Delta \overline{p}) .
   \label{eqn:delataqcl}

Comparing to (`[eq:deltaqc_exp2] <#eq:deltaqc_exp2>`__) and
(`[eq:dqcldt_discrete] <#eq:dqcldt_discrete>`__) we see that
:math:`\Delta \overline{q_{cl}}` is the same as if we had applied the
homogeneous forcing technique using :math:`\Delta \overline{q}`,
:math:`\Delta \overline{T}` and :math:`\Delta \overline{p}` as forcings,
except multiplied by a factor of :math:`\frac{1}{1-C_l}`.

We can calculate :math:`\Delta C` in a similar way. From
(`[eq:deltac] <#eq:deltac>`__)

.. math:: \Delta C_l = G(-Q_c) \Delta Q_c

and hence, using our value of :math:`\Delta Q_c` from
(`[eq:deltaqc_exp2] <#eq:deltaqc_exp2>`__) and
:math:`\Delta \overline{q_{cl}}` from
(`[eqn:delataqcl] <#eqn:delataqcl>`__)

.. math::

   \Delta C_l = G(-Q_c) (a_L ( \Delta \overline{q} - \alpha \Delta \overline{T} )
   + \frac{1}{1-C_l} C_l a_L ( \Delta \overline{q} - \alpha \Delta \overline{T} ) )

which rearranges to

.. math::

   \Delta C_l = \frac{1}{1-C_l} G(-Q_c) a_L ( \Delta \overline{q}
   - \alpha \Delta \overline{T} -\beta \Delta \overline{p}) .
   \label{eqn:c1mc}

This is also a factor of :math:`\frac{1}{1-C_l}` different from using
:math:`\Delta \overline{q}`, :math:`\Delta \overline{T}` and
:math:`\Delta \overline{p}` directly as forcings (the factor must be the
same, as we are still using the homogeneous forcing hypothesis). This
equation forms the basis for the more advanced technique discussed in
this section. However, it is undefined at :math:`C_l=1` and becomes
ill-conditioned near :math:`C_l=1`, hence there must be care taken when
this expression is solved numerically.

Numerical solution
------------------

The timestepping applied is picked as a result of numerical tests
forcing a single gridbox with uniform increments. Many numerical
techniques were tested, this gives a fast but reasonably well behaved
solution.

Initially, we calculate :math:`G(-Qc)` and :math:`\Delta Q_c` from the
input fields, as in the homogeneous forcing technique (section
`3.2 <#sec:homog>`__) and (`[eq:deltaqc_exp2] <#eq:deltaqc_exp2>`__).

An initial increment, :math:`\Delta C_l^1` is estimated directly using
the basic equation

.. math:: \Delta C_l^1 = \frac{1}{1-C_l^(n)} G(-Q_c) \Delta Q_c .

We then recalculate this expression, using a mid-timestep estimate for
:math:`\Delta C_l`;

.. math::

   \Delta C_l =  \frac{1}{1-(C_l^{[n]} + \frac{1}{2} \Delta C_l^1)}
   G(-Q_c) \Delta Q_c

where the term :math:`C_l^{[n]} + \frac{1}{2} \Delta C_l^1` is limited
to be no more than 0.9999 to avoid divide by zero problems. The final,
updated value of cloud fraction, :math:`C_l^{[n+1]}`, is then

.. math:: C_l^{[n+1]} = C_l^{[n]} + \Delta C_l

and this value is limited to 0 or 1.

The liquid water term simply uses the final version of :math:`C_l` in
its calculation.

.. math:: \Delta \overline{q_{cl}} = \frac{1}{1-C_l^{[n+1]}} C_l^{[n+1]} \Delta Q_c

and will be set to 0 if :math:`C^{[n+1]}` is 0. There is an additional
limit, see below, applied to the liquid water term, which will prevent
the value of :math:`\Delta \overline{q_{cl}}` increasing to a large
number if :math:`C_l^{[n+1]}` is very close to 1.

Limit on the liquid water content
---------------------------------

We will choose a limit on :math:`\overline{q_{cl}}` to be equal to its
value when the underlying PDF just corresponds to total cloud cover.
Therefore, from (`[eq:qclbar=int] <#eq:qclbar=int>`__)

.. math:: \overline{q_{cl \, max}} = \int_{s=-b_s}^{\infty} G(s) (b_s + s) ds .

We will use the current value of :math:`Q_c` (which won’t in general to
be equal to :math:`b_s`) to split the integral into two ranges of s:

.. math::

   \overline{q_{cl \, max}} = \int_{s=-b_s}^{-Q_c} G(s) (b_s + s) ds
   + \int_{s=-Q_c}^{\infty} G(s) (b_s + s) ds .

For the moment we write the first of these integrals as :math:`I1`, and
split the second integral whilst introducing a :math:`(+ Q_c - Q_c)`
term to the integrand:

.. math::

   \overline{q_{cl \, max}} = I1 + \int_{s=-Q_c}^{\infty} G(s) (b_s - Q_c) ds
   + \int_{s=-Q_c}^{\infty} G(s) (s + Q_c) ds .

The last of the integrals is now the current liquid water content,
:math:`\overline{q_{cl}}`. The second integral is proportional to the
liquid cloud fraction :math:`C_l`.

.. math:: \overline{q_{cl \, max}} = I1 + C_l (b_s - Q_c) + \overline{q_{cl}}

or

.. math::

   \overline{\Delta q_{cl \, max}} = I1 + C_l (b_s - Q_c) .
   \label{eqn:deltaqclmax}

Now consider the expression for the saturation deficit, which we have
defined, from (`[SD] <#SD>`__) as

.. math:: SD = \int_{-b_s}^{-Q_c} G(s) (-Q_c - s) ds .

Splitting and adding the term :math:`(+b_s - b_s)` to the integrand in a
similar way to above gives

.. math::

   \begin{aligned}
   SD = \int_{-b_s}^{-Q_c} G(s) (-Q_c + b_s) ds + \int_{-b_s}^{-Q_c}
   G(s) (-s - b_s) ds \nonumber \\
   = (-Q_c + b_s) (1 - C_l) - I1 ,
   \end{aligned}

and hence :math:`I1` in terms of :math:`SD`. Using this value of
:math:`I1` in (`[eqn:deltaqclmax] <#eqn:deltaqclmax>`__) and cancelling
the :math:`C_l` terms gives :math:`\Delta \overline{q_{cl \, max}}` as

.. math::

   \Delta \overline{q_{cl \, max}} = (-Q_c + b_s) - SD .
   \label{eqn:delta2}

This is a general expression, it is not fixed for a particular PDF. To
complete the analysis, we need to estimate :math:`-Q_c+b_s`. To do this,
we now make the *assumption* of a power-law type PDF, as in section
`3.4 <#sec:init>`__. If we start from the equivalent of
(`[eqn19] <#eqn19>`__) but at the :math:`s=-bs` end of the distribution,
equation (B.3) in :raw-latex:`\cite{wg03}` can be equivalently written
for :math:`(1-C_l)` as:

.. math::

   (1-C_l) = \frac{ A (-Q_c + b_s)^{n+1} }{n+1} .
   \label{eqn:1mc}

To derive this from (B.3) note that :math:`C_l` is swapped for
:math:`1-C_l` and :math:`(b_s - (-Q_c))` is swapped for
:math:`(-Qc - (-b_s))`, as in section `3.4.2 <#sec:numapp_init>`__.
Similarly, noting that :math:`\overline{q_{cl}}` can be swapped with
:math:`SD`, gives the equivalent to (B.4) in :raw-latex:`\cite{wg03}` as

.. math::

   SD = \frac{ A (-Q_c + b_s)^{n+2} }{(n+1)(n+2)}.
   \label{eqn:sd}

Using the value :math:`(1-C_l)` from (`[eqn:1mc] <#eqn:1mc>`__) in
(`[eqn:sd] <#eqn:sd>`__) gives

.. math:: \frac{SD}{1-C_l} = \frac {-Q_c + b_s}{n+2} .

Finally, we use this expression for :math:`(-Q_c + b_s)` in
(`[eqn:delta2] <#eqn:delta2>`__) to parametrize
:math:`\Delta \overline{q_{cl \, max}}` in terms of the saturation
deficit

.. math::

   \Delta \overline{q_{cl \, max}} = SD ( \frac{n+2}{1-C_l} - 1 ) .
   \label{eqn:sdr1mc}

This is the expression that is used for the limit on
:math:`\overline{q_{cl}}`. We subsequently apply a second limit, since
numerically this expression is still not well behaved when :math:`C_l`
is close to 1. Here we note that just at complete cloud cover for a
symmetric PDF we have :math:`\overline{q_{cl}} = b_s`. Hence we estimate
:math:`b_s` as in :raw-latex:`\cite{smith90}`,

.. math::

   b_s = a_L ( 1 - RH_{crit} ) q_{sat}(\overline{T_L}) ,
   \label{eqn:bs}

and take the smaller value for of (`[eqn:sdr1mc] <#eqn:sdr1mc>`__) and
(`[eqn:bs] <#eqn:bs>`__) for :math:`\Delta \overline{q_{cl \, max}}`.

Initiation from :math:`C_l=1`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The equations are not defined when :math:`C_l=1`. (Note that when
:math:`C_l=0` we will calculate :math:`G(-Q_c)=0` so there is no change
in cloud fraction or liquid water content in this case). The
assimilation is capable of lowering :math:`\overline{q}` below
:math:`q_{sat}(\overline{T})` and hence there should be a corresponding
change in cloud fraction and liquid water content. In theory, we can use
the expression for :math:`\Delta \overline{q_{cl \, max}}` and assume
that the initial liquid water is equal to :math:`b_s`. However, this
produces a tricky set of simulataneous equations, which are not easily
solvable except in the case where :math:`n=0`. We proceed by making this
assumption for :math:`n`, acknowledging that this is not necessarily
entirely consistent with the rest of the model (although it is in the
PC2:64 formulation).

We have (equivalent to B.6 from :raw-latex:`\cite{wg03}`)

.. math:: \frac{ (1-C_l)^2 }{SD} = G(-Q_c) \frac{n+2}{n+1}.

If n=0 (i.e. a ‘top-hat’ function) then
:math:`G(-Q_c) = \frac{1}{2 b_s}` and we can write

.. math:: C_l = 1 - \sqrt{ \frac{SD}{b_s} } .

We now assume :math:`b_s` is equal to our current value of
:math:`\overline{q_{cl}}` and hence

.. math::

   C_l^{[n+1]} = 1 - \sqrt{ \frac{SD^{[n+1]}}{\overline{q_{cl}^{[n]}}} }
   \label{eqn:1msqrt}

where :math:`C_l^{[n+1]}` and :math:`SD^{[n+1]}` are the values of
:math:`C_l` and :math:`SD` after this initiation has been applied. Using
our previous expression (`[eqn:sdr1mc] <#eqn:sdr1mc>`__) for
:math:`\Delta \overline{q_{cl max}}` gives (remembering that we are
considering the reverse process, so the sign is opposite),

.. math:: \Delta \overline{q_{cl}} = - SD^{[n+1]}  ( \frac{2}{1-C_l^{[n+1]}} - 1 )

(remembering that :math:`n=0` is assumed). Hence, replacing
:math:`C_l^{[n+1]}` by (`[eqn:1msqrt] <#eqn:1msqrt>`__) we have

.. math::

   \Delta \overline{q_{cl}} = SD^{[n+1]} - 2 \sqrt{ SD^{[n+1]}
   \overline{q_{cl}}^{[n]} } .

This is the expression we use, :math:`SD^{[n+1]}` is calculated after
the ssimilation increments have been applied, using (`[SD2] <#SD2>`__):

.. math::

   SD^{(n+1)} = a_L^{[n+1]} ( q_{sat}(\overline{T}^{[n+1]},
   \overline{p}^{[n+1]}) - \overline{q}^{[n+1]} ).

Results
-------

Results demonstrate a problem in that there is a distinct asymmetry
between changes when :math:`\Delta Q_c` is large and positive and when
:math:`\Delta Q_c` is large and negative, when cloud fractions start
near 1. In the former case, the limit to the amount of liquid and cloud
fraction that can be created means that changes must be kept relatively
small, whereas in the latter case, all the cloud and liquid water can be
removed easily. (The :math:`1/(1-C_l)` term allows this to be done
relatively quickly). Hence this assimilation method has a net tendency
to remove cloud from the simulation, which, at the moment, gives poorer
results than simply using the homogeneous forcing method.

Further work will be required to enable the implementation of this
:math:`\overline{q}` and :math:`\overline{T}` preserving method.

.. _`sec:code-development`:

Appendix: Essentials of PC2 for code developers
===============================================

This section provides some guidance to code developers on the treatment
of PC2. Code developers are advised to read the relevant part of section
`4 <#sec:app_um>`__ to understand the way in which the current PC2
scheme interacts with their section of code.

The essence of a prognostic cloud scheme is that each physical part of
the model is able to calculate increments to the cloud fractions and
condensate contents. These form an integral part of each physics scheme
and should be considered by code owners as such, hence any alteration to
a scheme *must* consider also the impact on :math:`q_{cl}`,
:math:`q_{cf}`, :math:`C_t`, :math:`C_l` and :math:`C_f`, as well as on
the more traditional :math:`T`, :math:`q` and wind prognostics. Often
there should be no impact, but this cannot be assumed without
consideration. There is no diagnostic cloud fraction and condensation
scheme which can be run in PC2, since this would reset any effect of the
cloud prognostics used elsewhere in the model. (The diagnostic scheme
can still be used for model *diagnostics*, such as visibility and fog
fraction, and will be kept in later versions of the UM).

Since this places a significant burden on code developers, the PC2
developers have produced two generic representations which can take
increments to :math:`q` and :math:`T` etc. and produce an estimate of
the condensation and cloud fraction changes associated with the
increments. These are the homogeneous forcing and injection forcing (or
inhomogeneous forcing) methods.

Homogeneous forcing
-------------------

This is described fully in section `3.2 <#sec:homog>`__. This assumes
that the distribution of :math:`q_T - q_{sat}(T_L)` about its gridbox
mean is unchanged when a process acts. (The mean will change of course,
but we assume that the variations in each part of the gridbox from the
mean do not). Since this is equivalent to every part of the gridbox
receiving the same :math:`q_T` and :math:`T_L` increment, we call this
‘Homogeneous Forcing’. We have provided a subroutine
*pc2-homog-plus-turb*, in deck *pc2-homo* in order to provide the
necessary updates.

Injection forcing
-----------------

This is described fully in section `3.5 <#sec:inhomog>`__. We assume
that we already know a condensate increment :math:`q_{cl}` or
:math:`q_{cf}` and that a corresponding cloud fraction increment
:math:`C_l` or :math:`C_f` (and :math:`C_t`) remains to be estimated.
The injection forcing assumes that new cloud randomly displaces existing
cloud in a gridbox, and is designed with detrainment from deep
convection in mind, although it is also used elsewhere. It will require
as an input an estimate of the ‘in-cloud’ water content of the new cloud
that is produced.

If you consider that both the homogeneous and injection forcing
representations are both poor assumptions for your scheme, you will need
to provide another method for calculating the condensation and cloud
fraction changes. The PC2 team can advise, but you should not expect
them to do the work. You can, of course, replace existing homogeneous
and inhomogeneous forcing calls with new representations of changes to
the prognostics if you think you have improved representations
available. This is part of the development of any prognostic variable
representation.

Do I need to modify anything when I change a parametrization scheme?
--------------------------------------------------------------------

Here we assume that you wish to do the minimum work possible to get PC2
to work, rather than a full reconsideration of the physics of the PC2
increment terms.

If your scheme is currently using the homogeneous forcing then there is
no need to update the cloud part of the scheme, *provided that you do
not alter values of :math:`T` and :math:`q` after the homogeneous
forcing section is called* and that the physical interpretation of your
:math:`q` and :math:`T` increments does not change. You need to be
careful if you are moving code from one subroutine to another that you
don’t inadvertently do this, although the forcing usually sits at the
end of the control subroutine.

If your scheme is currently using the injection forcing *subroutine*,
which necessitates that condensate increments are already calculated by
the scheme, then there is also no need to update the cloud part of the
scheme. This currently applies to the boundary layer, where
:math:`q_{cf}` is altered by tracer mixing. Like for the homogeneous
schemes, this is provided that you *do not alter :math:`T`, :math:`q` or
condensate values after the injection forcing subroutine is called* and
that the physical interpretation of your :math:`q` and :math:`T`
increments does not change.

Changes to winds do *not* need to have a condensation or cloud fraction
increment associated with them. There may be future scope for developing
an orographic cloud representation (probably diagnostic), but this is
not an essential part of the scheme as it stands.

If your scheme uses hardwired assumptions about what is happening e.g.
convection or microphysics, then you *do* need to be careful that
:math:`T`, :math:`q` and condensates are still calculated correctly
after you have performed your changes. Currently there are many PC2
assumptions hard-wired into the mass-flux convection scheme:

- Any change to the scientific basis by which changes to :math:`T`,
  :math:`q`, :math:`q_{cl}` and :math:`q_{cf}` are calculated requires
  careful consideration

- Simple changes to convective parameters, such as detrainment rates,
  should not require a change to the PC2 code

- Be particularly careful if you move code around, *especially the
  calculation of convective cloud fractions*, since PC2 incorporates a
  set-to-zero in the code. This will need to be replicated or there is a
  risk that the diagnostic cloud fraction is no longer set to zero
  correctly by PC2.

Each microphysics transfer term has been considered individually for PC2
and this should remain the case.

Be especially careful when you do anything in the atmphy and atmstep
levels of the code that includes additional changes :math:`T`,
:math:`q`, :math:`q_{cl}` or :math:`q_{cf}`, since they may need cloud
fraction or condensation changes to go along with them.

In summary, changes to existing increments of :math:`T`, :math:`q` etc.
within the current UM structure are unlikely to necessitate a
modification for PC2 if their physical interpretation has not changed.
However, new methods of generating :math:`T` and :math:`q` increments
will require new code to be added for PC2.

Further PC2 development work
----------------------------

There are a number of areas in which the PC2:66 formulation can be
developed further, and many of these have been mentioned in the
documentation above. Some of these are simple sensitivity studies which
have not been fully explored in development, others are more complex
alterations. It is fair to say that the link to the convection has
proved the most problematic issue so far with PC2 development.

PC2 cloud erosion
~~~~~~~~~~~~~~~~~

The cloud erosion is a critical term for the simulation of shallow
convective cloud. A large amount of erosion is required to keep the
cloud fractions relatively low in shallow convection, which is why we
have linked the erosion to the relative humidity. We recognise, however,
that this is more an empirical choice than a physically informed choice.
In particular, a low relative humidity (e.g. in the stratosphere) would
imply a very high erosion rate - although the net effect is to remove
any cloud, which is a reasonable thing to do, there is an implication of
the parametrization that mixing within the stratosphere is high, which
is clearly incorrect. We have also seen relatively low cloud fractions
in the mid-levels of deep convection in PC2, and presume that this is
influenced by the erosion formulation. A link to mass flux has also been
proposed, but tests with CRMs do not support a clear link. Perhaps it is
more natural to compare the erosion with the turbulent kinetic energy.
This should be available within the boundary layer and convection
schemes, but not outside of these in the current UM.

The erosion formulation in PC2:66 is one where the width of the PDF is
always narrowed (developed following :raw-latex:`\cite{sg03}`). It may
be advantageous to think whether there are unmodelled processes in the
atmosphere that result in an increase in width. Clearly convection is
likely to be one, but this is already represented in PC2. There may be
other models entirely for the way in which the PDF changes as a result
of mixing of air within a gridbox or within the column, these may prove
fruitful to explore.

Another issue is whether width-narrowing (or widening) is really an
effective way of representing the erosion process. CRM evidence suggests
that the required erosion rates to balance convective cloud generation
are larger for liquid cloud fraction than liquid water (by up to a
factor of 2), suggesting that the real atmospheric erosion favours
removal of cloud fraction over liquid water more strongly than the
model.

The in-cloud condensate that is detrained from convective plumes is
high. We might think that the mixing in of environmental air in reality
is likely to lead to more cloud around the plumes and lower condensate
within the plumes. However, the width narrowing scheme is not a good
model of mixing in this situation, always reducing the amount of cloud
because it is incorrectly assumed that much of the detrained plume has
condensate contents only just above zero and that the shape of the
moisture PDF remains unchanged. This may have a bearing on the problem
of the lack of mid-level cloud in the model (although I think there are
many reasons for this). A different mixing method may give significantly
different results for the areas around convective plumes.

Narrowing of the moisture PDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the parametrized terms in PC2 act to reduce the width of the
moisture PDF. The only terms that can increase the width are the
convection, and the initiation (which can reset the width). This may not
be the best way to describe the way in which the PDF evolves, in
particular it is sensible to ask whether the erosion term should
actually increase the width in the presence of large vertical gradients
of moisture.

Convective cloud increments in the mass-flux framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As discussed in section `3.5.3 <#sec:conv_imp_note>`__, it would be
useful to code up the convective cloud fraction changes to link directly
to the mass-flux convection scheme, and not to estimate them from the
values of :math:`Q4`, which can introduce errors.

.. _`sec:tbcs`:

Turbulence based convection scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will need to properly consider the links between PC2 and the
turbulence based convection scheme. In essence, we can use the diagnosed
cloud fraction and condensate values from the convection scheme to start
off the cloud again when convection has ceased. This has been tested to
some degree but will need proper analysis. The difficult decision comes
in choosing what to do with the condensate and cloud fraction that is
present *before* the convection starts, since we must ensure
conservation of moisture. This is not helped by the traditional view of
convective parametrization that ignores the existence of the condensate
phase in the atmosphere (i.e. it is only concerned with transport of
:math:`q` and :math:`\theta`, not of :math:`q_{cl}` and :math:`q_{cf}`)
despite the phase changes forming an integral part of the convection
scheme.

Detailed convective comparisons with CRM/LEM data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This work is already underway at the Met Office, in order to properly
evaluate the performance of the convective cloud parametrization in PC2
against high resolution research models.

Choice of PDF parameters
~~~~~~~~~~~~~~~~~~~~~~~~

Work by Dan Tang at Leeds University has highlighted an interesting and
undesirable property of the choice of :math:`m` and :math:`n` parameters
in the homogeneous forcing formulation. If a distribution is
homogeneously forced to :math:`C_l = 0`, then we do not necessarily get
:math:`\overline{q_{cl}}` tending to zero. This is because there is
enough influence from the :math:`\frac{{(1-C_l)}^2}{SD}` term in the
combination (`[eqn22] <#eqn22>`__) to stop the natural convergence of
the :math:`\frac{{C_l}^2}{\overline{q_{cl}}}` term to :math:`C_l =0` and
:math:`\overline{q_{cl}}=0`. Increasing the power of :math:`m` should
help. However, we note that the tests that have been done on the chosen
:math:`n` and :math:`m` values (0 and 0.5 respectively) do not show
particularly poor behaviour, and we do not pick up substantial evidence
of problems from this in the full model. This remains something to be
investigated.

.. _`sec:homog_improve`:

Homogeneous forcing section improvements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although the homogeneous forcing provides a convenient method to
calculate increments to :math:`C_l` and :math:`\overline{q_{cl}}`, it is
clearly not the best representation possible of the processes that use
it. For example, although the clear-sky radiative heating may perhaps
best be considered as a homogeneous process, the part of the radiative
heating influenced by clouds should, ideally, be applied to the cloudy
part of the gridbox and not the clear part. Vertical advection is likely
to be correlated with where there is already cloud, rather than being
uniform throughout the gridbox. There is no reason that a process that
uses homogeneous forcing as its condensation model should not be looked
at with a view to using something better. This is one of the strengths
of the PC2 framework and is an intention of the project.

Overlap of ice and liquid cloud changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have assumed within PC2 that ice and liquid cloud changes are
minimally overlapped with each other (within the same gridbox) in order
to maintain as much supercooled liquid water as possible. Although there
is good observational evidence to say that the two condensate phases
tend not to coexist together in a cloud, it may be possible to
characterise and apply this overlap in a more quantiative way.

Parameter tuning
~~~~~~~~~~~~~~~~

The sensitivity of some of the parameters in PC2 have not been properly
tested, mainly due to a lack of resources rather than a physical reason.
We have seen that the most effective method of tuning cloud is with the
erosion term, which has been increased to high values in order to remove
enough cloud and is probably as high as we reasonably wish to take it
given the length of the timestep.

- The phase change temperature (between liquid and ice) in the
  convective plume, TICE, is known to influence the strength of the
  convection through the latent heat differences. It also impacts on the
  amount of supercooled liquid water in the model. The quantitative
  impact of altering this could be explored. We note that CRM
  simulations of deep convection suggest that some supercooled liquid
  water exists within the plumes to :math:`-40 ^{\circ} C` and that a
  representation with partial liquid and partial ice phase would be more
  appropriate, based possibly on the current diagnosed convective cloud
  phase in the non-PC2 model. Although the theoretical work has been
  done to allow partial phases, we repeat the caution that care must be
  taken when doing the work and appropriate testing done to ensure that
  heat and moisture are properly conserved within the convection scheme.

- The growth of :math:`C_f` due to the fall-out of ice term in the
  microphysics is parametrized with a dependence on windshear. We have
  never linked this directly to the windshear, instead we have used
  estimated the windshear as a fixed value. There is no reason why the
  actual model windshear cannot be passed into the scheme in order to
  properly calculate this term.

- :math:`RH_{crit}` remains a tunable parameter. Although its impact is
  less than in a non-PC2 simulation, it is still significant in
  initiating cloud and in determining the evolution of the ice cloud.
  There is also an implicit overlap assumption regarding the ice cloud
  fractions, again this might be improved upon.

- :math:`n` and :math:`m` values in the homogeneous forcing have not
  been thoroughly investigated for a long time now, and may yield some
  sensitivities.

Cloud inhomogeneities
~~~~~~~~~~~~~~~~~~~~~

A cloud generator approach to cloud inhomogeneities is currently being
developed. However we note two particular issues that relate to PC2.

- The first is that in the diagnostic scheme, the two cloud fractions
  (convective and large-scale) allows, to some degree, a representation
  of cloud inhomogeneity. This is absent from PC2, although we note that
  the convective cloud fraction variable has not been removed from the
  radiative transfer code for PC2, it is merely set to zero, so it is
  easy to put back.

- The generation of inhomogeneities using a cloud generator requires
  some estimate of the variance (and possibly skewness) of the
  condensate in the gridbox. It is possible to back out the full
  moisture PDF at each grid point by homogeneous forcing (providing
  :math:`C_l` is not equal to 0 or 1), but this is very expensive and
  cannot be done on-line. Is there a quick *estimate* of the variance or
  skewness that it is possible to obtain from knowledge only of
  :math:`\overline{q}`, :math:`q_{sat}`, :math:`\overline{q_{cl}}` and
  :math:`C_l` etc.?

.. _`sec:timestepping`:

Time-stepping
~~~~~~~~~~~~~

A proper analysis of timestep sensitivities of PC2 (as opposed to
microphysics, convection etc) in the full UM or SCM has not been done
for a long time. In the early development stages much effort was placed
in developing good numerical techniques for each of the terms in PC2,
and to explore the way in which they coupled together. An example is the
homogeneous forcing timestep investigated by :raw-latex:`\cite{wg03}`.
We note that in shallow convection at 30 minutes timestep the erosion
term is trying to remove most of the cloud that the convective
detrainment places into the model. Since the erosion is limited by the
amount of cloud fraction and condensate present, what ends up happening
is that the ‘equilibrium’ that is achieved is actually one where the
cloud fraction and condensate at the end of the timestep are simply the
values that were detrained by the convection scheme (and hence depend on
the timestep). The CRM suggests a cycling time of around 15 minutes for
liquid water content and just less than half and hour for the cloud
fraction, so we would expect timestep dependency to occur from around a
timestep of 15 minutes upwards. We might just about get away with the 30
minute step of the climate model, but it is not a good situation to try
to model. This is demonstrating the difficulty of modelling shallow
convective cloud by a prognostic scheme, where the physical lifetime of
the clouds is of order the timestep - ideally we wouldn’t want to try to
model anything prognostically when the cycling time is less than the
timestep.

As discussed in section `4.3.4 <#sec:erosion_numerics>`__, the timestep
sensitivity of cloud amounts in shallow cumulus regimes can be addressed
by using a more accurate numerical method to solve the erosion term.
Several options are available under the UM namelist switch
**i_pc2_erosion_numerics**.

In the early development of PC2 we chose to incorporate the PC2 cloud
and condensation increments in the same location where the increments
were calculated (e.g. the microphysics cloud fraction increments get
added along with the microphysics :math:`\overline{T}` and
:math:`\overline{q}` increments). This choice was made in order not to
confuse the timestepping method in the UM, which has been carefully
developed over a number of years to achieve numerical accuracy. However,
we note that the rapidly varying nature (in space and time) of variables
such as :math:`\overline{q_{cl}}` and :math:`C_l` is very different from
the smooth fields of :math:`\overline{q_T}` and :math:`\overline{T}`,
for which the timestepping was developed, and it may not be appropriate
to implement these in the same locations. In particular, we might wish
to store the increments through the timestep and update values of
:math:`\overline{q_{cl}}` and :math:`C_l` etc. at the end of the
timestep, where many of the balances can be cancelled.

One issue is that we are calculating the increments due to condensation
associated with the adiabatic response to pressure changes after the
Helmholtz solver. Pragmatically, we need to do it here since we do not
know the arrival value of pressure until after the Helmholtz solver has
been used. However, in order to achieve balanced dynamical fields, it is
useful the Helmholtz solver to be called after all the latent heating
terms have been calculated (which not only includes the adiabatic
response to lifting but the cloud initiation term). We have shown that
PC2 can run with the two terms switched over, but this implies that we
are missing part of the pressure change following the parcel (the time
changing part rather than the spatially changing adiabatic part).
Although the adiabatic change is usually likely to dominate, it may be a
significant loss. Under the UM namelist switch **l_pc2_sl_advection**,
we can call the PC2 response twice, once before the Helmholtz solver and
once afterwards in order to pick up most of the latent heat change
before the solver, but not to have PC2 miss some of the pressure change.
The call for the advective part (before the Helmholtz solver) is
actually done before the call to atmos_physics2 as well, and so results
in more realistic, saturation-adjusted, profiles being passed to the
convection scheme.

We have placed the initiation at the end of the timestep, but it is
sensible to ask whether this could ideally be located elsewhere.

Initiation formulation
~~~~~~~~~~~~~~~~~~~~~~

Ideally this should be a relatively infrequent part of the model but
remains an essential part of the code. It is reasonable to ask whether
the initiation is optimal, particular in the diagnosis of when it is
applied. For example, we note that the initiation is currently
symmetrical, with initiation from :math:`C_l=1` occuring with the same
:math:`RH_{crit}` value as from :math:`C_l=0`. However, the
:raw-latex:`\cite{wf00}` observations hint that a higher
:math:`RH_{crit}` might be more appropriate for initiation from
:math:`C_l=1`.

70-levels performance
~~~~~~~~~~~~~~~~~~~~~

The performance of PC2:66 in the 70-levels model is not good as far as
shallow convective cloud is concerned (there is far too much of it in
the trade regions). It may be that PC2 is latching onto a convection
sensitivity that is present on going from L38 to L70 but had little
effect in a non-PC2 simulation. It may also be related to a reduction in
timestep from 30 minutes to 20 minutes. Investigations have not made
much progress in identifying the reasons for the differences, or
producing effective tunings to counter the problem.

High horizontal resolution performace
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PC2 has only been tested once at 4 km horizontal resolution. This
produced excessive shallow convective cloud (this may or may not be
related to the 70-levels problem above). Since this simulation the
erosion term has been increased dramatically, which may help. We note
that one of the main advantages of PC2, that of a prognostic link of
cloud to convection, is reduced at high resolution, as convection
becomes more explicit rather than diagnosed. We hence see a resolution
limit beyond which it is no longer appropriate to use PC2. Results look
acceptable at 12 km resolution, but we have not quantitatively explored
this limit.

Diagnostic evaluation
~~~~~~~~~~~~~~~~~~~~~

One of the principal areas for future cloud scheme development work
planned in the future is in the area of detailed evaluation against a
number of data sources, such as CloudSat, ground based radar, or case
study campaigns. The quantitative evaluation has been lacking to a
significant degree in the development of the scheme, as the focus has
been on tackling qualitatively poor results. Hence new sources of
evaluation work on PC2 would be very welcome.

Moisture distribution within the deposition/sublimation term
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The liquid cloud changes in PC2 (or in a non-PC2 run) are based upon a
moisture PDF, as are the deposition/sublimation changes. However, it is
not the same PDF. It has always been the case with the prognostic ice
microphysics term that its PDF, whether explicit or implicit, has not
been rigorously consistent with the PDF used in the calculation of
liquid water, because it was most easily developed that way and produced
reasonable results. It may be useful to investigate whether the two PDF
representations can be brought together in a rigourous way, both for the
PC2 scheme and the :raw-latex:`\cite{smith90}` scheme.

We have similarly noted potential inconsistencies in the parametrization
of cloud fraction changes between the evaporation of rain term and the
riming (or accretion) term. Again, it might be possible to bring
together these formulations into a single consistent framework.

Area cloud fraction representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current area cloud fraction representation is not used when
convection is taking place (signified by the *cumulus* logical). This
inevitably leads to a potential switching between two different values
of the cloud fields if the convective boundary layer (not whether the
convection is shallow or deep) switches on and off, which is
undesirable, although not as bad as switching cloud on and off
completely (as for the current convective cloud formulation).
Additionally, it is reasonable to argue that having an area cloud
fraction for cirrus cloud depend upon whether the boundary layer is well
mixed or has shallow convection occuring is not a reasonable link.

Work in Australia on a TWP-ICE single column model case study using PC2
suggests the area cloud fraction scheme over estimates the area cloud
coverage for tropical anvil clouds (which exist long after the
convection itself has ceased). This is perhaps not surprising since the
:raw-latex:`\cite{bhi05}` area cloud fraction scheme was evaluated
against mid-latitude cloud and it is known that tropical clouds have
greater vertical coherence. Tuning the parameters in
:math:`large_scale_cloud/ls_acf_brooks.F90` may be beneficial.

.. container:: float
   :name: fig:schematic

   .. container:: center

      |image1|

.. container:: float
   :name: fig:tstep_diag

   .. container:: center

      |image2|

.. container:: float
   :name: fig:tstep_prog

   .. container:: center

      |image3|

.. |image1| image:: pc2_process_explanation
.. |image2| image:: Timestepping_ctl66.epsi
.. |image3| image:: Timestepping_pc266.epsi
