.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _coordinates_and_vectors_rotation_and_stretching:

General Rotated and Stretched Meshes
====================================

Mesh rotation and stretching are defined as separate operations.
Rotation is specified through the location of the new North Pole,
:math:`(\lambda_N, \varphi_N)`, while stretching is described by a Schmidt
transform with stretching factor :math:`s`.

Rotation
--------

Mesh rotation is performed in two stages:

1. rotation by :math:`\vartheta` about the original North Pole;
2. rotation of the North Pole itself.

The first stage is represented by

.. math::

    \mathbf{R}^\dagger =
    \begin{bmatrix}
    \cos\vartheta & -\sin\vartheta & 0 \\
    \sin\vartheta & \cos\vartheta & 0 \\
    0 & 0 & 1
    \end{bmatrix}.

The second stage is a Rodrigues rotation. Writing :math:`\omega = \pi/2 -
\varphi_N` and
:math:`\widehat{\mathbf{\varOmega}} = (-\sin\lambda_N, \cos\lambda_N, 0)^T`,
the rotated coordinates are

.. math::

    \widehat{\mathbf{X}} = \mathbf{X}\cos\omega +
    (\widehat{\mathbf{\varOmega}}\times\mathbf{X})\sin\omega +
    \widehat{\mathbf{\varOmega}}(\widehat{\mathbf{\varOmega}}\cdot\mathbf{X})(1-\cos\omega).

Using :math:`\cos\omega = \sin\varphi_N` and :math:`\sin\omega =
\cos\varphi_N`, the source notes expand the Cartesian components explicitly.
The corresponding matrix form is

.. math::

    \widehat{\mathbf{R}} =
    \begin{bmatrix}
    \sin\varphi_N + (1-\sin\varphi_N)\sin^2\lambda_N &
    -\sin\lambda_N\cos\lambda_N(1-\sin\varphi_N) &
    \cos\lambda_N\cos\varphi_N \\
    -\sin\lambda_N\cos\lambda_N(1-\sin\varphi_N) &
    \sin\varphi_N + (1-\sin\varphi_N)\cos^2\lambda_N &
    \sin\lambda_N\cos\varphi_N \\
    -\cos\lambda_N\cos\varphi_N &
    -\sin\lambda_N\cos\varphi_N &
    \sin\varphi_N
    \end{bmatrix}.

In LFRic, the azimuthal rotation about the original North Pole is

.. math::

    \vartheta = \left\{
    \begin{matrix}
      \lambda_N + \pi & |\lambda_N| \geq \lambda_s \\
      0 & |\lambda_N| < \lambda_s
    \end{matrix}\right.,

where :math:`\lambda_s` is a small threshold, currently taken as
``10^{-4}`` degrees. This gives

.. math::

    \mathbf{R}^\dagger =
    \begin{bmatrix}
    -\cos\lambda_N & \sin\lambda_N & 0 \\
    -\sin\lambda_N & -\cos\lambda_N & 0 \\
    0 & 0 & 1
    \end{bmatrix}

for :math:`|\lambda_N| \geq \lambda_s`, and the identity matrix otherwise.

The overall rotation is then

.. math::

    \mathbf{R} = \widehat{\mathbf{R}}\mathbf{R}^\dagger.

Note that there are two current limitations to this implementation of rotation:

1. Rotation around the pole is currently by an angle :math:`\vartheta = \lambda_N + \pi`
   unless :math:`|\lambda_N| < \lambda_s`, which leads to unsmooth behaviour as
   the mesh rotation is discontinuous at :math:`|\lambda_N| = \lambda_s`.
2. Specification of the North Pole alone does not fully determine the rotation,
   so many possible rotations are not currently supported.

Stretching
----------

Stretching is applied towards the North or South Pole using a Schmidt
transform. The stretching factor :math:`s` satisfies :math:`s=1` for no
stretching, :math:`s<1` for stretching towards the South Pole and :math:`s>1`
for stretching towards the North Pole.

The transformed latitude is

.. math::

    S_\varphi = \sin^{-1}\left[
    \frac{1-s^2 + (1+s^2)\sin\varphi}{1+s^2 + (1-s^2)\sin\varphi}
    \right],

with inverse

.. math::

    \varphi = \sin^{-1}\left[
    \frac{(1+s^2)\sin S_\varphi - (1-s^2)}{(1+s^2) - (1-s^2)\sin S_\varphi}
    \right].

The transformation to stretched Cartesian coordinates is

.. math::

    S_X = X\sqrt{\frac{r^2 - S_Z^2}{r^2 - Z^2}}, \qquad
    S_Y = Y\sqrt{\frac{r^2 - S_Z^2}{r^2 - Z^2}}, \qquad
    S_Z = r\left(\frac{1-s^2 + (1+s^2)Z/r}{1+s^2 + (1-s^2)Z/r}\right),

with inverse

.. math::

    X = S_X\sqrt{\frac{r^2 - Z^2}{r^2 - S_Z^2}}, \qquad
    Y = S_Y\sqrt{\frac{r^2 - Z^2}{r^2 - S_Z^2}}, \qquad
    Z = r\left(\frac{(1+s^2)S_Z/r - (1-s^2)}{(1+s^2) - (1-s^2)S_Z/r}\right).

The Jacobian for the stretching transformation is built from three pieces:

* the Jacobian from geocentric Cartesian to spherical polar coordinates,
* the Jacobian for stretching in spherical polar coordinates,
* the Jacobian from stretched spherical polar to stretched geocentric Cartesian
  coordinates.

The stretch Jacobian in spherical polar coordinates is

.. math::

    \mathbf{J}^{(\lambda S_\varphi r,\lambda\varphi r)} =
    \begin{bmatrix}
    1 & 0 & 0 \\
    0 & \psi & 0 \\
    0 & 0 & 1
    \end{bmatrix}, \qquad
    \psi = \frac{\partial S_\varphi}{\partial \varphi} =
    \frac{2s}{1+s^2 + (1-s^2)\sin\varphi}.

The full stretching Jacobian is therefore

.. math::

    \mathbf{J}_S =
    \mathbf{J}^{(S_X S_Y S_Z,\lambda S_\varphi r)}
    \mathbf{J}^{(\lambda S_\varphi r,\lambda\varphi r)}
    \mathbf{J}^{(\lambda\varphi r,XYZ)}.

General rotated and stretched meshes
------------------------------------

For general rotated or stretched meshes, the native primal and dual basis
vectors are obtained through the chain rule. These can be expressed as

.. math::

    \left[\dot{\mathbf{e}}_{(\dot{f}_i)}\right]_j =
    \frac{\partial \chi^{(XYZ)}_j}{\partial \dot{\chi}^{(f)}_i}
    = \sum_{kl} R_{jk}\frac{\partial S_k}{\partial F_l}
      \frac{\partial F_l}{\partial \dot{\chi}^{(f)}_i}
    = \sum_{kl} R_{jk}\frac{\partial S_k}{\partial F_l}
      \left[\mathbf{e}_{(f_i)}\right]_l,

and

.. math::

    \left[\dot{\mathbf{e}}^{(\dot{f}_i)}\right]_j =
    \frac{\partial \dot{\chi}^{(f)}_i}{\partial \chi^{(XYZ)}_j}
    = \sum_{kl} \left[\mathbf{e}^{(f_i)}\right]_k
      \frac{\partial S^{-1}_k}{\partial \chi^{(XYZ)}_l} R^{-1}_{lj}.

The Jacobian from general native coordinates to geocentric Cartesian
coordinates is

.. math::

    \mathbf{J}^{(XYZ,\dot{f})} = \mathbf{R}\,\mathbf{J}_S\,\mathbf{J}^{(\dot{X}\dot{Y}\dot{Z},\dot{f})}.
