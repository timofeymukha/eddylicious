Available inflow generation methods
===================================

Common notation
---------------

Here we define the notation for several important physical quantities.

    * :math:`U_e` --- the edge, or freestream, velocity.

    * Measures of the thickness of a boundary layer.
        * :math:`\delta_{99}` --- the location where the mean velocity is equal
          to :math:`0.99U_e`.
        * :math:`\delta_*` --- the displacement thickness.
        * :math:`\theta` --- the momentum thickness.
    * The velocity field.
        * :math:`u, \: v, \: w` -- the streamwise, wall-normal and spanwise
          components of the velocity field.
        * :math:`U, \: V, \: W` -- the mean streamwise, wall-normal, and
          spanwise components of the velocity field.
        * :math:`u', \: v', \: w'` -- the streamwise, wall-normal and spanwise
          components of the fluctuations of the velocity field.


Lund's rescaling
----------------

Theory
______

This method has been presented in a paper by Lund, Wu and Squires
:cite:`Lund1998`.
The method is suitable for generating an accurate inflow field for a simulation
which involves a turbulent boundary layer approaching from upstream.
For example, it can be used in a simulation of flow around a smooth ramp or
a backwards-facing step.

In order to introduce the method, let us consider a zero-pressure flat plate
turbulent boundary layer (TBL).
At the inlet of the domain the desired  characteristics of the boundary
layer, such as :math:`\text{Re}_\tau` and :math:`\text{Re}_{\delta^{99}}`, are
known.

The central idea of the method is to obtain the values at the inlet by
rescaling the velocity field from a plane located some place downstream, see
:numref:`fig-lund-rescaling`.
The rescaling procedure involved insures that the target characteristics of
the boundary layer at the inlet are met.
The plane located downstream is referred to as the recycling plane.

.. _fig-lund-rescaling:

.. figure:: /figures/lund_rescaling.svg
    :align: center

    Schematic of the rescaling from a plane located downstream proposed by
    Lund et al :cite:`Lund1998`.

The rescaling procedure is based on self-similarity.
For TBL it is customary to define several regions



