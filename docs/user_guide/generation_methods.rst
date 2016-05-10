.. _generation_methods:

Available inflow generation methods
===================================

.. _lund_rescaling:

Lund's rescaling
----------------

Theory
______

This method was presented in a paper by Lund, Wu and Squires
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



