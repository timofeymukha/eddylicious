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

.. figure:: /figures/lund_rescaling.*
   :align: center

   Schematic of the rescaling from a plane located downstream proposed by
   Lund et al :cite:`Lund1998`.

Since eddylicious does not interact directly with an ongoing simulation,
the method has to be reformulated.
Instead of using a recycling plane located downstream of the inlet in the main
simulation, the proposed approach is to have a separate, precursor simulation
dedicated to generating a database of two-dimensional velocity distributions.
This database then serves as input for eddylicious, which applies the rescaling
procedure as defined in :cite:`Lund1998`.
The rescaled velocity fields are then saved to the hard-drive as serve as
inflow fields in the main simulation.
The whole processes is schematicaly illustrated in
:numref:`fig-lund-rescaling-eddylicious`

.. _fig-lund-rescaling-eddylicious:

.. figure:: /figures/lund_rescaling_eddylicious.*
   :align: center

   Schematic showing how the rescaling proposed in :cite:`Lund1998` is
   implemented in eddylicious.

From now on the subscript prec will be used to refer to the values obtained
in the precursor simulation.
The subscript infl will be used to refer to values at the inlet of the main
simulation

The rescaling procedure is based on the assumption of existence of similarity
solutions for the mean stream-wise velocity profile in the inner and outer
layers of a TBL.
The following relationships stem from this assumption.

.. math::

   \begin{align}
        & U^{\text{inner}}(y^+) = u_\tau f_1(y^+),\\
        & U_e - U^{\text{outer}}(\eta) = u_\tau f_2(\eta).
   \end{align}

Another assumption, that is fulfilled automatically in the setting proposed
by Lund et al, but not within the framework of eddylicious, is that the
relationships above are valid for both the precursor simulation and the main
simulation.
Strictly speaking, this requires the precursor simulation to be a TBL itself.
However, a flow sufficiently similar to a TBL, like channel flow, can also be
used with success.

Let :math:`\gamma = u_{\tau, \text{infl}}/u_{\tau, \text{prec}}`.
Then, if the assumption above is fulfilled, the rescaling procedure for the mean
streamwise velocity is

.. math::

   \begin{align}
      & U^\text{inner}_\text{infl}(y^+) =
      \gamma U^\text{inner}_\text{prec}(y^+),\\
      & U^\text{outer}_\text{infl}(\eta) =
      \gamma U^\text{outer}_\text{prec}(\eta) + U_{e, \text{infl}} -
      \gamma U_{e, \text{prec}}.
   \end{align}

The rescaling for the fluctuations is defined as

.. math::

   \begin{align}
      & (u'_i)^\text{inner}_\text{infl}(y^+) =
      \gamma (u'_i)^\text{inner}(y^+),\\
      & (u'_i)^\text{outer}\text{infl}(\eta) =
      \gamma (u'_i)^\text{outer}(\eta).
   \end{align}

The inner and outer components are blended together using a weighted average:

.. math::
   \begin{align}
      u_{i, \text{infl}} = u_{i, \text{infl}}^\text{inner}[1-W(\eta)] +
      u_{i, \text{infl}}^\text{outer}W(\eta).
   \end{align}

The weight function :math:`W(\eta)` is defined as

.. math::
   \begin{equation}
      W(\eta) = \frac{1}{2} \left\{ 1+ \dfrac{\tanh \left( \frac{\alpha(\eta - b)}{(1-2b)\eta +b}\right)}{\tan(\alpha)} \right\},
   \end{equation}

where :math:`\alpha =4` and :math:`b=0.2`.

Usage
_____

The `runLundRescaling.py` script should be used to generate the fields.
The script is paralleled using MPI, so it is possible to take advantage of all
the available cores present on the machine.

The configuration file should define the following parameters.

   * ``reader`` --- the file format that is used to store the precursor
     database, see :ref:`input_formats` for available options.

   * ``inflowGeometryReader`` --- the file format that is used to store the
     points defining the inflow surface, see :ref:`input_formats` for available
     options.

   * Other parameters that are needed for the ``reader`` and
     ``inflowGeometryReader``, this depends on the formats you use.

   * ``writer`` --- the file format that will be used to store the generated
     fields, see :ref:`output_formats` for available options.

   * Other parameters that are needed for the ``writer``, this depends on the
     format you use.

   * ``xOrigin`` --- the x-coordinate of the origin of the inflow plane.
     This value will be used as the x-coordinate of the inlet.

   * ``yOrigin`` --- the wall-normal coordinate of the origin of the inflow
     plane.
     This is used when evaluating non-dimensional coordinates like :math:`y^+`.

   * ``nuInflow`` --- the viscosity value in the main simulation.

   * ``nuPrecursor`` --- the viscosity value in the precursor simulation.

   * ``uTauInflow`` --- the friction velocity at he inlet of the main
     simulation, at the inlet.

   * ``uTauPrecursor`` --- the friction velocity at the sampling plane in the
     precursor simulation.

   * ``delta99`` --- desired :math:`\delta_{99}` at the inlet of the main
     simulation.

   * ``Ue`` --- desired freestream velocity at the inlet of the main simulation.

   * ``dt``--- the time-step in the main simulation.

   * ``t0`` --- the start-time of the main simulation.

   * ``tEnd`` --- the end-time of the simulation.

   * ``tPrecision`` --- write precision for time values.
     Should be chosen according to dt.

Example config files can be found in the tutorial :ref:`tut_of_channel_lund`.
