.. _generation_methods:

Inflow generation methods
=========================

.. _using_generators:

Using the generators
--------------------

Each generation method provided by the library has an executable Python script
associated with it.
Running the appropriate script executes the generation procedure.
All inflow generation methods depend on a certain amount of parameters.
These parameters are communicated to the script via a configuration file, which
is passed as a command-line argument. ::

   nameOfTheScript --config=configurationFileName

The configuration file is a simple text file, with the following
layout ::

   # This is a comment, it can describe the parameter below
   parameterOne    valueOfParameter1

   # Another comment
   parameterTwo    valueOfParameter2

Which parameters should be present depends on the method that is being
used.
They are therefore described for each method individually below.

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
a backward-facing step.

In order to introduce the method, let us consider a zero-pressure flat plate
turbulent boundary layer (TBL).
At the inlet of the domain the desired characteristics of the boundary
layer, such as :math:`\text{Re}_\tau` and :math:`\text{Re}_{\delta_{99}}`, are
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
The rescaled velocity fields are then saved to the hard-drive and serve as
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

   & U^{\text{inner}}(y^+) = u_\tau f_1(y^+),\\
   & U_0 - U^{\text{outer}}(\eta) = u_\tau f_2(\eta).

Another assumption, that is fulfilled automatically in the setting proposed
by Lund et al :cite:`Lund1998`, but not within the framework of eddylicious, is
that the relationships above are valid for both the precursor simulation and
the main simulation.
Strictly speaking, this requires the precursor simulation to be a TBL itself.
However, a flow sufficiently similar to a TBL, like channel flow, can also be
used with success.

Let :math:`\gamma = u_{\tau, \text{infl}}/u_{\tau, \text{prec}}`.
Then, if the assumption above is fulfilled, the rescaling procedure for the mean
streamwise velocity is

.. math::

   &  U^\text{inner}_\text{infl}(y^+_\text{infl}) =
   \gamma U^\text{inner}_\text{prec}(y^+_\text{infl}),\\
   &  U^\text{outer}_\text{infl}(\eta_\text{infl}) =
   \gamma U^\text{outer}_\text{prec}(\eta_\text{infl}) + U_{0, \text{infl}} -
   \gamma U_{0, \text{prec}}.

The rescaling for the mean wall-normal velocity is defined simpler, and is
not as rigorously based on any physical assumption.

.. math::

   &  V^\text{inner}_\text{infl}(y^+_\text{infl}) =
   V^\text{inner}_\text{prec}(y^+_\text{infl}),\\
   &  V^\text{outer}_\text{infl}(\eta_\text{infl}) =
   V^\text{outer}_\text{prec}(\eta_\text{infl}).

The rescaling for the fluctuations is defined as

.. math::

   & (u'_i)^\text{inner}_\text{infl}(y^+_\text{infl}) =
   \gamma (u'_i)^\text{inner}(y^+_\text{infl}),\\
   & (u'_i)^\text{outer}_\text{infl}(\eta_\text{infl}) =
   \gamma (u'_i)^\text{outer}(\eta_\text{infl}).

The inner and outer components are blended together using a weighted average:

.. math::

   u_{i, \text{infl}} = u_{i, \text{infl}}^\text{inner}[1-W(\eta_\text{infl})] +
   u_{i, \text{infl}}^\text{outer}W(\eta_\text{infl}).

The weight function :math:`W(\eta)` is defined as

.. math::

   W(\eta) = \frac{1}{2} \left\{ 1+ \dfrac{\tanh \left( \frac{\alpha(\eta - b)}{(1-2b)\eta +b}\right)}{\tan(\alpha)} \right\},

where :math:`\alpha=4` and :math:`b=0.2`.



Usage and practical information
_______________________________

The `runLundRescaling` script should be used to generate the fields.
The script is parallelized using MPI, so it is possible to take advantage of
all the available cores present on the machine.

Depending on what data is available for the TBL desired at the inlet it may
be convinient to either use :math:`\delta_{99}` or :math:`\theta` as the outer
scale (that is the length used to normalize :math:`y` to obtain :math:`\eta`).
Eddylicous can work with both and will use the scale which is provided in the
config file, i.e. one of the two should be present:

   * ``delta99`` --- desired :math:`\delta_{99}` at the inlet of the main
     simulation.

   * ``theta`` --- desired momentum thickness  at the inlet of the main
     simulation.

Note that using :math:`\theta` requires to scale :math:`\eta` before it can be
plugged into function :math:`W(\eta)`.
The value of 8 is used, based on the fact that :math:`\theta` is around 8 times
less than :math:`\delta_{99}` for a wide range of Reynolds numbers.

As evident from the equations, defining the rescaling procedure,
the value of the friction velocity at the inlet, :math:`u_{\tau, \text{infl}}`,
is needed for the procedure.
To this end, two options are available to the user.
One is to simply provide the value of the friction velocity directly.
The other is to let eddylicious compute it using the skin friction coefficient,
:math:`c_f`, and an empirical estimate connecting it to either
:math:`\text{Re}_{\delta_{99}}` or :math:`\text{Re}_\theta`.

.. math::

   & c_f = 0.02 \text{Re}^{-1/6}_{\delta_{99}}, \\
   & c_f = 0.013435(\text{Re}_\theta - 373.83)^{-2/11}.

The friction velocity is then obtained as :math:`U_0 \sqrt{c_f/2}`.
The related parameter in the configuration file is

   * ``uTauInflow`` --- the friction velocity at the inlet of the main
     simulation. Either the value of the velocity or ``compute``, which
     tells eddylicious to use one of the correlations above.

Another important feature is that eddylicious will always use only half of the
datapoints in the wall-normal direction available from the precursor
simulation.
This is natural if the precursor is channel flow, but is in fact unnecessary
when it comes to rescaling from another TBL simulation.
Basically, this demands that the boundary layer used as a precusor does not
occupy more than half of the computational domain in the wall-normal direction.

It is possible to choose which half of the precursor plane to consider, the
bottom or the top.
The following parameter in the configuration file controls this choice.

   * ``half`` --- which half of the precursor plane to grab the data from.
     Either ``bottom`` or ``top``.

Note, that this means that a single channel flow precursor actually contains
two independent precursor datasets.

The rescaling formulas involve the velocity from the precursor simulation
evaluated for the values of :math:`y^+` and :math:`\eta` defined by the
TBL at the inflow of the main simulation.
These values are obtained using linear interpolation.
This means that the values of :math:`\text{Re}_\tau` and
:math:`\text{Re}_\theta` for the precursor simulation must be higher than that
at the inflow of the main simulation.
Applied to rescaling from a precursor TBL this means that one can only rescale
from "downstream".

In the current implementation, eddylicious will compute the highest value of
:math:`\eta` available for the precursor simulation.
Then it will pick the points in the main simulation for which :math:`\eta` is
lower than this computed value.
This ensures that interpolation is possible for the outer part of the profile.
These chosen points will be considered as containing the inflow TBL.
In all points above, the freestream velocity will be prescribed.
If the range of :math:`\eta` in the precusor is not sufficient to cover the
whole inflow TBL, a jump in the mean streamwise velocity will be observed.

Note, that no similar procedure is performed for :math:`y^+`.
Therefore, if the range of :math:`y^+` in the precursor does not cover that in
the inflow TBL, eddylicious will simply crash.

Besides for the parameters mentioned above, the configuration file should also
define the following parameters.

   * All parameters associated with the chosen input and output formats.
     Refer to the associated parts of the User guide for information.

   * ``yOrigin`` --- the wall-normal coordinate of wall which the boundary
     layer is attached to in the main simulation.
     This is used when evaluating non-dimensional coordinates like :math:`y^+`.
     Also this is used to determine the "orientation" of the TBL with respect
     to the wall-normal coordinate.

   * ``nuInflow`` --- the kinematic viscosity value in the main simulation.

   * ``nuPrecursor`` --- the kinematic viscosity value in the precursor
     simulation.

   * ``U0`` --- desired freestream velocity at the inlet of the main simulation.

   * ``dt``--- the time-step in the main simulation.

   * ``t0`` --- the start-time of the main simulation.

   * ``tEnd`` --- the end-time of the simulation.

   * ``tPrecision`` --- write precision for time values.
     Should be chosen according to ``dt``.

Example configuration files can be found in the tutorial
:ref:`tut_of_channel_lund`.
