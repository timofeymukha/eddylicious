Channel flow in OpenFOAM using Lund's rescaling
===============================================

.. important::

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software and owner of the OPENFOAM®  and
    OpenCFD®  trade marks.


In this tutorial, channel flow will be simulated.
First, a precursor simulation will be setup, which will compute channel flow
using periodic boundary conditions in both the streamwise and spanwise
directions.

The velocity fields created by the precursor will be used as input for
the rescaling procedure developed by Lund et al :cite:`Lund1998`.
It will be used to generate the boundary conditions for the main simulation,
which is also channel flow but with inlet and outlet boundaries in the
streamwise direction.

After completing this tutorial you will be able to do the following.

    * Set-up a precusor simulation in OpenFOAM.

