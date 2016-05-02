Supported input formats
=======================

Some of the methods for inflow generation require input.
Typically these are precursor-based methods, such as the one based on the recaling procedure presented in :cite:`Lund1998`.

FoamFile format
---------------

This file format is associated with the CFD solver OpenFOAM.
Data from a user-defined sampling surface can be saved in this format.
OpenFoam creates a catalog each time the data is output, named with the value of the simulation time.
At the root level of the catalog the data related to the mesh is stored.
That includes the coordinates of the points, a list of faces each defined as list of points, and the coordiantes of the face centres.
Since the data resides in the faces centres, they are used in eddylicious to represent the geometry of the surface.
A folder `vectorField` is created to store field with vector-valued data.
In particular, the field `U` which represents the velocity field, will be located there.

Format of the face centres file
_______________________________

Let :math:`N_p` be the total number of face centres.
The face centres file has the following simple format.

:math:`N_p`
(
(:math:`x_0`, :math:`y_0`, :math:`z_0`)
(:math:`x_1`, :math:`y_1`, :math:`z_1`)
(:math:`x_2`, :math:`y_0`, :math:`z_2`)
.
.
.
(:math:`x_{p-1}`, :math:`y_{p-1}`, :math:`z_{p-1}`)

)





