.. _input_formats:

Supported input formats
=======================

Some of the methods for inflow generation require input.
Typically these are precursor-based methods, such as the one based on the
rescaling procedure presented in :cite:`Lund1998`.

Also, the geometry of the inlet is defined via reading the appropriate data
from a file.

.. warning::
   At this point eddylicious only supports rectangular inlets meshed with
   a rectilinear grid.
   The inlet plane is assumed to be perpendicular to the flow direction,
   which, in turn, is assummed to be aligned with the :math:`x` axis.

.. _foamfile_format:

The foamFile format
-------------------

This file format is associated with the CFD solver OpenFOAM.
Data from a user-defined sampling surface can be saved in this format.
OpenFOAM creates a catalog each time the data is output, named with the value
of the simulation time.
Inside this catalog yet another folder is created, named identically to the
name of the sample-surface as defined by the user.
At the root level of the catalog the data related to the mesh is stored.
That includes the coordinates of the points, a list of faces each defined as
list of points, and the coordinates of the face centres.
Since the data resides in the face centres, they are used in eddylicious to
represent the geometry of the surface.
A folder ``vectorField`` is created to store fields with vector-valued data.
In particular, the field ``U`` which represents the velocity field, will be
located there.

Let :math:`N_p` be the total number of face centres.
The file containing the face centres, called ``faceCentres``, has the following
simple format.

.. admonition:: The format of the face centres file

   .. math::

       & N_p\\
       & (\\
       & (x_0, \: y_0, \: z_0)\\
       & (x_1, \: y_1, \: z_1)\\
       & (x_2, \: y_2, \: z_2)\\
       & \vdots\\
       & (x_{p-1}, \:  y_{p-1}, \: z_{p-1})\\
       & )\\

The format of the file containing a sampled vector field is identical,
however instead of the coordinates each row contains the three components of
the vector.

The order in which the data is written corresponds to the order in which the
face centres are written to ``faceCentres``.

In order for eddylicious to read in the geometry of the inlet stored as a list
of face centres in the foamFile format the following should be added to the
configuration file. ::

   inflowGeometryReader    foamFile
   inflowGeometryPath      /path/to/the/faceCentres/file
   xOrigin                 the coordinate of the inlet on the axis parallel to
                           the main streamwise direction

In order for eddylicious to read in previously sampled velocity fields stored
in the foamFile format the following should be added to the
configuration file. ::

   reader                  foamFile
   readPath                /path/to/OpenFOAM/case
   sampleSurfaceName       name of sample surface as defined in controlDict

Eddylicious will search for the catalogs containing the data for different
time steps in
``readPath/postProcessing/sampledSurface/*time_value*/sampleSurfaceName``.

.. important::

   This offering is not approved or endorsed by OpenCFD Limited, producer
   and distributor of the OpenFOAM software and owner of the OPENFOAM®  and
   OpenCFD®  trade marks.

.. _input_hdf5_file_format:

The HDF5 format
---------------

HDF5 is a file format specifically developed for storing large scientific
datasets.
More details regarding HDF5 can be found `here <https://www.hdfgroup.org/>`_.

In HDF5, data can be sorted into groups.
The data itself is stored in the form of datasets.
One can think of datasets as of multidimensional arrays.

In eddylicious the file is expected to contain two groups: ``points`` and
``velocity``.

Let :math:`N_y` and :math:`N_z` be the number of points in inlet plane in the
wall-normal and spanwise direction respectively.

The ``points`` group contains the following two-dimensional datasets:

   * ``pointsY``, :math:`N_y \times N_z` --- dataset with the wall-normal
     coordinates of the points.

   * ``pointsZ``, :math:`N_y \times N_z` --- dataset with the spanwise
     coordinates of the points.

This structure implies, that all the columns of ``pointsY`` are identical, as
well as all the rows of ``pointsZ``.

Let :math:`N_t` be the amount of time steps for which velocity data is
available.

The ``velocity`` group contains the following three-dimensional datasets:

   * ``uX``, :math:`N_t \times N_y \times N_z` --- dataset with the values of
     the streamwise component of the velocity field.

   * ``uY``, :math:`N_t \times N_y \times N_z` --- dataset with the values of
     the wall-normal component of the velocity field.

   * ``uZ``, :math:`N_t \times N_y \times N_z` --- dataset with the values of
     the spanwise component of the velocity field.

The values located at position ``[k, i, j]`` in these arrays correspond to
the point with coordinates (``pointsY[i, j]``, ``pointsZ[i,j]``).

Also, the following one-dimensional arrays are stored in the ``velocity``
group:

   * ``uMeanX``, :math:`N_y` --- the values of the mean streamwise velocity.

   * ``uMeanY``, :math:`N_y` --- the values of the mean wall-normal velocity.

   * ``times``, :math:`N_t` --- the time values associated with the velocity
     fields.

The way the data is stored in the HDF5 file coincides with how it is
represented internally.
This implies that the overhead from reading the data is minimal.
HDF5 supports parallel processing of the data via MPI.
Different processes can therefore read in the required data simultaneously.

Therefore, this file format can be considered optimal.
Since solvers will not typically support output in this particular format,
utilities for converting a precursor database saved in a different format
into the HDF5 format are part of eddylicious.

In order for eddylicious to read in previously generated velocity fields stored
as an HDF5 file, the following should be added to the configuration file. ::

   reader                  hdf5
   readPath                /path/to/hdf5/file
