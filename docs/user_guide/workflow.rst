.. _workflow:

===================
Workflow guidelines
===================

Here the suggested workflows for using eddylicious in conduction with various
solvers are presented.
Basically, whatever solver is used, the following steps have to be performed.

   * Specifying the geometry of the inlet for eddylicious.

   * Choosing an output format that is compatible with the used solver.

   * Generating the inflow fields by running the python script associated
     with the chosen inflow generation method.

   * Setting up the solver to read in boundary data from the hard drive.

.. _workflow_openfoam:

Using eddylicious with OpenFOAM
-------------------------------

.. important::

   This offering is not approved or endorsed by OpenCFD Limited, producer
   and distributor of the OpenFOAM software and owner of the OPENFOAM®  and
   OpenCFD®  trade marks.

.. _inlet_geometry_openfoam:

Specifying the geometry of the inlet
____________________________________

Specifying the geometry boils down to producing the list of face centres
located at the inlet boundary.
The coordinates of the face centres can be used using the ``sample`` utility
shipped with OpenFOAM.

In the ``system/sampleDict`` file, create a sampling surface with the type
``patch``, and specify the inlet patch as the basis for the surface.
It is better to turn off triangulation to preserve the original geometry.
Choose ``foamFile`` as the write format for surfaces.

Even though we are only interested in the geometry, a field for sampling has
to be chosen.
Any field can be chosen, besides for the velocity field ``U``.
This is because  the sample utility will attempt to read in the field,
and, since we didn't generate it yet, the field-values simply don't exist yet.

A ``sampleDict`` for a case with the inlet patch named ``inlet`` might look
something like this. ::

   *--------------------------------*- C++ -*----------------------------------*\
   | =========                 |                                                 |
   | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
   |  \\    /   O peration     | Version:  2.3.1                                 |
   |   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
   |    \\/     M anipulation  |                                                 |
   \*---------------------------------------------------------------------------*/
   FoamFile
   {
       version     2.0;
       format      ascii;
       class       dictionary;
       object      sampleDict;
   }
   surfaceFormat foamFile;

   fields
   (
      p
   );

   surfaces
   (
       inlet
       {
           type patch;
           patches (inlet);
           interpolate false;
           triangulate false;
       }

   );

The produced file can be read using the ``foamFile`` inflow geometry reader.
The path to the ``faceCenters`` file should also be provided.
This is done by adding the following lines to the configuration file for the
inflow generation script. ::

   inflowGeometryReader    foamFile
   inflowGeometryPath      "/path/to/faceCenters/file"

.. _reading_fileds_openfoam:

Reading the inflow fields from OpenFOAM
_______________________________________

OpenFOAM has a special boundary condition that allows reading boundary data
from a file, it is called ``timeVaryingMappedFixedValue``.
A tutorial, which takes advantage of this boundary condition, is shipped
with OpenFOAM.
It can be found under
``tutorials/incompressible/simpleFoam/pitzDailyExptInlet/``.

The boundary condition is quite flexible.
If needed, interpolation in space will be used to obtain the vales at the face
centres from the values at the provided points.
Linear interpolation in time is also supported.

Let ``inlet`` be the name of the patch for which the inflow fields are
generated.
Then the following entry should be found in the ``U`` file. ::

   inlet
   {
       type            timeVaryingMappedFixedValue;
       offset          (0 0 0);
       setAverage      off;
       perturb         0;
   }

Setting ``perturb`` to 0 is important, since this option perturbs the location
of the points.

In order to generate the inflow fields the :ref:`of_native_format` should be
used for writing the velocity fields to the hard drive.

Note that, for a large time-span, the amount of files written to disk become
extremely large.
To rectify this issue, a modified version of ``timeVaryingMappedFixedValue``
that reads all the data from a single HDF5 file is available.
For more information regarding the structure of the file see
:ref:`hdf5_file_format`.

The modified boundary condition is called ``timeVaryingMappedHDF5FixedValue``
and can be downloaded at
https://bitbucket.org/lesituu/timevaryingmappedhdf5fixedvalue

If this boundary condition is used the entry in the ``U`` file should look
as follows. ::

   inlet
   {
       type            timeVaryingMappedHDF5FixedValue;
       setAverage      false;
       perturb         0;
       offset          (0 0 0);
       hdf5FileName    nameofthehdf5file.hdf5;
       hdf5PointsDatasetName    points;
       hdf5SampleTimesDatasetName    time;
       hdf5FieldValuesDatasetName    velocity;
   }

