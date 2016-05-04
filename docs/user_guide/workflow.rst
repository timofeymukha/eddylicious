===================
Workflow guidelines
===================

Here the suggested workflows for using eddylicious in conjuction with various
solvers are presented.
Basically, whatever solver is used, the following steps have to be performed.

    * Specifying the geometry of the inlet for eddylicious.

    * Choosing an output format that is compatible with the used solver.

    * Setting up the solver to read in boundary data from the hard drive.

    * In case a precursor database is being built, the velocity field should
      be saved or later converted to a format that eddylicious can read.


Using eddylicious with OpenFOAM
-------------------------------

.. important::

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software and owner of the OPENFOAM®  and
    OpenCFD®  trade marks.

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
The path to the ``faceCenters`` file should also be also provided.
This is done by adding the following lines to the configuration file for the
inflow generation script. ::

    inflowGeometryReader    foamFile
    inflowGeometryPath      "/path/to/faceCenters/file"

Reading the inflow fields from OpenFOAM
_______________________________________



Creating a precursor database
_____________________________

