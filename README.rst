================
Eddylicious
================

Eddylicious is a python package, designed for implementing methods for generating inflow boundary fields for Large Eddy Simulation (LES).
Currently the existing functionality is heavily coupled with OpenFOAM, so the users of that software will probably find the package to be most attractive.
However, by design the package is generic.

The package is divided into three main parts: readers, writers and generators.
Generators are functions that actually implement the generation of the boundary fields.
Writers are functions that implement output of the generated fields into a specific format, so that they can be used by a CFD solver.
Readers are functions that implement the input of either geometric data (such as point or face-center locations of the inflow plane) and/or flow field data, which can be needed for methods using a precursor-database.

Additionally, conifgurable scripts that perform the generation and output are provided.

Note, that this structure of the package means that independently of the generation method, the fields are stored in the form of a database.
This is somewhat counterintuitive in the case of a syntetic method, but the rationale is that this package can be used for quick development of methods (thanks to python) and serve as a hub uniting several methods under a shared framework. 

---------------
 Implemented Generation Methods
---------------

Currently, only one generation method is implemented in the package: the rescaling procedure presented in [lund]_.

---------------
 Implemented Output Formats
---------------
Two output formats are currently supported.

* HDF5. This is a file-format specifically designed to be a container for large datasets.
  More information can be found at the `offical homepage <https://www.hdfgroup.org/HDF5/>`_.
  HDF5 allows to store large arrays of data in the form of datasets.

  Three datasets are created by the writer.

  :points: 
    A 2d array of size N-by-3, where N is the number of points. 
    The three colums are the x, y and z coordinates of the points, respectively.
    The points represent the location of the datapoints for the generated velocity fields in space.
    In the finite volume setting, it is best that this corresponds to the face centers of the faces defining the inflow boundary.
    
  :time:
    A 1d array, the size corresponds to the number of time-values that the inflow boundary fields are generated for.
    Each element contains a time value.

  :velocity:
   A 3d array. The first dimension corresponds to time and has the same size as the time array. 
   The second dimension is of size N and the third is of size 3, therefore for each timestep there is a N-by-3 array that contains the values of the velocity field.
   The order corresponds to the points array.

  APIs for loading HDF5 datasets exist in most languages actively used in scientific computing.
  For using this output format with OpenFOAM see `this repository <https://bitbucket.org/lesituu/timevaryingmappedhdf5fixedvalue>`_.

* TimeVaryingMappedFixedValue (tvmf). This format corresponds to the TimeVaryingMappedFixedValue boundary condition in OpenFOAM, which allows to read in the values at the boundary from a file and perform interpolation in both space and time.
  Basically, this is a text file, but with some specific formatting required by OpenFOAM. 
  An example can be found `here <https://github.com/OpenFOAM/OpenFOAM-2.4.x/blob/master/tutorials/incompressible/simpleFoam/pitzDailyExptInlet/constant/boundaryData/inlet/0/U>`_.

  The advantage of this format is that it is native to OpenFOAM, but a huge disadvantage is that a separate file has to be created for each time-step.

---------------
 Implemented Input Formats
---------------   

Currently one reader is implemented -- for the foamFile format, the native output format of OpenFOAM.
This format is almost identical to tvmf, and is used by many utilities included in OpenFOAM.

---------------
 Executables
---------------  

runLundRescaling.py
=============
Command line options:

--config The configuration file, obligatory.

This script is currently the main workhorse of the package.
It performs the generation of the velocity field by rescaling a precursor database using the method of Lund et al [lund]_.
The script is configured by a config file which should be provided as an argument to the script.

Example config files can be found in the examples directory.
The obligatory parameters in the config file are also dependent on the choice of reader and writer.

Below is a discription of all the config paramters.

reader
    The reader that will be used to read in the data from the precursor simulation.

    *Valid choices:*

    - foamFile

inflowReader
    The reader that will be used to read in the points of the inflow surface.

    *Valid choices:*

    - foamFile

readPath
    A path provided for the reader.
    
    - foamFile reader: the location should be the path to the root of a OpenFOAM case.
      The reader will then look into the postProcessing/surfaces/*time*/*sampleSurfaceName*/vectorFields/U to get the velocity fields, where sampleSurfaceName is a config file parameter (see below).

inflowReadPath
    A path provided for the inflowReader.

    - foamFile inflowReader: location of the directory, where it will expect to find a directory with a name determined by the inflowPatchName parameter (see below), and there a file called faceCentres.
      The file is commonly generated by the sample utility in OpenFOAM, then the location is postProcessing/surfaces/0/.
 
inletPatchName
    Used by the tvmf writer and the foamFile inflowReader in order to read/write to the correct directory.
    This should be the name of the patch that you use as the inlet in your OpenFOAM case.

sampleSurfaceName
    Used by the foamFile reader to read in the values of velocity from the precursor.
    Should be the name that you used for the sampling surface in your sampleDict.

writer 
    The writer that will be used to output the generated fields.

    *Valid choices:*

    - hdf5

    - tvmf

writePath
    Path provided for the writer.

    - hdf5 writer. Can be an arbitrary path.

    - tvmf writer. Expects the root of an OpenFOAM case.
      The constatnt/boundaryData/*inflowPatchName* directory should exist.

hdf5FileName
    Used by the hdf5 writer only.
    The name of the created hdf5 file. 
    The file will be created in the writePath. 

xOrigin
    The x-coordinate of the origin of the inflow plane.
    This value will be added to the x-coordinates of the points array.        
yOrigin
    The y-value of the origin of the inflow plane.
    This is used when evaluating non-dimensional coordinates like y+.        

nuInflow
    The viscosity value in the simulation

nuPrecursor
    The viscosity value in the precursor

uTauInflow
    The friction velocity in the simulation, at the inlet. One can either provide the value of write "compute" withou the quotes.
    In the latter case the friction velocity will be computed based on delta99 and Ue.

uTauPrecursor
    The friction velocity in the precursor.

delta99
    Desired boundary layer thickness.

Ue
    Desired free-stream velocity.

dt              0.0005
    The time-step in the simulation

t0
    The start-time in the simulation.

tEnd
    The end-time in the simulation.

tPrecision
    Write precision for time values.
    Should be chosen according to dt. 


[lund] T. S. Lund, X. Wu, and K. D. Squires. On the Generation of Turbulent Infow Conditions for Boundary Layer Simulations. Journal of Computational Physics, 140:233-258, 1998.