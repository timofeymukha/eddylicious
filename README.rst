================
Eddylicious
================

Eddylicious is a python package, designed for implementing methods for generating inflow boundary fields for Large Eddy Simulation (LES).

The package is divided into three main parts: readers, writers and generators.
Generators are functions that actually implement the generation of the boundary fields.
Writers are functions that implement output of the generated fields into a specific format.
Readers are functions that implement the input of either geometric data (such as point or face-center locations of the inflow plane) and/or flow field data, which can be needed for methods using a precursor-database.