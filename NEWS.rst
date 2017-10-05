NEWS
====

v. 0.0.1
--------

New functionality
_________________

 * New Interpolation generator for interpolating data between to grids.

 * New foamFile reader for unstructured points.


Bug fixes
_________

 * Corrected output of cf, u_tau and y+_1 for the Lund rescaling generator.

API changes
___________

 * Reader functions re-named following the pattern
   read_points_from_foamfile -> read_structured_points_foamfile.
