
HV Characterization Study
=========================

This repository contains programs for the study and understanding of bounding
volume hierarchy construction (BVH) around high-valence mesh features. The Mesh
Oriented DataBase [MOAB](https://bitbucket.org/fathomteam/moab/) is used to do
so.

This study centers around the construction of a high-valence test model. This
model begins as the surface faceting of a cube (2 triangles per face). The
faceting of one face is then deleted and replaced with a manually generated
high-valence region which can be altered in size based on the fractional area of
the cube face. The valence of the region can also be controlled.

Several programs exist in this repository to examine MOAB's bounding volume
hierarchy construction around the high-valence test model and in some cases
other models as well.

hv_cube
-------

Generates a high valence test model with area fraction, $a_f$, and valence,
$v_n$. It then saves the file as "hvcube.h5m" for further examination.

`./hv_cube <A_f> <v_n>`

sweep
-----

Sweeps the area fraction and/or valence of the high-valence model from 0 to the
specified value. It will then fire 600k randomly oriented rays using the BVH
into the high valence surface of the cube and time the result. It will then
write this information for each model generated to a file called "data.dat".

`./sweep <A_f> <v_n> `

The number of intervals used to approach the specified area fraction and valence
are input using the options `--af_ints` and `--v_ints`, respectively.

A specific parameter of MOAB's BVH hierarchy constructor can also be swept if
desired. The wosrt splitting ratio (wsr) can be set and iterated for all area
fractions and valencies by setting the `--min_wsr`, `--max_wsr`, and
`--wsr_ints` options.


Other Handy Programs
====================

write_obbs
----------

This program will load a MOAB model into a DAGMC instance, generate the BVH, and
then write the BVH into a .vtk database for visualization puroposes. Each file
in the vtk databse represents all of the boxes for a given level in the BVH.

`./write_obbs <filename>`

If desired, the triangles of any leaf nodes can be written along with their .vtk
file for more information using the `--with-tris` option. A nice demonstration
of visualzation using the resulting .vtk database can be seen [here](https://www.youtube.com/watch?v=w16oiYxFJJc).

count_leaves
------------

This program counts the number of leaf nodes in the BVH for a given model. It
also provides some information about the average number of leaves in a model's
BVH.

`./count_leaves <filename>`

This is a nice way to ensure that a well-formed BVH is being created for the file.


