==============================
 Single Ray Traversal Program
==============================

This program is designed to fire a single ray through the OBB tree of a volume in the model (indicated in the command line by its GLOBAL_ID).

Upon successful completion, it will produce a .vtk file containing the ray that was specified (also in the command line) and all of the oriented boxes encountered upon traversal of the tree for visualization.


=======
 Notes
=======

This program relies on a version of MOAB in which the OrientedBox class header file is installed and available for use.

See `here <https://bitbucket.org/pshriwise/moab/branch/expose_ob>` for more information on that.
