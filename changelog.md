 Update README.md

Version 1.1.0

a) Fixed an issue when attempting to perform anisotropic volume changes which deforms the shape of the simulation box;

b) Added box distortion parameters and lattice reduction algorithms to prevent non-physical deformations of the simulation box;

c) Fixed the overlap check algorithm allowing overlapping configurations (the disk-rim test has been changed);

d) Added an instruction which undoes any rotation of the simulation box after performing a lattice reduction;

e) Added an option to enable the calculation of the intermolecular potential (perturbed system only). Molecular configurations are still defined by the purely repulsive hard-core potential.
