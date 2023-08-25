<p align="right"><b><sub>Version: 1.1.0</sub></b></p>

<p align="center"><b>Authors</b></p>
<p align="center">
Nathan Barros de Souza<br>
Luís Fernando Mercier Franco<br></p>

**Version 1.1.0**

a) Fixed an issue when attempting to perform anisotropic volume changes which deforms the shape of the simulation box;

b) Added box distortion parameters and lattice reduction algorithms to prevent non-physical deformations of the simulation box;

c) Fixed the overlap check algorithm allowing overlapping configurations (the disk-rim test has been changed);

d) Added an instruction which undoes any rotation of the simulation box after performing a lattice reduction;

e) Added an option to enable the calculation of the intermolecular potential (perturbed system only). Molecular configurations are still defined by the purely repulsive hard-core potential.
