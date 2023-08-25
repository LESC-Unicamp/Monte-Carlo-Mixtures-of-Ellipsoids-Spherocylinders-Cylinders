<p align="right"><b><sub>Version: 1.1.0</sub></b></p>

<p align="center"><b>Authors</b></p>
<p align="center">
Nathan Barros de Souza<br>
Luís Fernando Mercier Franco<br></p>

**Version 1.1.0**

- Fixed an issue when attempting to perform anisotropic volume changes;

- Added box distortion parameters and lattice reduction algorithms to prevent non-physical deformations of the simulation box (anisotropic only);

- Fixed the overlap check algorithm allowing overlapping configurations (the disk-rim test has been severely changed);

- Added an instruction which undoes any rotation of the simulation box after performing a lattice reduction;

- Added an option to enable the calculation of the intermolecular potential (perturbed system only). Molecular configurations are still defined by a purely repulsive hard-core potential.

**Version 1.0.0**

- Creation of repository