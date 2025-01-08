# <p align="center">NPT/NVT-Monte Carlo Simulation of Mixtures of <br>Anisomorphic Hard Convex Bodies: <br>Ellipsoids of revolution, Spherocylinders, and Cylinders</p>
<p align="right"><b><sub>Version: 1.4.0</sub></b></p>

<p align="center"><b>Authors</b></p>
<p align="center">
Nathan Barros de Souza<br>
Luís Fernando Mercier Franco<br></p>

**Version 1.4.0**

- Added compatibility with OpenMP API (parallelization);
- Enhanced trajectory file handling: Box dimensions (cubic and noncubic) and particle properties are now correctly written to .xyz files, ensuring proper readability by OVITO.

**Version 1.3.1**

- Improved readability and performance of the code;

- The intermolecular potentials can now be applied to the reference system.

**Version 1.3.0**

- Added an integer for the adjustment frequency of volume changes;

- Removed the LCG generator of random numbers. User can now choose between Fortran's own random number generator (intrinsic routine called <code>RANDOM_NUMBER</code>) or a generator based on bitwise operations (Golden sampling scheme). Check <a href="https://github.com/MaginnGroup/Cassandra/blob/master/Src/random_generators.f90">Cassandra's repository</a> for more information.

**Version 1.2.1**

- Added a cell lists feature to optimize simulation performance. Check this <a href="https://github.com/Allen-Tildesley/examples/blob/master/link_list_module.f90">repository</a> for more information;

- Added a backup feature to the simulation segment.

**Version 1.1.1**

- Improved readability of the code.

- The <code>Box/</code> subfolder is now a proper folder that holds information on box properties, such as dimensions, volume, and distortion. Valid only for the <i>NPT</i>-simulation. 

**Version 1.1.0**

- Fixed an issue when attempting to perform anisotropic volume changes;

- Added box distortion parameters and lattice reduction algorithms to prevent non-physical deformations of the simulation box (anisotropic only);

- Fixed the overlap check algorithm allowing overlapping configurations (the disk-rim test has been severely changed);

- Added an instruction which undoes any rotation of the simulation box after performing a lattice reduction;

- Added an option to enable the calculation of the intermolecular potential (perturbed system only). Molecular configurations are still defined by a purely repulsive hard-core potential;

- Added an overlap check between cylindrical and spherical particles;

- Spherical particles now must be properly identified by the user.

**Version 1.0.0**

- Creation of repository.
