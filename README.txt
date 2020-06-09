NUP Model

This is a coarse-grained Brownian (overdamped Langevin) dynamics simulation code
with optional implicit hydrodynamic interactions using the Rotne-Prager-Yamakawa
tensor. The model is explained in detail in Chapter 2 Section 3 of [1]. The
program is not parallelized, but includes GPU device code which can speed up the
operations involving the mobility matrix considerably, although the Makefile
will need to be altered to enable it.

The program is currently set up to simulate a single polymer but can be altered
to simulate systems with multiple polymers and nanoparticles. The individual
monomers/particles can have different properties such as LJ interaction radius,
hydrodynamic radius, cohesiveness, and charge. Other properties of the
simulation can be modified via a settings file.

Dependencies:
GSL
Intel MKL
(OPTIONAL) SFML - For visualization

Compile with:
make nup_model_sfml - with visualization included (requires SFML)
make nup_model_temp - with only file output

To run:
nup_model_* [FILE]
- where file is the main settings file (see WorldSettings.h)
-if no command line argument is given, the default settings file
(WorldSettings.conf) in the same directory as the executable will be loaded

Brief descriptions of the classes/files:
main.cpp - initialization and FILE/SCREEN output
AA - load Particle properties
CubeSpace - cell list to find neighbouring Particles for interactions
Force - definitions of bond, cohesive, repulsive, charge, noise, external forces
Matrix3 - 3x3 Matrix
NormalDistribution - wrapper for random number generation
Particle - stores position and properties of individual Particles
PeriodicBoundary - 
Polymer - initialize a polymer conformation and handle bond forces
Vector3D - 3 Component Vector
World - setup, simulation update, and measurement functions
WorldSettings - load main settings file

Keyboard Controls (in SCREEN mode):
Q - quit
Arrow Keys - move the polymer/camera
WASD Keys - rotate the polymer
Mouse Wheel - zoom in/out

Signals:
TODO

References:
[1] Vovk, A. (2019). Coarse Grained Modeling of Intrinsically Disordered Protein
Structures and Dynamics (Doctoral dissertation). University of Toronto, Canada.
