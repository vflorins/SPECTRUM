# SPECTRUM

Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes (SPECTRUM) is a suite of scientific numerical simulation codes designed to efficiently perform test-particle simulations, classical or relativistic, in virtually any astrophysical environment. In addition, ongoing development seeks to extend the utility of SPECTRUM by equipping it with MHD simulation capabilities. 

**Test-particle Simulations**

Regarding test-particle simulations, SPECTRUM was designed with two main purposes in mind. The first is to compute real, physical trajectories of charged particles by specifying their equations of motion subject to known background fields. The second is to solve partial differential equations (PDEs) describing the diffusive transport of particle distribution functions via the *method of stochastic characteristics*. In a nutshell, this is a Monte Carlo (MC) approach that consists in finding a suitable stochastic differential equation (SDE) whose solution has an expected value that solves the PDE in question when weighed appropriately to account for the initial/boundary conditions, sources/sinks, and linear growth terms.

**MHD Simulations**

This functionality of SPECTRUM is currently under development.

**General File Structure**

The 'common' directory contains fundamental classes and routines used throughout the software. The 'src' directory contains the source code used for the test-particle simulations. The 'fluid', 'geodesic', and 'geometry' directories hold the code relevant to the MHD simulations (under development). The 'benchmarks' directory houses numerous testing and benchmarking files used to verify certain modules are working properly and gauge the software's performance. Finally, the 'runs', and 'tests' directories are meant for the user to implement their own simulations and feature some example files to help them get started.

**Compiling and Running the Code**

SPECTRUM was designed primarily for devices running popular Linux distros, such as Debian or Fedora, and uses the GNU Automake tool for compilation. Note that a C++ compiler (c++17 or above), a properly installed MPI library (e.g. Open MPI or MPICH), and the GNU Scientic Library (GSL) of mathematical tools are pre-requisites for compiling and running SPECTRUM codes.

To configure the code from a fresh download navigate to the main SPECTRUM directory and execute the following commands in terminal:

`autoreconf`

`automake --add-missing`

`./configure <configure-options>`

The `<configure-options>` are used to specify parameters such as which type of trajectory transport should be simulated or what integration method should be used. Consult the documentation ('SPECTRUM.pdf') for more details.
