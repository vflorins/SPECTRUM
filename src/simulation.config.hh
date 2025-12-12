/*!
\file simulation.config.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM test particle trajectory simulation
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SIMULATION_CONFIG_HH
#define SPECTRUM_SIMULATION_CONFIG_HH

#include "common/compiletime_lists.hh"
#include "background.config.hh"
#include "trajectory.config.hh"
#include "diffusion.config.hh"

namespace Spectrum {

/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM test particle trajectory distribution solver simulation
\author Lucius Schoenbaum
\date 11/25/2025
*/
struct SimulationConfig {

   static constexpr auto specieid = SpecieId::proton_core;

   static constexpr auto build_mode = BuildMode::debug;

   // Whether to print the last trajectory
   static constexpr auto print_last_trajectory = false;

// Whether there is a supervisor process. This does not require there to be server processes.
   static constexpr auto supervisor = false;

// These can also be set from the command line (argv) at runtime.
   static constexpr auto num_trajectories = 1;
   static constexpr auto batch_size = 1;
   static constexpr auto max_trajectories_per_worker = 1;

};

}

#endif