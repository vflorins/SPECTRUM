/*!
\file simulation.config.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM test particle trajectory simulation
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/


#ifndef SPECTRUM_SIMULATION_CONFIG_HH
#define SPECTRUM_SIMULATION_CONFIG_HH

#include "common/compiletime_lists.hh"
#include "trajectory.config.hh"
#include "background.config.hh"
#include "diffusion.config.hh"
#include "distribution.config.hh"

namespace Spectrum {

/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM test particle trajectory distribution solver simulation
\author Lucius Schoenbaum
\date 09/05/2025
*/
template <
      SpecieId specieid_,
      TrajectoryId trajectoryid_,
      BuildMode build_mode_ = BuildMode::release,
      int num_trajectories_ = 1,
      int batch_size_ = 1,
      typename TrajectoryConfig_ = Default,
      typename BackgroundConfig_ = Default,
      typename DiffusionConfig_ = Default
>
struct SimulationConfig {
   using TrajectoryConfig = TrajectoryConfig_;
   using BackgroundConfig = BackgroundConfig_;
   using DiffusionConfig = DiffusionConfig_;

// The specieid is only needed at compile time, during generation of default config types.
   static constexpr auto specieid = specieid_;
   static constexpr auto trajectoryid = trajectoryid_;
   static constexpr auto build_mode = build_mode_;
   static constexpr auto num_trajectories = num_trajectories_;
   static constexpr auto batch_size = batch_size_;

   static constexpr auto specie = Specie<specieid>();

};

}

#endif
