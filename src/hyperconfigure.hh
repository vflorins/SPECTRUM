/*!
\file hyperconfigure.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM test particle trajectory simulation
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_HYPERCONFIGURE_HH
#define SPECTRUM_HYPERCONFIGURE_HH

#include "common/compiletime_lists.hh"
#include "simulation.config.hh"
#include "background.config.hh"
#include "trajectory.config.hh"
#include "diffusion.config.hh"


namespace Spectrum {

/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM test particle trajectory distribution solver simulation
\author Lucius Schoenbaum
\date 11/25/2025
*/
template <
      BuildMode build_mode_ = BuildMode::release,
      SpecieId specieid_ = SpecieId::proton_core,
      Config::Background background_ = Config::Background::Dipole,
      Config::Trajectory trajectory_ = Config::Trajectory::Lorentz,
      Config::Diffusion diffusion_ = Config::Diffusion::None,
      typename SimulationConfig_ = Default,
      typename BackgroundConfig_ = Default,
      typename TrajectoryConfig_ = Default,
      typename DiffusionConfig_ = Default
>
struct HyperConfigure {

   using SimulationConfig = Cond<std::same_as<SimulationConfig_, Default>, SimulationConfig<>, SimulationConfig_>;
   using BackgroundConfig = Cond<std::same_as<BackgroundConfig_, Default>, BackgroundConfig<background_, specieid_>, BackgroundConfig_>;
   using TrajectoryConfig = Cond<std::same_as<TrajectoryConfig_, Default>, TrajectoryConfig<trajectory_, specieid_>, TrajectoryConfig_>;
   using DiffusionConfig = Cond<std::same_as<DiffusionConfig_, Default>, DiffusionConfig<diffusion_, specieid_>, DiffusionConfig_>;

   static constexpr auto build_mode = build_mode_;
   static constexpr auto specie = Specie<specieid_>();

   static constexpr bool MPI_enabled() {
#ifdef USE_MPI
      return true;
#else
      return false;
#endif
   }

   static constexpr bool data_background() {
      return (background_ == Config::Background::DataCartesian || background_ == Config::Background::DataBATL);
   }

   static constexpr bool need_server() {
      if constexpr (data_background()) return BackgroundConfig::allow_server_worker;
      else return false;
   }

   static constexpr bool numeric_derivatives() {
      bool out = BackgroundConfig::derivative_method == BackgroundOptions::DerivativeMethod::numeric;
      if constexpr (data_background()) {
// If interp_order is -1, values are read directly from data (reader), not interpolated.
// In this case, no gradient computing method except numeric is possible and the derivative_method is ignored.
         out &= BackgroundConfig::server_interp_order == -1;
      }
      return out;
   }

/*
 * The coordinates for ALL backgrounds are: Pos_t, Time_t, with cartesian coordinates for position.
 * However, these coordinates are unused in the *current* implementation.
 * This information is stored here but this is only done in order to make an update easier if this assumption changes.
 */
   using BackgroundCoordinates = Fields<FConfig<specieid_, CoordinateSystem::cartesian>, Time_t, Pos_t>;

};

}


#endif
