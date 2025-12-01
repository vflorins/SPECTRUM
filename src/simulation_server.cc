/*!
\file simulation_server.cc
\brief Implements application level classes for server processes
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "simulation_server.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationServer methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/08/2022
*/
template <typename HConfig, typename Trajectory>
SimulationServer<HConfig, Trajectory>::SimulationServer(void)
      : SimulationWorker()
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/26/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationServer<HConfig, Trajectory>::ServerStart(void)
{
// Set the number of workers that are initially active in our node and start the server backend
   if constexpr (HConfig::need_server()) {
      active_local_workers = MPI::workers_in_node;
      server->ServerStart();
   }
   else {
// Signal the master that this CPU is available to do work if server is worker
      if (MPI::is_worker) WorkerStart();
   }
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationServer<HConfig, Trajectory>::ServerFinish(void)
{
// Stop the server backend
   if constexpr (HConfig::need_server()) {
      server->ServerFinish();
      if constexpr (HConfig::build_mode == BuildMode::debug) {
// Print status message that server left simulation
         std::cerr << "Server with rank " << MPI::server_comm_rank << " exited simulation." << std::endl;
      }
   }
   else {
// Report partial cumulatives if server is worker
      if (MPI::is_worker) WorkerFinish();
   }
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
template <typename HConfig, typename Trajectory>
void SimulationServer<HConfig, Trajectory>::ServerDuties(void)
{
   if constexpr (HConfig::need_server()) {
      active_local_workers -= server->ServerFunctions();
   }
};

/*!
\author Juan G Alonso Guzman
\date 10/25/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationServer<HConfig, Trajectory>::MainLoop(void)
{
   ServerStart();

   if constexpr (HConfig::need_server()) {
      while (active_local_workers) ServerDuties();
   }
   else {
      if (MPI::is_worker) {
         while (current_batch_size) WorkerDuties();
      };
   }

   ServerFinish();
};


};
