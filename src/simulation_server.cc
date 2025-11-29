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
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/25/2022
\param[in] background_in    Background object for type recognition
\param[in] container_in     Data container for initializating the background object
\param[in] fname_pattern_in File naming pattern for the server
*/
template <typename HConfig, typename Trajectory>
void SimulationServer<HConfig, Trajectory>::AddBackground(const Background& background_in, const DataContainer& container_in, const std::string& fname_pattern_in)
{
// Create a unique server backend object based on the user preference stored in "server_config.hh".
#ifdef NEED_SERVER
   server_back = std::make_unique<ServerBackType>(fname_pattern_in);
#else
// This is required if the server is also a worker
   SimulationWorker::AddBackground(background_in, container_in);
#endif
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
#ifdef NEED_SERVER
   active_local_workers = MPI::workers_in_node;
   server_back->ServerStart();
#else
// Signal the master that this CPU is available to do work if server is worker
   if (MPI::is_worker) WorkerStart();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationServer<HConfig, Trajectory>::ServerFinish(void)
{
// Stop the server backend
#ifdef NEED_SERVER
   server_back->ServerFinish();

// Print status message that server left simulation
#ifdef GEO_DEBUG
   std::cerr << "Server with rank " << MPI::server_comm_rank << " exited simulation." << std::endl;
#endif

#else
// Report partial cumulatives if server is worker
   if (MPI::is_worker) WorkerFinish();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
template <typename HConfig, typename Trajectory>
void SimulationServer<HConfig, Trajectory>::ServerDuties(void)
{
#ifdef NEED_SERVER // This needs to be here to avoid a compilation error when NEED_SERVER is not defined
   int workers_stopped;
   workers_stopped = server_back->ServerFunctions();
   active_local_workers -= workers_stopped;
#endif
};

/*!
\author Juan G Alonso Guzman
\date 10/25/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationServer<HConfig, Trajectory>::MainLoop(void)
{
   ServerStart();

#ifdef NEED_SERVER
   while (active_local_workers) ServerDuties();
#else
   if (MPI::is_worker) {
      while (current_batch_size) WorkerDuties();
   };
#endif

   ServerFinish();
};


};
