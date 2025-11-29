/*!
\file simulation_server.hh
\brief Declares application level class for server processes
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SIMULATION_SERVER_HH
#define SPECTRUM_SIMULATION_SERVER_HH

#include "simulation_worker.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationServer (derived) class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A server class
\author Juan G Alonso Guzman
*/
// todo if !servers_are_workers, does this need to inherit from worker?
template<typename HConfig_, typename Trajectory_>
class SimulationServer : public SimulationWorker<HConfig_, Trajectory_> {
public:

   using HConfig = HConfig_;
   using Trajectory = Trajectory_;
   using Background = Trajectory::Background;
   using MPI = MPI<HConfig>;

   using DistributionBase = DistributionBase<HConfig>;
   using DiffusionBase = DiffusionBase<HConfig>;
   using SimulationWorker = SimulationWorker<HConfig, Trajectory>;
   using SimulationWorker::current_batch_size;
   using SimulationWorker::is_parallel;
   using SimulationWorker::specie;
   // methods:
   using SimulationWorker::WorkerStart;
   using SimulationWorker::WorkerFinish;
   using SimulationWorker::WorkerDuties;

protected:

//! Number of active workers in node
// todo constexpr compute
   int active_local_workers;

#ifdef NEED_SERVER

   //! Server backend object
   std::shared_ptr<ServerBaseBack> server_back = nullptr;

#endif

//! Set up for server prior to main loop
   void ServerStart(void);

//! Server finish tasks
   void ServerFinish(void);

//! Server duties
   void ServerDuties(void);

public:

//! Default constructor
   SimulationServer(void);

//! Add a background object
   void AddBackground(const Background &background_in, const DataContainer &container_in,
                      const std::string &fname_pattern_in = "");

//! Main simulation loop
   void MainLoop(void) override;

};

};

#include "simulation_server.cc"

#endif