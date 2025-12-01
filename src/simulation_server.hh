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
#include "server_batl.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationServer (derived) class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Application level class for server processes
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
*/
template<typename HConfig_, typename Trajectory_>
class SimulationServer : public SimulationWorker<HConfig_, Trajectory_> {
public:

   using HConfig = HConfig_;
   using Trajectory = Trajectory_;
   using Server = Cond<HConfig::background == Config::Background::DataCartesian, ServerCartesian<HConfig>,
         Cond<HConfig::background == Config::Background::DataBATL, ServerBATL<HConfig>,
         ServerNone<HConfig>>>;
   using MPI = MPI<HConfig, HConfig::MPI_enabled>;

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
   int active_local_workers;

   //! Server backend object
   std::shared_ptr<Server> server = std::make_unique<Server>();

//! Set up for server prior to main loop
   void ServerStart(void);

//! Server finish tasks
   void ServerFinish(void);

//! Server duties
   void ServerDuties(void);

public:

//! Default constructor
   SimulationServer(void);

//! Main simulation loop
   void MainLoop(void) override;

};

};

#include "simulation_server.cc"

#endif