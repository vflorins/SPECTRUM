/*!
\file simulation.hh
\brief Declares application level classes to perform a complete simulation from start to finish
\author Juan G Alonso Guzman
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SIMULATION_HH
#define SPECTRUM_SIMULATION_HH

#include "common/mpi_config.hh"
#include "common/workload_manager.hh"
#include "server_config.hh"
//#include "trajectory.hh"
#include "background_base.hh"
#include "diffusion_base.hh"
#include "distribution_base.hh"
#include "initial_base.hh"
#include "boundary_base.hh"
#include "trajectory_base.hh"
#include <memory>
#include <chrono>

namespace Spectrum {

//! Whether to print the last trajectory
const bool print_last_trajectory = false;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationWorker (base) class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Application level base class describing a compete simulation process controlled from without
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
template <typename Trajectory_>
class SimulationWorker {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using BackgroundBase = BackgroundBase<Fields>;

   using DistributionBase = DistributionBase<Trajectory>;
   using DiffusionBase = DiffusionBase<Trajectory>;
   using BoundaryBase = BoundaryBase<Trajectory>;
   using InitialBase = InitialBase<Trajectory>;

protected:

//! Whether this is a parallel or a serial run
   bool is_parallel = false;

//! Particle's specie
   unsigned int specie = 0;

//! Random number generator object
   std::shared_ptr<RNG> rng;

//! Local distribution objects
   std::vector<std::shared_ptr<DistributionBase>> local_distros;

//! Trajectory object
   std::unique_ptr<Trajectory> trajectory = nullptr;

//! Number of trajectories in this batch
   int current_batch_size;

//! Number of batches completed by this worker
   int jobsdone;

//! Shortest simulated trajectory time
   double shortest_sim_time;

//! Longest simulated trajectory time
   double longest_sim_time;

//! Elapsed time (both virtual and real)
   double elapsed_time;

//! Send data to master
   void SendDataToMaster(void);

//! Set up for worker prior to main loop
   void WorkerStart(void);

//! Worker finish tasks (reporting partial cumulative)
   void WorkerFinish(void);

//! Worker Duties
   void WorkerDuties(void);

public:

//! Default constructor
   SimulationWorker(void);

//! Destructor
   virtual ~SimulationWorker() = default;
//   virtual ~SimulationWorker() {std::cerr << "Destructor called\n";};

//! Set distro base filename
   virtual void DistroFileName(const std::string& file_name);

//! Set the particle count and the size of one batch
   virtual void SetTasks(int n_traj_in, int batch_size_in, int max_traj_per_worker_in = 0);

//! Get name (type) of trajectory
   std::string GetTrajectoryName(void) const;

//! Get number of completed batches
   int GetJobsDone(void) const {return jobsdone;};

//! Set the particle specie
   void SetSpecie(unsigned int specie_in);

//! Add a background object (passthrough to trajectory)
   void AddBackground(const BackgroundBase& background_in, const DataContainer& container_in, const std::string& fname_pattern_in = "");

   // TODO: experiment
//! Add a distribution object
   virtual void AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in) {
         local_distros.push_back(distribution_in.Clone());
         local_distros.back()->SetSpecie(specie);
         local_distros.back()->SetupObject(container_in);
         trajectory->ConnectDistribution(local_distros.back());
         PrintMessage(__FILE__, __LINE__, "Distribution object added", MPI_Config::is_master);
   }


   // TODO: experiment
//! Add boundary condition object (passthrough to trajectory)
   virtual void AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in) {
      trajectory->AddBoundary(boundary_in, container_in);
      PrintMessage(__FILE__, __LINE__, "Boundary condition added", MPI_Config::is_master);
   }

// TODO: experiment
//! Add initial condition object (passthrough to trajectory)
   virtual void AddInitial(const InitialBase& initial_in, const DataContainer& container_in) {
      trajectory->AddInitial(initial_in, container_in);
      PrintMessage(__FILE__, __LINE__, "Initial condition added", MPI_Config::is_master);
   }

// TODO: experiment
//! Add diffusion object (passthrough to trajectory)
   virtual void AddDiffusion(const DiffusionBase& diffusion_in, const DataContainer& container_in) {
      trajectory->AddDiffusion(diffusion_in, container_in);
      PrintMessage(__FILE__, __LINE__, "Diffusion model added", MPI_Config::is_master);
   }

//! Restore distribution (stub)
   virtual void RestoreDistro(int distro);

//! Print the reduced distribution (stub)
   virtual void PrintDistro1D(int distro, int ijk, const std::string& file_name, bool phys_units) const;

//! Print the reduced distribution in 2D (stub)
   virtual void PrintDistro2D(int distro, int ijk1, int ijk2, const std::string& file_name, bool phys_units) const;

//! Print the distribution records (stub)
   virtual void PrintRecords(int distro, const std::string& file_name, bool phys_units) const;

//! Main simulation loop
   virtual void MainLoop(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationBoss (derived) class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A boss class
\author Juan G Alonso Guzman
*/
template <typename Trajectory_>
class SimulationBoss : public SimulationWorker<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using BackgroundBase = BackgroundBase<Fields>;

   using DistributionBase = DistributionBase<Trajectory>;
   using DiffusionBase = DiffusionBase<Trajectory>;
   using SimulationWorker = SimulationWorker<Trajectory>;
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

#ifdef NEED_SERVER

//! Server backend object
   std::shared_ptr<ServerBaseBack> server_back = nullptr;

#endif

//! Set up for boss prior to main loop
   void BossStart(void);

//! Boss finish tasks
   void BossFinish(void);

//! Boss duties
   void BossDuties(void);

public:

//! Default constructor
   SimulationBoss(void);

//! Add a background object
   void AddBackground(const BackgroundBase& background_in, const DataContainer& container_in, const std::string& fname_pattern_in = "");

//! Main simulation loop
   void MainLoop(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationMaster (derived) class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A master class
\author Juan G Alonso Guzman
*/
template <typename Trajectory_>
class SimulationMaster : public SimulationBoss<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using BackgroundBase = BackgroundBase<Fields>;
   using TrajectoryBase = TrajectoryBase<Trajectory, Fields>;

   using DistributionBase = DistributionBase<Trajectory>;
   using DiffusionBase = DiffusionBase<Trajectory>;
   using SimulationWorker = SimulationWorker<Trajectory>;
   using SimulationBoss = SimulationBoss<Trajectory>;
   using SimulationWorker::current_batch_size;
   using SimulationWorker::is_parallel;
   using SimulationWorker::specie;
   using SimulationWorker::local_distros;
   using SimulationWorker::elapsed_time;
   using SimulationWorker::shortest_sim_time;
   using SimulationWorker::longest_sim_time;
   using SimulationWorker::trajectory;
   // methods:
   using SimulationWorker::WorkerDuties;

protected:

//! Total number of trajectories in the simulation
   int n_trajectories_total;

//! Number of remaining trajectories
   int n_trajectories;

//! Ratio between completed and total trajectory counts * 100
   int percentage_work_done;

//! Number of active workers in node (workers that have not finished all their tasks)
   int active_workers;

//! Maximum number of trajectories per worker. 0 (default) means no maximum.
   int max_traj_per_worker;

//! Number of trajectories already processed and assigned to each worker (for accounting/debugging purposes)
   std::vector <int> trajectories_assigned;

//! Number of trajectories currently being processed on each CPU
   std::vector <int> worker_processing;

//! Total time spent processing trajectories for each worker
   std::vector <double> time_spent_processing;

//! Base file name for partial distros
   std::string distro_file_name = "distribution_";

//! Flag to signal whether to restore distros or not
   std::vector <bool> restore_distros;

//! Local distribution object
// TODO make this a unique pointer
   std::vector<std::shared_ptr<DistributionBase>> partial_distros;

//! MPI request information for worker availability
   std::unique_ptr<MPI_Request_Info> req_cpuavail = nullptr;

//! Workload manager handler
   Workload_Manager_Handler workload_manager_handler;

//! Simulation start time
   std::chrono::system_clock::time_point sim_start_time;

//! Decrement the number of trajectories remaining
   void DecrementTrajectoryCount(void);

//! Receive data from worker
   void RecvDataFromWorker(int cpu);

//! Set up for master prior to main loop
   void MasterStart(void);

//! Master final tasks (gathering partial cumulatives)
   void MasterFinish(void);

//! Master Duties
   void MasterDuties(void);

public:

//! Default constructor
   SimulationMaster(void);

//! Set distro base filename
   void DistroFileName(const std::string& file_name) override;

//! Set the particle count and the size of one batch
   void SetTasks(int n_traj_in, int batch_size_in, int max_traj_per_worker_in = 0) override;

// TODO: experiment
//! Add a distribution object
   void AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in) override {
      // TODO this should use make_unique instead of shared
      partial_distros.push_back(distribution_in.Clone());
      partial_distros.back()->SetSpecie(specie);
      partial_distros.back()->SetupObject(container_in);
      SimulationWorker::AddDistribution(distribution_in, container_in);

// Preset all restore_distro flags to false
      restore_distros.push_back(false);
   }

//! Print simulation info
   void PrintMPICommsInfo(void);

//! Main simulation loop
   void MainLoop(void) override;

//! Restore distribution
   void RestoreDistro(int distro) override;

//! Print the reduced distribution
   void PrintDistro1D(int distro, int ijk, const std::string& file_name, bool phys_units) const override;

//! Print the reduced distribution in 2D
   void PrintDistro2D(int distro, int ijk1, int ijk2, const std::string& file_name, bool phys_units) const override;

//! Print the distribution records
   void PrintRecords(int distro, const std::string& file_name, bool phys_units) const override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Top level methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


//template <typename Trajectory>
//using Simulation = std::variant<std::unique_ptr<SimulationWorker<Trajectory>>, std::unique_ptr<SimulationBoss<Trajectory>>, std::unique_ptr<SimulationMaster<Trajectory>>>;

////! Generate a complete simulation object
//template <typename Trajectory>
//std::unique_ptr<SimulationWorker<Trajectory>> CreateSimulation(int argc, char** argv);




};

// Something like this is needed for templated classes
#include "simulation.cc"

#endif
