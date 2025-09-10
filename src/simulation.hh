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
#include "traj_config.hh"
#include <memory>
#include <chrono>

namespace Spectrum {

//! Whether or not to print each worker's last trajectory
const bool print_last_trajectory = false;

//! Whether or not to print each failed trajectory
const bool print_failed_trajectory = false;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationWorker (base) class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Application level base class describing a compete simulation process controlled from without
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
class SimulationWorker {

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
   std::unique_ptr<TrajectoryBase> trajectory = nullptr;

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

//! Add a distribution object
   virtual void AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in);

//! Add a background object (passthrough to trajectory)
   virtual void AddBackground(const BackgroundBase& background_in, const DataContainer& container_in, const std::string& fname_pattern_in = "");

//! Add boundary condition object (passthrough to trajectory)
   void AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in);

//! Add initial condition object (passthrough to trajectory)
   void AddInitial(const InitialBase& initial_in, const DataContainer& container_in);

//! Add diffusion object (passthrough to trajectory)
   void AddDiffusion(const DiffusionBase& diffusion_in, const DataContainer& container_in);

//! Add source object (passthrough to trajectory)
   void AddSource(const SourceBase& source_in, const DataContainer& container_in);

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
class SimulationBoss : public SimulationWorker {

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
   void AddBackground(const BackgroundBase& background_in, const DataContainer& container_in, const std::string& fname_pattern_in = "") override;

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
class SimulationMaster : public SimulationBoss {

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

//! Add a distribution object
   void AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in) override;

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

//! Generate a complete simulation object
std::unique_ptr<SimulationWorker> CreateSimulation(int argc, char** argv);

};

#endif
