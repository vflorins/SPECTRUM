/*!
\file simulation_master.hh
\brief Declares application level class for the master/supervisor process
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SIMULATION_MASTER_HH
#define SPECTRUM_SIMULATION_MASTER_HH

#include "simulation_server.hh"


namespace Spectrum {


//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationMaster (derived) class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A master class
\author Juan G Alonso Guzman
*/
template <typename HConfig_, typename Trajectory_>
class SimulationMaster : public SimulationServer<HConfig_, Trajectory_> {
public:

   using HConfig = HConfig_;
   using Trajectory = Trajectory_;
   using Background = Trajectory::Background;
   using MPI = MPI<HConfig>;

   using DistributionBase = DistributionBase<HConfig>;
   using DiffusionBase = DiffusionBase<HConfig>;
   using SimulationWorker = SimulationWorker<HConfig, Trajectory>;
   using SimulationServer = SimulationServer<HConfig, Trajectory>;
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
   std::unique_ptr<RequestInfo> req_cpuavail = nullptr;

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

};

#include "simulation_master.cc"

#endif
