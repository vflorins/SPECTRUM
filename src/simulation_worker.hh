/*!
\file simulation_worker.hh
\brief Declares application level class for all worker processes
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SIMULATION_WORKER_HH
#define SPECTRUM_SIMULATION_WORKER_HH

#include "common/mpi_config.hh"
#include "common/workload_manager.hh"
#include "trajectory_base.hh"
#include "utils_numerical_derivatives.hh"
#include "diffusion_base.hh"
#include "distribution_base.hh"
#include "initial_base.hh"
#include "boundary_base.hh"
#include "trajectory_base.hh"


namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationWorker (base) class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Application level base class describing a compete simulation process controlled from without
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
template <typename HConfig_, typename Trajectory_>
class SimulationWorker {
public:

   using HConfig = HConfig_;
   using Trajectory = Trajectory_;
   using Config = HConfig::SimulationConfig;
   using MPI = MPI<HConfig, HConfig::MPI_enabled>;

   using DistributionBase = DistributionBase<HConfig>;
   using BoundaryBase = BoundaryBase<HConfig>;
   using InitialBase = InitialBase<HConfig>;

   static constexpr bool print_last_trajectory = Config::print_last_trajectory;

   static constexpr TrajectoryOptions::TimeFlow timeflow = HConfig::TrajectoryConfig::timeflow;

protected:

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

//! Add a distribution object
   void AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in);

//! Add boundary condition object (passthrough to trajectory)
   void AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in);

//! Add initial condition object (passthrough to trajectory)
   void AddInitial(const InitialBase& initial_in, const DataContainer& container_in);

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

};

#include "simulation_worker.cc"

#endif
