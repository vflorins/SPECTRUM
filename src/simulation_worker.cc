/*!
\file simulation_worker.cc
\brief Declares application level class for all worker processes
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "simulation_worker.hh"

#include <memory>
#include <chrono>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationWorker methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Application level class for all worker processes
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/15/2024
*/
template <typename HConfig, typename Trajectory>
SimulationWorker<HConfig, Trajectory>::SimulationWorker(void)
{
// Check to make sure that there is at least one worker in the simulation. Otherwise, the simulation cannot happen and the program will exit immediately.
   if (MPI::n_workers < 1) {
      PrintError(__FILE__, __LINE__, "No worker processes available.", MPI::is_master);
      PrintError(__FILE__, __LINE__, "Simulation will exit immediately.", MPI::is_master);
      exit(1);
   }

// Create a common RNG object
   rng = std::make_shared<RNG>(time(NULL) + MPI::glob_comm_rank);

   if constexpr (HConfig::build_mode == BuildMode::debug) {
      std::cerr << "process " << MPI::glob_comm_rank << " random seed = " << time(NULL) + MPI::glob_comm_rank << std::endl;
   }

// Create a unique trajectory object based on the user preference stored in "traj_config.hh".
   trajectory = std::make_unique<Trajectory>();
   trajectory->ConnectRNG(rng);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/01/2022
\param[in] file_name
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::DistroFileName(const std::string& file_name)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/08/2022
\param[in] n_traj_in              Unused
\param[in] batch_size_in          Unused
\param[in] max_traj_per_worker_in Unused
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::SetTasks(int n_traj_in, int batch_size_in, int max_traj_per_worker_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/27/2022
\return name of trajectory type
*/
template <typename HConfig, typename Trajectory>
std::string SimulationWorker<HConfig, Trajectory>::GetTrajectoryName(void) const
{
   return trajectory->GetName();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/27/2022
\param[in] distribution_in Distribution object for type recognition
\param[in] container_in    Data container for initializating the distribution object
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in)
{
   local_distros.push_back(distribution_in.Clone());
   local_distros.back()->SetupObject(container_in);
   trajectory->ConnectDistribution(local_distros.back());
   PrintMessage(__FILE__, __LINE__, "Distribution object added", MPI::is_master);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/27/2022
\param[in] boundary_in  Boundary object for type recognition
\param[in] container_in Data container for initializating the boundary object
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in)
{
   trajectory->AddBoundary(boundary_in, container_in);
   PrintMessage(__FILE__, __LINE__, "Boundary condition added", MPI::is_master);
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] initial_in   Initial object for type recognition
\param[in] container_in Data container for initializating the initial object
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::AddInitial(const InitialBase& initial_in, const DataContainer& container_in)
{
   trajectory->AddInitial(initial_in, container_in);
   PrintMessage(__FILE__, __LINE__, "Initial condition added", MPI::is_master);
};

/*!
\author Juan G Alonso Guzman
\date 09/14/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::SendDataToMaster(void)
{
   int n_events_local, n_records_local;
   size_t distro_size, w_records_size;
   void * distro_addr, * w_records_addr;

// Send distribution data to master and reset distribution
   for (int distro = 0; distro < local_distros.size(); distro++) {
      MPI_Send(local_distros[distro]->GetCountsAddress(), local_distros[distro]->NBins().Prod(), MPI_INT, 0, MPI::tag::distrdata, MPI::work_comm);
      distro_addr = local_distros[distro]->GetDistroAddress(distro_size);
      MPI_Send(distro_addr, distro_size * local_distros[distro]->NBins().Prod(), MPI_BYTE, 0, MPI::tag::distrdata, MPI::work_comm);
      n_events_local = local_distros[distro]->NEvents();
      MPI_Send(&n_events_local, 1, MPI_INT, 0, MPI::tag::distrdata, MPI::work_comm);
      local_distros[distro]->ResetDistribution();

// Send records if they are being kept
      if (local_distros[distro]->GetKeepRecords()) {
         n_records_local = local_distros[distro]->NRecords();
         MPI_Send(&n_records_local, 1, MPI_INT, 0, MPI::tag::distrdata, MPI::work_comm);
         MPI_Send(local_distros[distro]->GetValuesRecordAddress(), 3*n_records_local, MPI_DOUBLE, 0, MPI::tag::distrdata, MPI::work_comm);
         w_records_addr = local_distros[distro]->GetWeightsRecordAddress(w_records_size);
         MPI_Send(w_records_addr, w_records_size * n_records_local, MPI_BYTE, 0, MPI::tag::distrdata, MPI::work_comm);
         local_distros[distro]->ResetRecords();
      };
   };

// Send min/max time data to master (no need to reset)
   MPI_Send(&shortest_sim_time, 1, MPI_DOUBLE, 0, MPI::tag::distrdata, MPI::work_comm);
   MPI_Send(&longest_sim_time , 1, MPI_DOUBLE, 0, MPI::tag::distrdata, MPI::work_comm);
   MPI_Send(&elapsed_time     , 1, MPI_DOUBLE, 0, MPI::tag::distrdata, MPI::work_comm);
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::WorkerStart(void)
{
   trajectory->StartBackground();
// Reset quantities
   for (int distro = 0; distro < local_distros.size(); distro++) local_distros[distro]->ResetDistribution();
   jobsdone = 0;
   if constexpr (timeflow == TrajectoryOptions::TimeFlow::forward) {
      shortest_sim_time = 1.0E300;
   }
   else {
      shortest_sim_time = -1.0E300;
   }
   longest_sim_time = 0.0;
   elapsed_time = 0.0;

// Signal the master (with an empty message) that this CPU is available and receive confirmation to do more work.
   if (MPI::is_parallel()) {
      MPI_Send(NULL, 0, MPI_INT, 0, MPI::tag::cpuavail, MPI::work_comm);
      SendDataToMaster();
      MPI_Recv(&current_batch_size, 1, MPI_INT, 0, MPI::tag::needmore_MW, MPI::work_comm, MPI_STATUS_IGNORE);
   };
};

/*!
\author Juan G Alonso Guzman
\date 05/25/2021
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::WorkerFinish(void)
{
// Print last trajectory
   if (print_last_trajectory) {
      std::string rank_str = std::to_string(MPI::work_comm_rank);
      std::string size_str = std::to_string(MPI::work_comm_size);
      rank_str.insert(0, size_str.size() - rank_str.size(), '0');
      std::string traj_name = "trajectory_rank_" + rank_str + ".lines";
//      trajectory->records->PrintTrajectory(traj_name, true, 0x01 | 0x02 | 0x04 | 0x08, 0, 1.0 / unit_time_fluid);
      trajectory->records->PrintCSV(traj_name, false, 1);
      trajectory->InterpretStatus();
   };

   trajectory->StopBackground();

// Print status message that worker left simulation
   if constexpr (HConfig::build_mode == BuildMode::debug) {
      std::cerr << "Worker with rank " << MPI::work_comm_rank << " exited simulation." << std::endl;
   }
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/09/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::WorkerDuties(void)
{
// Iterate over all trajectories in the batch and record average time per trajectory
   int traj_count = 0;
   double traj_elapsed_time, batch_elapsed_time;
   auto start = std::chrono::system_clock::now();
   while (traj_count < current_batch_size) {
      try {
         trajectory->SetStart();
         trajectory->Integrate();
         traj_elapsed_time = trajectory->ElapsedTime();
         if constexpr (timeflow == TrajectoryOptions::TimeFlow::forward) {
            if (shortest_sim_time > traj_elapsed_time) shortest_sim_time = traj_elapsed_time;
            if (longest_sim_time < traj_elapsed_time) longest_sim_time = traj_elapsed_time;
         }
         else {
            if (shortest_sim_time < traj_elapsed_time) shortest_sim_time = traj_elapsed_time;
            if (longest_sim_time > traj_elapsed_time) longest_sim_time = traj_elapsed_time;
         }
         traj_count++;
      }
      catch (std::exception& exception) {
         std::cerr << "Trajectory discarded by worker with rank " << MPI::work_comm_rank
                   << ": " << exception.what() << std::endl;
      }
   };
   auto end = std::chrono::system_clock::now();
   batch_elapsed_time = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

// Process batch time depending on whether the code is parallel or not
   if (MPI::is_parallel()) elapsed_time = batch_elapsed_time;
   else elapsed_time += batch_elapsed_time;

// Increment counter of jobs done by this process
   jobsdone++;

// Signal master that this CPU is available to do work, send batch data, and receive confirmation that more data is needed
   if (MPI::is_parallel()) {
      MPI_Send(NULL, 0, MPI_INT, 0, MPI::tag::cpuavail, MPI::work_comm);
      SendDataToMaster();
      MPI_Recv(&current_batch_size, 1, MPI_INT, 0, MPI::tag::needmore_MW, MPI::work_comm, MPI_STATUS_IGNORE);
   };
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::MainLoop(void)
{
   WorkerStart();
   while (current_batch_size) WorkerDuties();
   WorkerFinish();
};

/*!
\author Juan G Alonso Guzman
\date 11/03/2022
\param[in] distro which distribution to restore
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::RestoreDistro(int distro)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
\param[in] distro which distribution to print
\param[in] ijk which dimension to print (remaining dimensions are collapsed)
\param[in] file_name filename
\param[in] phys_units whether to print in physical units or not
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::PrintDistro1D(int distro, int ijk, const std::string& file_name, bool phys_units) const
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/29/2022
\param[in] distro which distribution to print
\param[in] ijk1 first dimension to print
\param[in] ijk2 second dimension to print (remaining dimension is collapsed)
\param[in] file_name filename
\param[in] phys_units whether to print in physical units or not
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::PrintDistro2D(int distro, int ijk1, int ijk2, const std::string& file_name, bool phys_units) const
{
};

/*!
\author Juan G Alonso Guzman
\date 12/02/2022
\param[in] distro which distribution to print
\param[in] file_name filename
\param[in] phys_units whether to print in physical units or not
*/
template <typename HConfig, typename Trajectory>
void SimulationWorker<HConfig, Trajectory>::PrintRecords(int distro, const std::string& file_name, bool phys_units) const
{
};



}