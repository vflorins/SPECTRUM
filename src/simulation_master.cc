/*!
\file simulation_master.cc
\brief Declares application level class for the master/supervisor processes
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "simulation_master.hh"

#include <memory>
#include "common/print_warn.hh"
#include <numeric>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationMaster methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
template <typename HConfig, typename Trajectory>
SimulationMaster<HConfig, Trajectory>::SimulationMaster(void)
      : SimulationServer()
{
// SetTasks parameters can be set via template args or explicitly via SetTasks.
// This takes care of the former case.
   SetTasks(HConfig::num_trajectories, HConfig::batch_size, HConfig::max_trajectories_per_worker);
// If simulation is parallel, initialize batches assigned array and cpu available requests array
   if (MPI::is_parallel()) {
      trajectories_assigned.assign(MPI::work_comm_size, 0);
      time_spent_processing.assign(MPI::work_comm_size, 0.0);
      worker_processing.assign(MPI::work_comm_size, 0);

      req_cpuavail = std::make_unique<RequestInfo>(MPI::work_comm_size);
      for (int cpu = 0; cpu < MPI::work_comm_size; cpu++) req_cpuavail->mpi_req[cpu] = MPI_REQUEST_NULL;
   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/01/2022
\param[in] file_name
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::DistroFileName(const std::string& file_name)
{
   distro_file_name = file_name;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 01/04/2024
\param[in] n_traj_in              Total trajectories to simulate
\param[in] batch_size_in          Size (trajectory count) of one batch
\param[in] max_traj_per_worker_in Maximum trajectories per worker, 0 (default) means no maximum
These can 'task' parameters can optionally be set via the constant expression config structure.
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::SetTasks(int n_traj_in, int batch_size_in, int max_traj_per_worker_in)
{
   if (n_traj_in < 0) n_traj_in = 0;
   if ((batch_size_in < 1) || (batch_size_in > n_traj_in)) batch_size_in = fmin(1, n_traj_in);
   if (max_traj_per_worker_in * MPI::n_workers < n_traj_in) max_traj_per_worker_in = 0;

   n_trajectories_total = n_trajectories = n_traj_in;
   current_batch_size = batch_size_in;
   max_traj_per_worker = max_traj_per_worker_in;
   percentage_work_done = 0;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 06/08/2022
\param[in] distribution_in Distribution object for type recognition
\param[in] container_in    Data container for initializating the distribution object
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in)
{
   partial_distros.push_back(distribution_in.Clone());
   partial_distros.back()->SetSpecie(specie);
   partial_distros.back()->SetupObject(container_in);
   SimulationWorker::AddDistribution(distribution_in, container_in);

// Preset all restore_distro flags to false
   restore_distros.push_back(false);
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::PrintMPICommsInfo(void)
{
   std::cerr << "========== MPI Communicators Info ==========" << std::endl;
   std::cerr << "Number of nodes: " << MPI::n_nodes << std::endl;
   std::cerr << "Size of global communicator: " << MPI::glob_comm_size << std::endl;
   std::cerr << "Size of server communicator: " << MPI::server_comm_size << std::endl;
   std::cerr << "Size of worker communicator: " << MPI::work_comm_size << std::endl;
   std::cerr << "Total number of workers: " << MPI::n_workers << std::endl;
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
Only master has this function so only master can decrement batch counter
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::DecrementTrajectoryCount(void)
{
   long int n_trajectories_remaining, percentage_work_new, distro;
   long int remaining_alloc_time, sim_time_left;

   n_trajectories -= current_batch_size;
   if (n_trajectories < current_batch_size) current_batch_size = n_trajectories;

   if constexpr (HConfig::build_mode == BuildMode::debug) {
      std::cerr << "Trajectories left unassigned: " + std::to_string(n_trajectories) << std::endl;
   }

   if (MPI::is_parallel()) {
      n_trajectories_remaining = n_trajectories + std::accumulate(worker_processing.begin(), worker_processing.end(), 0);
      percentage_work_new = 100 - (100 * n_trajectories_remaining) / n_trajectories_total;
   }
   else {
      percentage_work_new = 100 - (100 * n_trajectories) / n_trajectories_total;
   };

   if (percentage_work_new > percentage_work_done) {
      percentage_work_done = percentage_work_new;

// Report the percentage of work done
      std::cerr << std::endl;
      std::cerr << "Checkpoint reached: " << percentage_work_done << "% of work completed.\n";
      if (MPI::is_parallel()) {
         std::cerr << "Best performing process: "  << *std::max_element(trajectories_assigned.begin() + 1, trajectories_assigned.end())
                   << " trajectories\n";
         std::cerr << "Worst performing process: " << *std::min_element(trajectories_assigned.begin() + 1, trajectories_assigned.end())
                   << " trajectories\n";
      };

// Save the partial distributions
      for (distro = 0; distro < local_distros.size(); distro++) {
         local_distros[distro]->Dump(distro_file_name + std::to_string(distro) + ".out");
      };

// Estimate the remaining silumation time
      if (MPI::is_parallel()) {
// time_left [s] = n_traj_rem [traj] x total_time_spent_integ [ms] x 0.001 [s/ms] / traj_completed [traj] / n_active_workers
         sim_time_left = n_trajectories_remaining * std::accumulate(time_spent_processing.begin(), time_spent_processing.end(), 0.0) * 0.001
                         / (std::accumulate(trajectories_assigned.begin(), trajectories_assigned.end(), 0)
                            - std::accumulate(worker_processing.begin(), worker_processing.end(), 0))
                         / (active_workers + (active_workers == 0 ? 1 : 0));
      }
      else {
         sim_time_left = n_trajectories * elapsed_time * 0.001 / (n_trajectories_total - n_trajectories);
      };
      std::cerr << "Approx time left to complete simulation: " << std::setw(20) << sim_time_left << " seconds." << std::endl;

// Estimate whether remaining allocated time is enough for the rest of the work
      remaining_alloc_time = workload_manager_handler.GetRemAllocTime();
      if (remaining_alloc_time == -1) std::cerr << "No workload manager. Unbounded time allocation." << std::endl;
      else std::cerr << "Remaining time allocated for simulation: " << std::setw(20) << remaining_alloc_time << " seconds." << std::endl;
      std::cerr << std::endl;
   };
};

/*!
\author Juan G Alonso Guzman
\date 08/11/2024
\param[in] cpu which cpu from which to receive data
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::RecvDataFromWorker(int cpu)
{
   int n_events_partial, n_records_partial;
   double shortest_sim_time_cpu, longest_sim_time_cpu;
   size_t distro_size, w_records_size;
   void * distro_addr, * w_records_addr;

// Receive partial distros and add it to cumulative distros
   for (int distro = 0; distro < local_distros.size(); distro++) {
      MPI_Recv(partial_distros[distro]->GetCountsAddress(), partial_distros[distro]->NBins().Prod(),
               MPI_INT, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
      distro_addr = partial_distros[distro]->GetDistroAddress(distro_size);
      MPI_Recv(distro_addr, distro_size * partial_distros[distro]->NBins().Prod(),
               MPI_BYTE, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
      MPI_Recv(&n_events_partial, 1, MPI_INT, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
      partial_distros[distro]->SetNEvents(n_events_partial);
      *local_distros[distro] += *partial_distros[distro];

// Receive records if they are being kept
      if (local_distros[distro]->GetKeepRecords()) {
         MPI_Recv(&n_records_partial , 1, MPI_INT, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
         partial_distros[distro]->SetNRecords(n_records_partial);
         MPI_Recv(partial_distros[distro]->GetValuesRecordAddress(), 3*n_records_partial,
                  MPI_DOUBLE, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
         w_records_addr = partial_distros[distro]->GetWeightsRecordAddress(w_records_size);
         MPI_Recv(w_records_addr, w_records_size * n_records_partial,
                  MPI_BYTE, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
         local_distros[distro]->CopyRecords(*partial_distros[distro]);
      };
   };

// Receive min/max simulated time data and process
   MPI_Recv(&shortest_sim_time_cpu, 1, MPI_DOUBLE, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
   MPI_Recv(&longest_sim_time_cpu , 1, MPI_DOUBLE, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
   if constexpr (HConfig::TrajectoryConfig::timeflow == TrajectoryOptions::TimeFlow::forward) {
      if (shortest_sim_time_cpu < shortest_sim_time) shortest_sim_time = shortest_sim_time_cpu;
      if (longest_sim_time_cpu > longest_sim_time) longest_sim_time = longest_sim_time_cpu;
   }
   else {
      if (shortest_sim_time_cpu > shortest_sim_time) shortest_sim_time = shortest_sim_time_cpu;
      if (longest_sim_time_cpu < longest_sim_time) longest_sim_time = longest_sim_time_cpu;
   }
   MPI_Recv(&elapsed_time, 1, MPI_DOUBLE, cpu, MPI::tag::distrdata, MPI::work_comm, MPI_STATUS_IGNORE);
   time_spent_processing[cpu] += elapsed_time;
};

/*!
\author Juan G Alonso Guzman
\date 08/11/2024
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::MasterStart(void)
{
// Print info message
   PrintMPICommsInfo();

// Set the number of workers that are initially active
   active_workers = MPI::n_workers;

// Detect workload manager
   workload_manager_handler.DetectManager();

// Get simulation start time
   sim_start_time = std::chrono::system_clock::now();

// Reset quantities
   for (int distro = 0; distro < local_distros.size(); distro++) {
      local_distros[distro]->ResetDistribution();
      partial_distros[distro]->ResetDistribution();

// Restore distros if requested by user
      if (restore_distros[distro]) local_distros[distro]->Restore(distro_file_name + std::to_string(distro) + ".out");
   };
   if constexpr (HConfig::TrajectoryConfig::timeflow == TrajectoryOptions::TimeFlow::forward) {
      shortest_sim_time = 1.0E300;
   }
   else {
      shortest_sim_time = -1.0E300;
   }
   longest_sim_time = 0.0;
   elapsed_time = 0.0;

// Post an initial receive for workers to respond with availability
   if (MPI::is_parallel()) {
      for (int cpu = 1; cpu < MPI::work_comm_size; cpu++) {
         MPI_Irecv(NULL, 0, MPI_INT, cpu, MPI::tag::cpuavail, MPI::work_comm, req_cpuavail->mpi_req + cpu);
      };
   };

   std::cerr << "Trajectories left: " + std::to_string(n_trajectories) << std::endl;
};

/*!
\author Juan G Alonso Guzman
\date 01/04/2024
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::MasterDuties(void)
{
   int cpu, cpu_ind, next_batch_size;

// Service the CPU available requests from all workers
   MPI_Testsome(MPI::work_comm_size, req_cpuavail->mpi_req, &req_cpuavail->count, req_cpuavail->cpu_rank, MPI_STATUSES_IGNORE);
   for (cpu_ind = 0; cpu_ind < req_cpuavail->count; cpu_ind++) {
      cpu = req_cpuavail->cpu_rank[cpu_ind];

// Get data from worker and tell the worker if more data is needed (i.e. assign a batch)
      RecvDataFromWorker(cpu);
      if (max_traj_per_worker) next_batch_size = (trajectories_assigned[cpu] < max_traj_per_worker ? current_batch_size : 0);
      else next_batch_size = current_batch_size;
      MPI_Send(&next_batch_size, 1, MPI_INT, cpu, MPI::tag::needmore_MW, MPI::work_comm);
      trajectories_assigned[cpu] += next_batch_size;
      worker_processing[cpu] = next_batch_size;

// If this CPU will start a batch, post an availability request for when it is finished.
      if (next_batch_size) MPI_Irecv(NULL, 0, MPI_INT, cpu, MPI::tag::cpuavail, MPI::work_comm, req_cpuavail->mpi_req + cpu);
// There are no unassigned batches - this worker will quit since we just sent it a zero "needmore" signal.
      else active_workers--;
// If there are still unassigned batches - decrement the counter.
      DecrementTrajectoryCount();
   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/09/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::MasterFinish(void)
{
   int distro, cpu;
   std::chrono::seconds sim_time_elapsed;
   std::chrono::system_clock::time_point sim_current_time;
   PrintMessage(__FILE__, __LINE__, "Simulation completed", MPI::is_master);

// Print last trajectory
   if (print_last_trajectory && MPI::is_worker) {
      std::string rank_str = std::to_string(MPI::work_comm_rank);
      std::string size_str = std::to_string(MPI::work_comm_size);
      rank_str.insert(0, size_str.size() - rank_str.size(), '0');
      std::string traj_name = "trajectory_rank_" + rank_str + ".lines";
      trajectory->PrintCSV(traj_name, false, 1);
      trajectory->InterpretStatus();
   };

// Save the final distributions
   for (distro = 0; distro < local_distros.size(); distro++) {
      local_distros[distro]->Dump(distro_file_name + std::to_string(distro) + ".out");
   };

// Print shortest and longest simulated time
   std::cerr << "Shortest simulated trajectory time = " << shortest_sim_time * unit_time_fluid << " s" << std::endl;
   std::cerr << "Longest simulated trajectory time = " << longest_sim_time * unit_time_fluid << " s" << std::endl;
   if (MPI::is_parallel()) {
      std::cerr << "Time per trajectory integration:" << std::endl;
      for (cpu = 1; cpu < MPI::work_comm_size; cpu++) {
         std::cerr << "\tcpu " << cpu << " = " << time_spent_processing[cpu] / trajectories_assigned[cpu] << " ms" << std::endl;
      };
   }
   else {
      std::cerr << "Time per trajectory integration = " << elapsed_time / n_trajectories_total << " ms" << std::endl;
   };

   sim_current_time = std::chrono::system_clock::now();
   sim_time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(sim_current_time - sim_start_time);
   std::cerr << "Total time of this simulation: " << std::setw(20) << sim_time_elapsed.count() << " seconds." << std::endl;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/09/2022
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::MainLoop(void)
{
   MasterStart();

// This is a parallel run in which the master does not do any simulation work, but assigns batches to worker processes.
   if (MPI::is_parallel()) {
      while (active_workers) MasterDuties();
   }
// This is a serial run in which the master process does the work. "active_workers" is checked because it could be 0 from an error in "mpi_config.hh" setup by the user.
   else if (active_workers) {
      while (current_batch_size) {
         WorkerDuties();
         DecrementTrajectoryCount();
      };
   };

   MasterFinish();
};

/*!
\author Juan G Alonso Guzman
\date 11/03/2022
\param[in] distro which distribution to restore
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::RestoreDistro(int distro)
{
   restore_distros[distro] = true;
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
void SimulationMaster<HConfig, Trajectory>::PrintDistro1D(int distro, int ijk, const std::string& file_name, bool phys_units) const
{
   local_distros[distro]->Print1D(ijk, file_name, phys_units);
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
void SimulationMaster<HConfig, Trajectory>::PrintDistro2D(int distro, int ijk1, int ijk2, const std::string& file_name, bool phys_units) const
{
   local_distros[distro]->Print2D(ijk1, ijk2, file_name, phys_units);
};

/*!
\author Juan G Alonso Guzman
\date 12/02/2022
\param[in] distro which distribution to print
\param[in] file_name filename
\param[in] phys_units whether to print in physical units or not
*/
template <typename HConfig, typename Trajectory>
void SimulationMaster<HConfig, Trajectory>::PrintRecords(int distro, const std::string& file_name, bool phys_units) const
{
   local_distros[distro]->PrintRecords(file_name, phys_units);
};


}