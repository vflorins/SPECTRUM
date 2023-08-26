/*!
\file simulation.cc
\brief Implements application level classes to perform a complete simulation from start to finish
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Keyvan Ghanbari

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "simulation.hh"
#include "common/print_warn.hh"
#include <numeric>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationWorker methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 06/09/2022
*/
SimulationWorker::SimulationWorker(std::shared_ptr<MPI_Config> mpi_config_in)
{
// Analyze our parallel setup and determine if this is a serial or parallel run.
   mpi_config = mpi_config_in;
   if(mpi_config->n_workers < 1) {
      PrintError(__FILE__, __LINE__, "No worker processes available", mpi_config->is_master);
      is_parallel = false;
   }
   else is_parallel = true;

// Some ambiguity exists for "n_workers" equal to 1 because there can be a separate master, which counts as a parallel run.
#ifdef ALLOW_MASTER_BOSS 
#ifdef ALLOW_BOSS_WORKER
   if(mpi_config->n_workers == 1) is_parallel = false;
#endif
#endif

// Create a common RNG object
   rng = std::make_shared<RNG>(time(NULL) + mpi_config->glob_comm_rank);
#ifdef GEO_DEBUG
   std::cerr << "process " << mpi_config->glob_comm_rank << " random seed = " << time(NULL) + mpi_config->glob_comm_rank << std::endl;
#endif

// Create a unique trajectory object based on the user preference stored in "traj_config.hh".
   trajectory = std::make_unique<TrajectoryType>();
   trajectory->ConnectRNG(rng);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/01/2022
\param[in] file_name
*/
void SimulationWorker::DistroFileName(const std::string& file_name)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/08/2022
\param[in] n_traj_in     Unused
\param[in] batch_size_in Unused
*/
void SimulationWorker::SetTasks(int n_traj_in, int batch_size_in)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2022
\return True is this process is the master
*/
bool SimulationWorker::IsMaster(void)
{
   return mpi_config->is_master;
};

/*!
\author Juan G Alonso Guzman
\date 10/27/2022
\return name of trajectory type
*/
std::string SimulationWorker::GetTrajectoryName(void) const
{
   return trajectory->GetName();
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
\return Number of batches processed by this worker
*/
int SimulationWorker::GetJobsDone(void) const
{
   return jobsdone;
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] specie_in Index of the particle species defined in physics.hh
*/
void SimulationWorker::SetSpecie(unsigned int specie_in)
{
   if((specie_in < 0) || (specie_in >= MAX_PARTICLE_SPECIES)) {
      PrintError(__FILE__, __LINE__, "Invalid particle specie", mpi_config->is_master);
      return;
   };

// The "trajectory" object will set the specie for its sub-classes
   specie = specie_in;
   trajectory->SetSpecie(specie_in);
   PrintMessage(__FILE__, __LINE__, "Particle specie added", mpi_config->is_master);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/27/2022
\param[in] distribution_in Distribution object for type recognition
\param[in] container_in    Data container for initializating the distribution object
*/
void SimulationWorker::AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in)
{
   local_distros.push_back(distribution_in.Clone());
   local_distros.back()->SetSpecie(specie);
   local_distros.back()->SetupObject(container_in);
   trajectory->ConnectDistribution(local_distros.back());
   PrintMessage(__FILE__, __LINE__, "Distribution object added", mpi_config->is_master);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/27/2022
\param[in] background_in    Background object for type recognition
\param[in] container_in     Data container for initializating the background object
\param[in] fname_pattern_in Unused
*/
void SimulationWorker::AddBackground(const BackgroundBase& background_in, const DataContainer& container_in, const std::string& fname_pattern_in)
{
   DataContainer container_mpi(container_in);

#ifdef NEED_SERVER
// Smuggle "mpi_config" in the container through trajectory and background to the server frontend.
   container_mpi.Insert(&mpi_config);
#endif

   trajectory->AddBackground(background_in, container_mpi);
   PrintMessage(__FILE__, __LINE__, "Background object added", mpi_config->is_master);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/27/2022
\param[in] boundary_in  Boundary object for type recognition
\param[in] container_in Data container for initializating the boundary object
*/
void SimulationWorker::AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in)
{
   trajectory->AddBoundary(boundary_in, container_in);
   PrintMessage(__FILE__, __LINE__, "Boundary condition added", mpi_config->is_master);
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] initial_in   Initial object for type recognition
\param[in] container_in Data container for initializating the initial object
*/
void SimulationWorker::AddInitial(const InitialBase& initial_in, const DataContainer& container_in)
{
   trajectory->AddInitial(initial_in, container_in);
   PrintMessage(__FILE__, __LINE__, "Initial condition added", mpi_config->is_master);
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] diffusion_in Diffusion object for type recognition
\param[in] container_in Data container for initializating the diffusion object
*/
void SimulationWorker::AddDiffusion(const DiffusionBase& diffusion_in, const DataContainer& container_in)
{
   trajectory->AddDiffusion(diffusion_in, container_in);
   PrintMessage(__FILE__, __LINE__, "Diffusion model added", mpi_config->is_master);
};

/*!
\author Juan G Alonso Guzman
\date 09/14/2022
*/
void SimulationWorker::SendDataToMaster(void)
{
   int n_events_local, n_records_local;
   size_t distro_size, w_records_size;
   void * distro_addr, * w_records_addr;

// Send distribution data to master and reset distribution
   for(int distro = 0; distro < local_distros.size(); distro++) {
      MPI_Send(local_distros[distro]->GetCountsAddress(), local_distros[distro]->NBins().Prod(), MPI_INT, 0, tag_distrdata, mpi_config->work_comm);
      distro_addr = local_distros[distro]->GetDistroAddress(distro_size);
      MPI_Send(distro_addr, distro_size * local_distros[distro]->NBins().Prod(), MPI_BYTE, 0, tag_distrdata, mpi_config->work_comm);
      n_events_local = local_distros[distro]->NEvents();
      MPI_Send(&n_events_local, 1, MPI_INT, 0, tag_distrdata, mpi_config->work_comm);
      local_distros[distro]->ResetDistribution();

// Send records if they are being kept
      if(local_distros[distro]->GetKeepRecords()) {
         n_records_local = local_distros[distro]->NRecords();
         MPI_Send(&n_records_local, 1, MPI_INT, 0, tag_distrdata, mpi_config->work_comm);
         MPI_Send(local_distros[distro]->GetValuesRecordAddress(), 3*n_records_local, MPI_DOUBLE, 0, tag_distrdata, mpi_config->work_comm);
         w_records_addr = local_distros[distro]->GetWeightsRecordAddress(w_records_size);
         MPI_Send(w_records_addr, w_records_size * n_records_local, MPI_BYTE, 0, tag_distrdata, mpi_config->work_comm);
         local_distros[distro]->ResetRecords();
      };
   };

// Send min/max time data to master (no need to reset)
   MPI_Send(&shortest_sim_time, 1, MPI_DOUBLE, 0, tag_distrdata, mpi_config->work_comm);
   MPI_Send(&longest_sim_time , 1, MPI_DOUBLE, 0, tag_distrdata, mpi_config->work_comm);
   MPI_Send(&elapsed_time     , 1, MPI_DOUBLE, 0, tag_distrdata, mpi_config->work_comm);
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
void SimulationWorker::WorkerStart(void)
{
// Reset quantities
   for(int distro = 0; distro < local_distros.size(); distro++) local_distros[distro]->ResetDistribution();
   jobsdone = 0;
   shortest_sim_time = 1.0E300;
   longest_sim_time = 0.0;
   elapsed_time = 0.0;

// Signal the master that this CPU is available to do work and receive confirmation that more data is needed. Sending data to master does not serve a purpose, but it must be called after the nodes advertise their availability.
   if(is_parallel) {
      MPI_Send(NULL, 0, MPI_INT, 0, tag_cpuavail, mpi_config->work_comm);
      SendDataToMaster();
      MPI_Recv(&current_batch_size, 1, MPI_INT, 0, tag_needmore_MW, mpi_config->work_comm, MPI_STATUS_IGNORE);
   };
};

/*!
\author Juan G Alonso Guzman
\date 05/25/2021
*/

void SimulationWorker::WorkerFinish(void)
{
// Print last trajectory
// TODO possibly pad the number with zeros
   if(print_last_trajectory) {
      std::string traj_name = "trajectory_rank" + std::to_string(mpi_config->work_comm_rank) + ".lines";
      trajectory->PrintTrajectory(traj_name, true, 0x01 | 0x08);
      // trajectory->PrintCSV(traj_name, false);
      // trajectory->InterpretStatus();
   };

#ifdef NEED_SERVER
   trajectory->StopBackground();
#endif

// Print status message that worker left simulation
#ifdef GEO_DEBUG
   std::cerr << "Worker with rank " << mpi_config->work_comm_rank << " exited simulation." << std::endl;
#endif
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/09/2022
*/
void SimulationWorker::WorkerDuties(void)
{
// Iterate over all trajectories in the batch and record average time per trajectory
   int traj_count = 0;
   double traj_elapsed_time, batch_elapsed_time;
   auto start = std::chrono::system_clock::now();
   while(traj_count < current_batch_size) {
      try {
         trajectory->SetStart();
         trajectory->Integrate();
         traj_elapsed_time = trajectory->ElapsedTime();
         if(shortest_sim_time > traj_elapsed_time) shortest_sim_time = traj_elapsed_time;
         if(longest_sim_time < traj_elapsed_time) longest_sim_time = traj_elapsed_time;
         traj_count++;
      }
// TODO Should we treat all exceptions the same?
      catch(std::exception& exception) {
         std::cerr << "Trajectory discarded by worker with rank " << mpi_config->work_comm_rank
                   << ": " << exception.what() << std::endl;
      }
   };
   auto end = std::chrono::system_clock::now();
   batch_elapsed_time = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

// Process batch time depending on whether the code is parallel or not
   if(is_parallel) elapsed_time = batch_elapsed_time;
   else elapsed_time += batch_elapsed_time;

// Increment counter of jobs done by this process
   jobsdone++;

// Signal master that this CPU is available to do work, send batch data, and receive confirmation that more data is needed
   if(is_parallel) {
      MPI_Send(NULL, 0, MPI_INT, 0, tag_cpuavail, mpi_config->work_comm);
      SendDataToMaster();
      MPI_Recv(&current_batch_size, 1, MPI_INT, 0, tag_needmore_MW, mpi_config->work_comm, MPI_STATUS_IGNORE);
   };
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
void SimulationWorker::MainLoop(void)
{
   WorkerStart();
   while(current_batch_size) WorkerDuties();
   WorkerFinish();
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
void SimulationWorker::PrintDistro1D(int distro, int ijk, const std::string& file_name, bool phys_units) const
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
void SimulationWorker::PrintDistro2D(int distro, int ijk1, int ijk2, const std::string& file_name, bool phys_units) const
{
};

/*!
\author Juan G Alonso Guzman
\date 12/02/2022
\param[in] distro which distribution to print
\param[in] file_name filename
\param[in] phys_units whether to print in physical units or not
*/
void SimulationWorker::PrintRecords(int distro, const std::string& file_name, bool phys_units) const
{
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationBoss methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/08/2022
\param[in] mpi_config_in MPI config object
*/
SimulationBoss::SimulationBoss(std::shared_ptr<MPI_Config> mpi_config_in)
              : SimulationWorker(mpi_config_in)
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
void SimulationBoss::AddBackground(const BackgroundBase& background_in, const DataContainer& container_in, const std::string& fname_pattern_in)
{
// Create a unique server backend object based on the user preference stored in "server_config.hh".
#ifdef NEED_SERVER
   server_back = std::make_unique<ServerBackType>(fname_pattern_in);
   server_back->ConnectMPIConfig(mpi_config);
#endif

// This is required if the boss is also a worker
   SimulationWorker::AddBackground(background_in, container_in);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/26/2022
*/
void SimulationBoss::BossStart(void)
{
// Set the number of workers that are initially active in our node and start the server backend
#ifdef NEED_SERVER
   active_local_workers = mpi_config->workers_in_node;
   server_back->ServerStart();
#else
// Signal the master that this CPU is available to do work if boss is worker
   if(mpi_config->is_worker) WorkerStart();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2022
*/
void SimulationBoss::BossFinish(void)
{
// Stop the server backend
#ifdef NEED_SERVER
   server_back->ServerFinish();
#else
// Report partial cumulatives if boss is worker
   if(mpi_config->is_worker) WorkerFinish();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
void SimulationBoss::BossDuties(void)
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
void SimulationBoss::MainLoop(void)
{
   BossStart();

#ifdef NEED_SERVER
   while(active_local_workers) BossDuties();
#else
   if(mpi_config->is_worker) {
      while(current_batch_size) WorkerDuties();
   };
#endif

   BossFinish();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimulationMaster methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
\param[in] mpi_config_in MPI config object
*/
SimulationMaster::SimulationMaster(std::shared_ptr<MPI_Config> mpi_config_in)
                : SimulationBoss(mpi_config_in)
{
// If simulation is parallel, initialize batches assigned array and cpu available requests array
   if(is_parallel) {
      trajectories_assigned.assign(mpi_config->work_comm_size, 0);
      time_spent_processing.assign(mpi_config->work_comm_size, 0.0);
      worker_processing.assign(mpi_config->work_comm_size, 0);

      req_cpuavail = std::make_unique<MPI_Request_Info>(mpi_config->work_comm_size);
      for(int cpu = 0; cpu < mpi_config->work_comm_size; cpu++) req_cpuavail->mpi_req[cpu] = MPI_REQUEST_NULL;
   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/01/2022
\param[in] file_name
*/
void SimulationMaster::DistroFileName(const std::string& file_name)
{
   distro_file_name = file_name;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/08/2022
\param[in] n_traj_in     Total trajectories to simulate
\param[in] batch_size_in Size (trajectory count) of one batch
*/
void SimulationMaster::SetTasks(int n_traj_in, int batch_size_in)
{
   if(n_traj_in < 0) n_traj_in = 0;
   if(batch_size_in < 1) batch_size_in = 1;
   if(batch_size_in > n_traj_in) batch_size_in = n_traj_in;

   n_trajectories_total = n_trajectories = n_traj_in;
   current_batch_size = batch_size_in;
   default_batch_size = batch_size_in;
   percentage_work_done = 0;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 06/08/2022
\param[in] distribution_in Distribution object for type recognition
\param[in] container_in    Data container for initializating the distribution object
*/
void SimulationMaster::AddDistribution(const DistributionBase& distribution_in, const DataContainer& container_in)
{
// TODO this should use make_unique instead of shared
   partial_distros.push_back(distribution_in.Clone());
   partial_distros.back()->SetSpecie(specie);
   partial_distros.back()->SetupObject(container_in);
   SimulationWorker::AddDistribution(distribution_in, container_in);
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2022
*/
void SimulationMaster::PrintMPICommsInfo(void)
{
   std::cerr << "========== MPI Communicators Info ==========" << std::endl;
   std::cerr << "Number of nodes: " << mpi_config->n_nodes << std::endl;
   std::cerr << "Size of global communicator: " << mpi_config->glob_comm_size << std::endl;
   std::cerr << "Size of boss communicator: " << mpi_config->boss_comm_size << std::endl;
   std::cerr << "Size of worker communicator: " << mpi_config->work_comm_size << std::endl;
   std::cerr << "Total number of workers: " << mpi_config->n_workers << std::endl;
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
Only master has this function so only master can decrement batch counter
*/
void SimulationMaster::DecrementTrajectoryCount(void)
{
   long int n_trajectories_remaining, percentage_work_new, distro;
   long int remaining_alloc_time, sim_time_left;

   n_trajectories -= current_batch_size;
   if(n_trajectories < current_batch_size) current_batch_size = n_trajectories;
   
#ifdef GEO_DEBUG
   std::cerr << "Trajectories left unassigned: " + std::to_string(n_trajectories) << std::endl;
#endif

   if(is_parallel) {
      n_trajectories_remaining = n_trajectories + std::accumulate(worker_processing.begin(), worker_processing.end(), 0);
      percentage_work_new = 100 - (100 * n_trajectories_remaining) / n_trajectories_total; 
   }
   else {
      percentage_work_new = 100 - (100 * n_trajectories) / n_trajectories_total;
   };   

   if(percentage_work_new > percentage_work_done) {
      percentage_work_done = percentage_work_new;

// Report the percentage of work done
      std::cerr << std::endl;
      std::cerr << "Checkpoint reached: " << percentage_work_done << "% of work completed.\n";
      if(is_parallel) {
         std::cerr << "Best performing process: "  << *std::max_element(trajectories_assigned.begin() + 1, trajectories_assigned.end())
                   << " trajectories\n";
         std::cerr << "Worst performing process: " << *std::min_element(trajectories_assigned.begin() + 1, trajectories_assigned.end())
                   << " trajectories\n";
      };

// Save the partial distributions
      for(distro = 0; distro < local_distros.size(); distro++) {
         local_distros[distro]->Dump(distro_file_name + std::to_string(distro) + ".out");
      };

// Estimate the remaining silumation time
      if(is_parallel) {
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
      if(remaining_alloc_time == -1) std::cerr << "No workload manager. Unbounded time allocation." << std::endl;
      else std::cerr << "Remaining time allocated for simulation: " << std::setw(20) << remaining_alloc_time << " seconds." << std::endl;
      std::cerr << std::endl;
   };
};

/*!
\author Juan G Alonso Guzman
\date 09/14/2022
\param[in] cpu which cpu from which to receive data
*/
void SimulationMaster::RecvDataFromWorker(int cpu)
{
   int n_events_partial, n_records_partial;
   double shortest_sim_time_cpu, longest_sim_time_cpu;
   size_t distro_size, w_records_size;
   void * distro_addr, * w_records_addr;

// Receive partial distros and add it to cumulative distros
   for(int distro = 0; distro < local_distros.size(); distro++) {
      MPI_Recv(partial_distros[distro]->GetCountsAddress(), partial_distros[distro]->NBins().Prod(),
               MPI_INT, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
      distro_addr = partial_distros[distro]->GetDistroAddress(distro_size);
      MPI_Recv(distro_addr, distro_size * partial_distros[distro]->NBins().Prod(),
               MPI_BYTE, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
      MPI_Recv(&n_events_partial, 1, MPI_INT, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
      partial_distros[distro]->SetNEvents(n_events_partial);
      *local_distros[distro] += *partial_distros[distro];

// Receive records if they are being kept
      if(local_distros[distro]->GetKeepRecords()) {
         MPI_Recv(&n_records_partial , 1, MPI_INT, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
         partial_distros[distro]->SetNRecords(n_records_partial);
         MPI_Recv(partial_distros[distro]->GetValuesRecordAddress(), 3*n_records_partial,
                  MPI_DOUBLE, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
         w_records_addr = partial_distros[distro]->GetWeightsRecordAddress(w_records_size);
         MPI_Recv(w_records_addr, w_records_size * n_records_partial,
                  MPI_BYTE, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
         local_distros[distro]->CopyRecords(*partial_distros[distro]);
      };
   };

// Receive min/max simulated time data and process
   MPI_Recv(&shortest_sim_time_cpu, 1, MPI_DOUBLE, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
   MPI_Recv(&longest_sim_time_cpu , 1, MPI_DOUBLE, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
   if(shortest_sim_time_cpu < shortest_sim_time) shortest_sim_time = shortest_sim_time_cpu;
   if(longest_sim_time_cpu > longest_sim_time) longest_sim_time = longest_sim_time_cpu;
   MPI_Recv(&elapsed_time, 1, MPI_DOUBLE, cpu, tag_distrdata, mpi_config->work_comm, MPI_STATUS_IGNORE);
   time_spent_processing[cpu] += elapsed_time;
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
void SimulationMaster::MasterStart(void)
{
// Print info message
   PrintMPICommsInfo();

// Set the number of workers that are initially active
   active_workers = mpi_config->n_workers;

// Detect workload manager
   workload_manager_handler.DetectManager();

// Get simulation start time
   sim_start_time = std::chrono::system_clock::now();

// Reset quantities
   for(int distro = 0; distro < local_distros.size(); distro++) {
      local_distros[distro]->ResetDistribution();
      partial_distros[distro]->ResetDistribution();
   };
   shortest_sim_time = 1.0E300;
   longest_sim_time = 0.0;
   elapsed_time = 0.0;

// Post an initial receive for workers to respond with availability
   if(is_parallel) {
      for(int cpu = 1; cpu < mpi_config->work_comm_size; cpu++) {
         MPI_Irecv(NULL, 0, MPI_INT, cpu, tag_cpuavail, mpi_config->work_comm, req_cpuavail->mpi_req + cpu);
      };
   };

   std::cerr << "Trajectories left: " + std::to_string(n_trajectories) << std::endl;
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
*/
void SimulationMaster::MasterDuties(void)
{
   int cpu, cpu_ind;

// Service the CPU available requests from all workers
   MPI_Testsome(mpi_config->work_comm_size, req_cpuavail->mpi_req, &req_cpuavail->count, req_cpuavail->cpu_rank, MPI_STATUSES_IGNORE);
   for(cpu_ind = 0; cpu_ind < req_cpuavail->count; cpu_ind++) {
      cpu = req_cpuavail->cpu_rank[cpu_ind];

// Get data from worker and tell the worker if more data is needed (i.e. assign a batch)
      RecvDataFromWorker(cpu);
      MPI_Send(&current_batch_size, 1, MPI_INT, cpu, tag_needmore_MW, mpi_config->work_comm);
      if(current_batch_size) trajectories_assigned[cpu] += current_batch_size;
      worker_processing[cpu] = current_batch_size;

// There are still unassigned batches - decrement the counter. We know that this CPU will start a batch, so post an availability request for when it is finished.
      if(n_trajectories) {
         MPI_Irecv(NULL, 0, MPI_INT, cpu, tag_cpuavail, mpi_config->work_comm, req_cpuavail->mpi_req + cpu);
      }

// There are no unassigned batches - this worker will quit since we just sent it a zero "needmore" signal.
      else active_workers--;

      DecrementTrajectoryCount();
   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/09/2022
*/
void SimulationMaster::MasterFinish(void)
{
   int distro, cpu;
   std::chrono::seconds sim_time_elapsed;
   std::chrono::system_clock::time_point sim_current_time;
   PrintMessage(__FILE__, __LINE__, "Simulation completed", mpi_config->is_master);

// Print last trajectory
// TODO possibly pad the number with zeros
   if(print_last_trajectory && mpi_config->is_worker) {
      std::string traj_name = "trajectory_rank" + std::to_string(mpi_config->work_comm_rank) + ".lines";
      trajectory->PrintCSV(traj_name, false);
      trajectory->InterpretStatus();
   };

// Save the final distributions
   for(distro = 0; distro < local_distros.size(); distro++) {
      local_distros[distro]->Dump(distro_file_name + std::to_string(distro) + ".out");
   };

// Print shortest and longest simulated time
   std::cerr << "Shortest simulated trajectory time = " << shortest_sim_time << std::endl;
   std::cerr << "Longest simulated trajectory time = " << longest_sim_time << std::endl;
   if(is_parallel) {
      std::cerr << "Time per trajectory integration:" << std::endl;
      for(cpu = 1; cpu < mpi_config->work_comm_size; cpu++) {
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
void SimulationMaster::MainLoop(void)
{
   MasterStart();

   if(is_parallel) {
      while(active_workers) MasterDuties();
   }
   else if(active_workers) {
      while(current_batch_size) {
         WorkerDuties();
         DecrementTrajectoryCount();
      };
   };

   MasterFinish();
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
void SimulationMaster::PrintDistro1D(int distro, int ijk, const std::string& file_name, bool phys_units) const
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
void SimulationMaster::PrintDistro2D(int distro, int ijk1, int ijk2, const std::string& file_name, bool phys_units) const
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
void SimulationMaster::PrintRecords(int distro, const std::string& file_name, bool phys_units) const
{
   local_distros[distro]->PrintRecords(file_name, phys_units);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Top level methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

std::unique_ptr<SimulationWorker> CreateSimulation(int argc, char** argv)
{
   std::shared_ptr<MPI_Config> mpi_config = std::make_shared<MPI_Config>(argc, argv);
   if(mpi_config->is_master) return std::make_unique<SimulationMaster>(mpi_config);
   else if(mpi_config->is_boss) return std::make_unique<SimulationBoss>(mpi_config);
   else return std::make_unique<SimulationWorker>(mpi_config);
};

};
