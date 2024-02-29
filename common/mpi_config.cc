/*!
\file mpi_config.cc
\brief Implements classes responsible for creating MPI various communicators and determining hardware configuration, as well as aiding in non-blocking communication
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "mpi_config.hh" 

#ifdef USE_MPI

#include "definitions.hh"
#include "print_warn.hh"
#include <algorithm>
#include <fstream>

namespace Spectrum {

// Example 1: 3 nodes, 4 proc/node, ALLOW_MASTER_BOSS, ALLOW_BOSS_WORKER
//
// node_comm:     0       1       2       3       4       5       6       7       8       9      10      11
//                |       |       |       |       |       |       |       |       |       |       |       |
//                |_______|_______|_______|       |_______|_______|_______|       |_______|_______|_______|
//
// boss_comm:     0                               4                               8
//                |                               |                               |
//                |_______________________________|_______________________________|
//
// work_comm:     0       1       2       3       4       5       6       7       8       9      10      11
//                |       |       |       |       |       |       |       |       |       |       |       |
//                |_______|_______|_______|_______|_______|_______|_______|_______|_______|_______|_______|



// Example 2: 3 nodes, 4 proc/node, ALLOW_MASTER_BOSS
//
// node_comm:     0       1       2       3       4       5       6       7       8       9      10      11
//                |       |       |       |       |       |       |       |       |       |       |       |
//                |_______|_______|_______|       |_______|_______|_______|       |_______|_______|_______|
//
// boss_comm:     0                               4                               8
//                |                               |                               |
//                |_______________________________|_______________________________|
//
// work_comm:     0       1       2       3               5       6       7               9      10      11
//                |       |       |       |               |       |       |               |       |       |
//                |_______|_______|_______|_______________|_______|_______|_______________|_______|_______|



// Example 3: 3 nodes, 4 proc/node, ALLOW_BOSS_WORKER
//
// node_comm:             1       2       3       4       5       6       7       8       9      10      11
//                        |       |       |       |       |       |       |       |       |       |       |
//                        |_______|_______|       |_______|_______|_______|       |_______|_______|_______|
//
// boss_comm:     0       1                       4                               8
//                |       |                       |                               |
//                |_______|_______________________|_______________________________|
//
// work_comm:     0       1       2       3       4       5       6       7       8       9      10      11
//                |       |       |       |       |       |       |       |       |       |       |       |
//                |_______|_______|_______|_______|_______|_______|_______|_______|_______|_______|_______|



// Example 4: 3 nodes, 4 proc/node
//
// node_comm:             1       2       3       4       5       6       7       8       9      10      11
//                        |       |       |       |       |       |       |       |       |       |       |
//                        |_______|_______|       |_______|_______|_______|       |_______|_______|_______|
//
// boss_comm:     0       1                       4                               8
//                |       |                       |                               |
//                |_______|_______________________|_______________________________|
//
// work_comm:     0               2       3               5       6       7               9      10      11
//                |               |       |               |       |       |               |       |       |
//                |_______________|_______|_______________|_______|_______|_______________|_______|_______|

/*!
\author Vladimir Florinski
\date 06/07/2022
\param[in] argc Number of command line arguments
\param[in] argv Array of command line arguments
*/
MPI_Config::MPI_Config(int argc, char** argv)
{
   int proc;

   MPI_Init(&argc, &argv);

// Global communicator functions
   MPI_Comm_dup(MPI_COMM_WORLD, &glob_comm);
   MPI_Comm_size(glob_comm, &glob_comm_size);
   MPI_Comm_rank(glob_comm, &glob_comm_rank);
   is_master = (glob_comm_rank ? false : true);

// Create the node communicator based on memory sharing. This is a preliminary partition that could be changed depending on whether the node contains the master process.
   MPI_Comm_split_type(glob_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
   MPI_Comm_size(node_comm, &node_comm_size);
   MPI_Comm_rank(node_comm, &node_comm_rank);

// Adjust for multiple bosses per node
#ifdef N_BOSSES_PER_NODE
#if N_BOSSES_PER_NODE > 1
   int sock_comm_size, sock_comm_rank;
   MPI_Comm sock_comm, sock_comm_tmp;

// First split based on socket. This is done to ensure that bosses only communicate with workers within their physical socket. 
   MPI_Comm_split_type(node_comm, OMPI_COMM_TYPE_SOCKET, 0, MPI_INFO_NULL, &sock_comm);
   MPI_Comm_size(sock_comm, &sock_comm_size);
   MPI_Comm_rank(sock_comm, &sock_comm_rank);

// Figure out the number of sockets per node by counting the number of processes with rank 0 in "sock_comm"
   int sock_count_msg = (sock_comm_rank ? 0 : 1);
   MPI_Allreduce(&sock_count_msg, &n_sockets_per_node, 1, MPI_INT, MPI_SUM, node_comm);

// Find number of bosses per socket. This number is forced to be at least 1 and rounded down to be an integer for simplicity.
   n_bosses_per_socket = N_BOSSES_PER_NODE / n_sockets_per_node;
// If there are more sockets than bosses, just pretend N_BOSSES_PER_NODE = n_bosses_per_socket * n_sockets_per_node
   if(n_bosses_per_socket < 1) {
      n_bosses_per_socket = 1;
   }
// If necessary, further split sock_comm for as many bosses as
   else if(n_bosses_per_socket > 1) {
      MPI_Comm_split(sock_comm, (sock_comm_rank % n_bosses_per_socket), 0, &sock_comm_tmp);
      MPI_Comm_free(&sock_comm);
      sock_comm = sock_comm_tmp;
   };

   MPI_Comm_free(&node_comm);
   node_comm = sock_comm;

// Update rank and size in new "node" communicators
   MPI_Comm_size(node_comm, &node_comm_size);
   MPI_Comm_rank(node_comm, &node_comm_rank);
#endif
#endif

// Figure out the number of nodes (or node-like units) by counting the number of processes with rank 0 in "node_comm"
   int node_count_msg = (node_comm_rank ? 0 : 1);
   MPI_Allreduce(&node_count_msg, &n_nodes, 1, MPI_INT, MPI_SUM, glob_comm);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Figure out the node
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int* glob_ranks = new int[n_nodes];

// Post requests on the master to receive global ranks. We know the number of sends (n_nodes), but not which processes will be sending. The answers are received in random order.
   MPI_Request* req_glob_rank;
   if(is_master) {
      req_glob_rank = new MPI_Request[n_nodes];
      for(proc = 0; proc < n_nodes; proc++) {
         MPI_Irecv(glob_ranks + proc, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, glob_comm, req_glob_rank + proc);
      };
   }

// Process 0 on each node sends its global rank to the master (including master to self). At this time it is not known whether process 0 is a boss because a master-boss may not be allowed.
   if(!node_comm_rank) MPI_Send(&glob_comm_rank, 1, MPI_INT, 0, 0, glob_comm);

// Sort the "glob_rank" array and broadcast the sorted array to all processes.
   if(is_master) {
      MPI_Waitall(n_nodes, req_glob_rank, MPI_STATUSES_IGNORE);
      delete[] req_glob_rank;
      std::sort(glob_ranks, glob_ranks + n_nodes);
   };
   MPI_Bcast(glob_ranks, n_nodes, MPI_INT, 0, glob_comm);

// Find the node on process 0 in that node (nodes are numbered based on the place of its rank 0 process in "glob_ranks[]").
   if(!node_comm_rank) {
      my_node = InList(n_nodes, glob_ranks, glob_comm_rank);
      if(my_node < 0) PrintError(__FILE__, __LINE__, "Error in determining node order", true);
   };

// Broadcast the node to all other proceses on the node
   MPI_Bcast(&my_node, 1, MPI_INT, 0, node_comm);

   delete[] glob_ranks;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Boss functions
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef ALLOW_MASTER_BOSS

// If the master cannot be a boss, remove the master from the respective node communicator. On the master "node_comm" will be set to MPI_COMM_NULL.
   MPI_Comm node_comm_tmp;
   if(!my_node) {
      if(node_comm_size < 2) {
         PrintError(__FILE__, __LINE__, "Cannot create a boss different from the master with less than 2 processes in the master's node.", is_master);
         PrintError(__FILE__, __LINE__, "Master will also be boss in its own node.", is_master);
      }
      else {
         MPI_Comm_split(node_comm, (is_master ? MPI_UNDEFINED : 1), node_comm_rank, &node_comm_tmp);
         MPI_Comm_free(&node_comm);
         node_comm = node_comm_tmp;
      };

// The size and rank of the node communicator have changed, so we need to update "node_comm_size" and "node_comm_rank". A non-boss master will set these to -1.
      if(node_comm != MPI_COMM_NULL) {
         MPI_Comm_size(node_comm, &node_comm_size);
         MPI_Comm_rank(node_comm, &node_comm_rank);
      }
      else {
         node_comm_size = -1;
         node_comm_rank = -1;
         my_node = -1;
      };
   };

#endif

   is_boss = (node_comm_rank ? false : true);

// Create a boss communicator. The order is the same as the node order, but if the master is separate, the boss ranks will be off by 1. The master is always rank 0 in "boss_comm" even if separate because it has a lower global rank than the boss on its node.
   MPI_Comm_split(glob_comm, ((is_boss || is_master) ? 1 : MPI_UNDEFINED), my_node, &boss_comm);
   if(boss_comm != MPI_COMM_NULL) {
      MPI_Comm_size(boss_comm, &boss_comm_size);
      MPI_Comm_rank(boss_comm, &boss_comm_rank);
   }
   else {
      boss_comm_size = -1;
      boss_comm_rank = -1;
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Worker functions
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Determine the number of workers and set the worker status. The number could be zero if there is only one process on the node and a boss-worker is not allowed.
   if(is_boss) {

#ifdef ALLOW_BOSS_WORKER
      workers_in_node = node_comm_size;
      is_worker = true;
#else
      workers_in_node = node_comm_size - 1;
      is_worker = false;
#endif
      MPI_Bcast(&workers_in_node, 1, MPI_INT, 0, node_comm);
   }
   else {
      if(node_comm == MPI_COMM_NULL) {
         is_worker = false;
         workers_in_node = -1;
      }
      else {
         is_worker = true;
         MPI_Bcast(&workers_in_node, 1, MPI_INT, 0, node_comm);
      };
   };

// Create a worker communicator that includes the workers and the master. If boss-workers are not allowed, boss processes will be excluded.
   MPI_Comm_split(glob_comm, ((is_master || is_worker) ? 1 : MPI_UNDEFINED), glob_comm_rank, &work_comm);
   if(work_comm != MPI_COMM_NULL) {
      MPI_Comm_size(work_comm, &work_comm_size);
      MPI_Comm_rank(work_comm, &work_comm_rank);
   }
   else {
      work_comm_size = -1;
      work_comm_rank = -1;
   };

// Post requests on the master to receive worker counts.
   MPI_Request* req_worker_count;
   if(is_master) {
      workers_per_node = new int[n_nodes];
      req_worker_count = new MPI_Request[n_nodes];
      for(proc = 0; proc < n_nodes; proc++) {
         MPI_Irecv(workers_per_node + proc, 1, MPI_INT, (is_boss ? proc : proc + 1), MPI_ANY_TAG, boss_comm, req_worker_count + proc);
      };
   };

// Each boss sends its number of workers to the master (including master to self if allowed)
   if(is_boss) MPI_Send(&workers_in_node, 1, MPI_INT, 0, 0, boss_comm);

// Wait for worker count to be received, count the total
   if(is_master) {
      MPI_Waitall(n_nodes, req_worker_count, MPI_STATUSES_IGNORE);
      delete[] req_worker_count;
      n_workers = 0;
      for(proc = 0; proc < n_nodes; proc++) n_workers += workers_per_node[proc];
   };

// Broadcast the number of workers to all
   MPI_Bcast(&n_workers, 1, MPI_INT, 0, glob_comm);
};

/*!
\author Vladimir Florinski
\date 06/07/2022
*/
MPI_Config::~MPI_Config()
{
   if(is_master) delete[] workers_per_node;
   if(work_comm != MPI_COMM_NULL) MPI_Comm_free(&work_comm);
   if(boss_comm != MPI_COMM_NULL) MPI_Comm_free(&boss_comm);
   if(node_comm != MPI_COMM_NULL) MPI_Comm_free(&node_comm);
   MPI_Comm_free(&glob_comm);

   MPI_Finalize();
};

/*!
\author Juan G Alonso Guzman
\date 04/27/2021
\param[in] size Size of communicator in which non-blocking receives will be called
*/
MPI_Request_Info::MPI_Request_Info(int size)
{
   mpi_req = new MPI_Request[size];
   cpu_rank = new int[size];
};

/*!
\author Juan G Alonso Guzman
\date 04/27/2021
*/
MPI_Request_Info::~MPI_Request_Info()
{
   delete[] mpi_req;
   delete[] cpu_rank;
}

/*!
\author Vladimir Florinski
\date 09/17/2020
*/
void TestMPIConfig(void)
{
   MPI_Config mpi_config(0, nullptr);

// Each process prints its config to a different file
   std::string filename = "MPI_CONFIG/mpi_config_" + std::to_string(mpi_config.glob_comm_rank);
   std::ofstream infofile;
   infofile.open(filename.c_str(), std::ofstream::out);

   infofile << "--------------------------------------------------------------------------------\n";
   infofile << "Size of global communicator: " << mpi_config.glob_comm_size << std::endl;
   infofile << "Rank in global communicator: " << mpi_config.glob_comm_rank << std::endl;
   infofile << "Master status: " << (mpi_config.is_master ? "YES" : "NO") << std::endl;
   infofile << std::endl;

   infofile << "Number of nodes: " << mpi_config.n_nodes << std::endl;
   infofile << "Size of boss communicator: " << mpi_config.boss_comm_size << std::endl;
   infofile << "Rank in boss communicator: " << mpi_config.boss_comm_rank << std::endl;
   infofile << "Boss status: " << (mpi_config.is_boss ? "YES" : "NO") << std::endl;
   infofile << std::endl;

   infofile << "This is node: " << mpi_config.my_node << std::endl;
   infofile << "Size of node communicator: " << mpi_config.node_comm_size << std::endl;
   infofile << "Rank in node communicator: " << mpi_config.node_comm_rank << std::endl;
   infofile << "Number of sockets per node: " << mpi_config.n_sockets_per_node << std::endl;
   infofile << "Number of bosses per socket: " << mpi_config.n_bosses_per_socket << std::endl;
   infofile << std::endl;

   infofile << "Total number of workers: " << mpi_config.n_workers << std::endl;
   infofile << "Number of workers in this node: " << mpi_config.workers_in_node << std::endl;
   infofile << "Size of worker communicator: " << mpi_config.work_comm_size << std::endl;
   infofile << "Rank in worker communicator: " << mpi_config.work_comm_rank << std::endl;
   infofile << "Worker status: " << (mpi_config.is_worker ? "YES" : "NO") << std::endl;
   infofile << std::endl;

   if(mpi_config.is_master) {
      for(auto node = 0; node < mpi_config.n_nodes; node++) {
         infofile << "Number of workers in node " << node << ": " << mpi_config.workers_per_node[node] << std::endl;
      };
   };
   infofile << "--------------------------------------------------------------------------------\n";
   
   infofile.close();
};

};

#endif
