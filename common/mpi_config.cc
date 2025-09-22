/*!
\file MPI_Config::cc
\brief Implements classes responsible for creating MPI various communicators and determining hardware configuration, as well as aiding in non-blocking communication
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "mpi_config.hh" 

#ifdef USE_MPI

#include <algorithm>
#include <fstream>

#include "common/definitions.hh"
#include "common/print_warn.hh"

namespace Spectrum {

// Example 1: 3 nodes, 4 proc/node, ALLOW_MASTER_BOSS, ALLOW_BOSS_WORKER
//
// node_comm:     0       1       2       3       4       5       6       7       8       9      10      11
//                |       |       |       |       |       |       |       |       |       |       |       |
//                |_______|_______|_______|       |_______|_______|_______|       |_______|_______|_______|
//
// server_comm:     0                               4                               8
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
// server_comm:     0                               4                               8
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
// server_comm:     0       1                       4                               8
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
// server_comm:     0       1                       4                               8
//                |       |                       |                               |
//                |_______|_______________________|_______________________________|
//
// work_comm:     0               2       3               5       6       7               9      10      11
//                |               |       |               |       |       |               |       |       |
//                |_______________|_______|_______________|_______|_______|_______________|_______|_______|

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/15/2025
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

// Adjust for multiple servers per node
#ifdef N_BOSSES_PER_NODE
#if N_BOSSES_PER_NODE > 1
   int sock_comm_size, sock_comm_rank;
   MPI_Comm sock_comm, sock_comm_tmp;

// First split based on socket. This is done to ensure that servers only communicate with workers within their physical socket.
   MPI_Comm_split_type(node_comm, OMPI_COMM_TYPE_SOCKET, 0, MPI_INFO_NULL, &sock_comm);
   MPI_Comm_size(sock_comm, &sock_comm_size);
   MPI_Comm_rank(sock_comm, &sock_comm_rank);

// Figure out the number of sockets per node by counting the number of processes with rank 0 in "sock_comm"
   int sock_count_msg = (sock_comm_rank ? 0 : 1);
   MPI_Allreduce(&sock_count_msg, &n_sockets_per_node, 1, MPI_INT, MPI_SUM, node_comm);

// Find number of servers per socket. This number is forced to be at least 1 and rounded down to be an integer for simplicity.
   n_servers_per_socket = N_BOSSES_PER_NODE / n_sockets_per_node;

// If there are more sockets than servers, just pretend N_BOSSES_PER_NODE = n_servers_per_socket * n_sockets_per_node
   if (n_servers_per_socket < 1) {
      n_servers_per_socket = 1;
   }

// If necessary, further split sock_comm for as many servers as
   else if (n_servers_per_socket > 1) {
      MPI_Comm_split(sock_comm, (sock_comm_rank % n_servers_per_socket), 0, &sock_comm_tmp);
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
   if (is_master) {
      req_glob_rank = new MPI_Request[n_nodes];
      for (proc = 0; proc < n_nodes; proc++) {
         MPI_Irecv(glob_ranks + proc, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, glob_comm, req_glob_rank + proc);
      };
   };

// Process 0 on each node sends its global rank to the master (including master to self). At this time it is not known whether process 0 is a server because a master-server may not be allowed.
   if (!node_comm_rank) MPI_Send(&glob_comm_rank, 1, MPI_INT, 0, 0, glob_comm);

// Sort the "glob_rank" array and broadcast the sorted array to all processes.
   if (is_master) {
      MPI_Waitall(n_nodes, req_glob_rank, MPI_STATUSES_IGNORE);
      delete[] req_glob_rank;
      std::sort(glob_ranks, glob_ranks + n_nodes);
   };
   MPI_Bcast(glob_ranks, n_nodes, MPI_INT, 0, glob_comm);

// Find the node on process 0 in that node (nodes are numbered based on the place of its rank 0 process in "glob_ranks[]").
   if (!node_comm_rank) {
      my_node = InList(n_nodes, glob_ranks, glob_comm_rank);
      if (my_node < 0) PrintError(__FILE__, __LINE__, "Error in determining node order", true);
   };

// Broadcast the node to all other proceses on the node
   MPI_Bcast(&my_node, 1, MPI_INT, 0, node_comm);

   delete[] glob_ranks;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Server functions
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef ALLOW_MASTER_BOSS

// If the master cannot be a server, remove the master from the respective node communicator. On the master "node_comm" will be set to MPI_COMM_NULL.
   MPI_Comm node_comm_tmp;
   if (!my_node) {
      if (node_comm_size < 2) {
         PrintError(__FILE__, __LINE__, "Less than 2 processes in the master's node.", is_master);
#ifdef NEED_SERVER
         PrintError(__FILE__, __LINE__, "Cannot create a server different from the master.", is_master);
         PrintError(__FILE__, __LINE__, "Master will also be server in its own node.", is_master);
#else
         PrintError(__FILE__, __LINE__, "Cannot create a worker different from the master.", is_master);
         PrintError(__FILE__, __LINE__, "Master will also be worker in its own node.", is_master);
#endif
      }
      else {
         MPI_Comm_split(node_comm, (is_master ? MPI_UNDEFINED : 1), node_comm_rank, &node_comm_tmp);
         MPI_Comm_free(&node_comm);
         node_comm = node_comm_tmp;
      };

// The size and rank of the node communicator have changed, so we need to update "node_comm_size" and "node_comm_rank". A non-server master will set these to -1.
      if (node_comm != MPI_COMM_NULL) {
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

   is_server = (node_comm_rank ? false : true);

// Create a server communicator. The order is the same as the node order, but if the master is separate, the server ranks will be off by 1. The master is always rank 0 in "server_comm" even if separate because it has a lower global rank than the server on its node.
   MPI_Comm_split(glob_comm, ((is_server || is_master) ? 1 : MPI_UNDEFINED), my_node, &server_comm);
   if (server_comm != MPI_COMM_NULL) {
      MPI_Comm_size(server_comm, &server_comm_size);
      MPI_Comm_rank(server_comm, &server_comm_rank);
   }
   else {
      server_comm_size = -1;
      server_comm_rank = -1;
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Worker functions
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Determine the number of workers and set the worker status. The number could be zero if there is only one process on the node and a server-worker is not allowed.
   if (is_server) {

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
      if (node_comm == MPI_COMM_NULL) {
         is_worker = false;
         workers_in_node = -1;
      }
      else {
         is_worker = true;
         MPI_Bcast(&workers_in_node, 1, MPI_INT, 0, node_comm);
      };
   };

// Create a worker communicator that includes the workers and the master. If server-workers are not allowed, server processes will be excluded.
   MPI_Comm_split(glob_comm, ((is_master || is_worker) ? 1 : MPI_UNDEFINED), glob_comm_rank, &work_comm);
   if (work_comm != MPI_COMM_NULL) {
      MPI_Comm_size(work_comm, &work_comm_size);
      MPI_Comm_rank(work_comm, &work_comm_rank);
   }
   else {
      work_comm_size = -1;
      work_comm_rank = -1;
   };

// Post requests on the master to receive worker counts.
   MPI_Request* req_worker_count;
   if (is_master) {
      workers_per_node = new int[n_nodes];
      req_worker_count = new MPI_Request[n_nodes];
      for (proc = 0; proc < n_nodes; proc++) {
         MPI_Irecv(workers_per_node + proc, 1, MPI_INT, (is_server ? proc : proc + 1), MPI_ANY_TAG, server_comm, req_worker_count + proc);
      };
   };

// Each server sends its number of workers to the master (including master to self if allowed)
   if (is_server) MPI_Send(&workers_in_node, 1, MPI_INT, 0, 0, server_comm);

// Wait for worker count to be received, count the total
   if (is_master) {
      MPI_Waitall(n_nodes, req_worker_count, MPI_STATUSES_IGNORE);
      delete[] req_worker_count;
      n_workers = 0;
      for (proc = 0; proc < n_nodes; proc++) n_workers += workers_per_node[proc];
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
   if (is_master) delete[] workers_per_node;
   if (work_comm != MPI_COMM_NULL) MPI_Comm_free(&work_comm);
   if (server_comm != MPI_COMM_NULL) MPI_Comm_free(&server_comm);
   if (node_comm != MPI_COMM_NULL) MPI_Comm_free(&node_comm);
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
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 12/02/2024
*/
void TestMPIConfig(void)
{
// Each process prints its config to a different file
   std::string filename = "MPI_CONFIG/mpi_config_" + std::to_string(MPI_Config::glob_comm_rank);
   std::ofstream infofile;
   infofile.open(filename.c_str(), std::ofstream::out);

   infofile << "--------------------------------------------------------------------------------\n";
   infofile << "Size of global communicator: " << MPI_Config::glob_comm_size << std::endl;
   infofile << "Rank in global communicator: " << MPI_Config::glob_comm_rank << std::endl;
   infofile << "Master status: " << (MPI_Config::is_master ? "YES" : "NO") << std::endl;
   infofile << std::endl;

   infofile << "Number of nodes: " << MPI_Config::n_nodes << std::endl;
   infofile << "Size of server communicator: " << MPI_Config::server_comm_size << std::endl;
   infofile << "Rank in server communicator: " << MPI_Config::server_comm_rank << std::endl;
   infofile << "Server status: " << (MPI_Config::is_server ? "YES" : "NO") << std::endl;
   infofile << std::endl;

   infofile << "This is node: " << MPI_Config::my_node << std::endl;
   infofile << "Size of node communicator: " << MPI_Config::node_comm_size << std::endl;
   infofile << "Rank in node communicator: " << MPI_Config::node_comm_rank << std::endl;
   infofile << "Number of sockets per node: " << MPI_Config::n_sockets_per_node << std::endl;
   infofile << "Number of servers per socket: " << MPI_Config::n_servers_per_socket << std::endl;
   infofile << std::endl;

   infofile << "Total number of workers: " << MPI_Config::n_workers << std::endl;
   infofile << "Number of workers in this node: " << MPI_Config::workers_in_node << std::endl;
   infofile << "Size of worker communicator: " << MPI_Config::work_comm_size << std::endl;
   infofile << "Rank in worker communicator: " << MPI_Config::work_comm_rank << std::endl;
   infofile << "Worker status: " << (MPI_Config::is_worker ? "YES" : "NO") << std::endl;
   infofile << std::endl;

   if (MPI_Config::is_master) {
      for (auto node = 0; node < MPI_Config::n_nodes; node++) {
         infofile << "Number of workers in node " << node << ": " << MPI_Config::workers_per_node[node] << std::endl;
      };
   };
   infofile << "--------------------------------------------------------------------------------\n";
   
   infofile.close();
};

#endif

};

#endif
