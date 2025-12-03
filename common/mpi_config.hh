/*!
\file mpi_config.hh
\brief Declares classes responsible for creating MPI various communicators and determining hardware configuration, as well as aiding in non-blocking communication
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MPI_CONFIG_HH
#define SPECTRUM_MPI_CONFIG_HH

#include "compiletime_lists.hh"

#ifdef USE_MPI
#include <mpi.h>
#endif


namespace Spectrum {


template <typename HConfig, bool MPI_enabled>
struct MPI;


/*!
\brief Class to store the information about the topology of the MPI environment
\author Vladimir Florinski
\author Lucius Schoenbaum
Instantiation when MPI is not enabled.
*/
template <typename HConfig_>
struct MPI<HConfig_, false> {
//! Number of processes in the global communicator
   static constexpr int glob_comm_size = 1;
//! Number of processes in the node communicator - not always equal to number of processes on the node
   static constexpr int node_comm_size = 1;
//! Number of processes in the server communicator - not always equal to the number of nodes
   static constexpr int server_comm_size = 1;
//! Number of processes in the worker communicator
   static constexpr int work_comm_size = 1;
//! Rank in the global communicator
   static constexpr int glob_comm_rank = 0;
//! Rank in the node communicator
   static constexpr int node_comm_rank = 0;
//! Rank in the server communicator
   static constexpr int server_comm_rank = 0;
//! Rank in the worker communicator
   static constexpr int work_comm_rank = 0;
//! Number of nodes
   static constexpr int n_nodes = 1;
//! Number of sockets per node
   static constexpr int n_sockets_per_node = 0;
//! Number of servers per socket
   static constexpr int n_servers_per_socket = 1;
//! This process's node
   static constexpr int my_node = 0;
//! Whether this process is a master
   static constexpr bool is_master = true;
//! Whether this process is a server
   static constexpr bool is_server = true;
//! Whether this process is a worker
   static constexpr bool is_worker = true;
//! Number of worker processes in this node
   static constexpr int workers_in_node = 1;
//! Total number of workers
   static constexpr int n_workers = 1;
//! Number of worker processes on each node (master only)
   static constexpr int workers_per_node[1] = {1};

   static constexpr bool is_parallel() {
      return false;
   }
};


/*!
\brief Class to store the information about the topology of the MPI environment
\author Vladimir Florinski
\author Lucius Schoenbaum

This class is used to partition the simulation into nodes (shared memory) and assigning servers for each node. It sets up three communicators, between the processes in a single node, between the server processes, and between the worker processes. It also identifies how many worker processes are available in the simulation. An application should create a single instance of this type and pass a pointer to other objects. As a convenience feature, the RNG is initialized using the global rank and timer.
*/
template <typename HConfig_>
struct MPI<HConfig_, true>
{

   using HConfig = HConfig_;
   using Config = HConfig::SimulationConfig;

   enum tag {
//! MPI tag for "need block" message (W->S)
      endblock = 1001,
//! MPI tag for "send block" message (S->W)
      sendblock = 1002,
//! MPI tag for "need stencil" message (W->S)
      needstencil = 1003,
//! MPI tag for "send stencil" message (S->W)
      sendstencil = 1004,
//! MPI tag for "need vars" message (W->S)
      needvars = 1005,
//! MPI tag for "send vars" message (S->W)
      sendvars = 1006,
//! MPI tag for "stop serve" message (W->S)
      stopserve = 1007,
//! MPI tag for "cpu_avail" message (W->M)
      cpuavail = 1011,
//! MPI tag for "distribution" message (W->M)
      distrdata = 1012,
//! MPI tag for "need more" message (M->W)
      needmore_MW = 1013,
//! MPI tag for "position" message (W->S)
      positdata = 1014,
//! MPI tag for "auxiliary" message (S->W)
      fielddata = 1015,
//! MPI tag for "need more" message (W->S)
      needmore_WS = 1016,
   };

//! Global communicator
   inline static MPI_Comm glob_comm;

//! Intra-node comminicator (all processes on this shared memory node)
   inline static MPI_Comm node_comm;

//! Server communicator
   inline static MPI_Comm server_comm;

//! Worker communicator
   inline static MPI_Comm work_comm;

//! Number of processes in the global communicator
   inline static int glob_comm_size;

//! Number of processes in the node communicator - not always equal to number of processes on the node
   inline static int node_comm_size;

//! Number of processes in the server communicator - not always equal to the number of nodes
   inline static int server_comm_size;

//! Number of processes in the worker communicator
   inline static int work_comm_size;

//! Rank in the global communicator
   inline static int glob_comm_rank;

//! Rank in the node communicator
   inline static int node_comm_rank;

//! Rank in the server communicator
   inline static int server_comm_rank;

//! Rank in the worker communicator
   inline static int work_comm_rank;

//! Number of nodes
   inline static int n_nodes;

//! Number of sockets per node
   inline static int n_sockets_per_node;

//! Number of servers per socket
   inline static int n_servers_per_socket;

//! This process's node
   inline static int my_node;

//! Whether this process is a master
   inline static bool is_master;

//! Whether this process is a server
   inline static bool is_server;

//! Whether this process is a worker
   inline static bool is_worker;

//! Number of worker processes in this node
   inline static int workers_in_node;

//! Total number of workers
   inline static int n_workers;

//! Number of worker processes on each node (master only)
   inline static int* workers_per_node;

//! Default constructor
   MPI(void) = delete;

//! Copy constructor - deleted, no copy allowed
   MPI(const MPI& other) = delete;

//! Constructor with arguments
   MPI(int argc, char** argv);

//! Destructor
   ~MPI();

   //! Self-test for the class
   void TestMPIConfig(void) const requires (HConfig::buildmode == BuildMode::debug);

//! Whether parallel, defined to mean that there is more than one worker process, or else there is a server that is not a worker.
   static bool is_parallel() {
      if constexpr (HConfig::data_background()) {
         return (work_comm_size > 1) && !HConfig::BackgroundConfig::allow_server_worker;
      }
      else {
         return work_comm_size > 1;
      }
   }

};

};

// Something like this is needed for templated classes
#include "mpi_config.cc"

#endif
