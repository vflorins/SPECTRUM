/*!
\file mpi_config.hh
\brief Declares classes responsible for creating MPI various communicators and determining hardware configuration, as well as aiding in non-blocking communication
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MPI_CONFIG_HH
#define SPECTRUM_MPI_CONFIG_HH

#include "config.h"

#ifdef USE_MPI

#include <mpi.h>


//! Do not compile if server is needed in a serial run. The autotools configure process should not allow this but config.h could still be edited manually.
#if (EXEC_TYPE == EXEC_SERIAL) && (SERVER_TYPE != SERVER_SELF)
#error "SERIAL execution can only used with a SELF server type"
#endif

//! Master/Boss/Worker policy - set by configure
#if EXEC_TYPE == EXEC_SERIAL
#define ALLOW_MASTER_BOSS
#endif

//! Number of bosses per physical node
#define N_BOSSES_PER_NODE 1

#if SERVER_TYPE == SERVER_SELF
#define ALLOW_BOSS_WORKER
#else
#define NEED_SERVER
#endif

namespace Spectrum {

//! MPI tag for "cpu_avail" message (W->M)
const int tag_cpuavail = 1011;

//! MPI tag for "distribution" message (W->M)
const int tag_distrdata = 1012;

//! MPI tag for "need more" message (M->W)
const int tag_needmore_MW = 1013;

//! MPI tag for "position" message (W->B)
const int tag_positdata = 1014;

//! MPI tag for "auxiliary" message (B->W)
const int tag_fielddata = 1015;

//! MPI tag for "need more" message (W->B)
const int tag_needmore_WB = 1016;

/*!
\brief Class to store the information about the topology of the MPI environment
\author Vladimir Florinski

This class is used to partition the simulation into nodes (shared memory) and assigning bosses for each node. It sets up three communicators, between the processes in a single node, between the boss processes, and between the worker processes. It also identifies how many worker processes are available in the simulation. An application should create a single instance of this type and pass a pointer to other objects. As a convenience feature, the RNG is initialized using the global rank and timer.
*/
struct MPI_Config
{
//! Global communicator
   inline static MPI_Comm glob_comm;

//! Intra-node comminicator (all processes on this shared memory node)
   inline static MPI_Comm node_comm;

//! Boss communicator
   inline static MPI_Comm boss_comm;

//! Worker communicator
   inline static MPI_Comm work_comm;

//! Number of processes in the global communicator
   inline static int glob_comm_size;

//! Number of processes in the node communicator - not always equal to number of processes on the node
   inline static int node_comm_size;

//! Number of processes in the boss communicator - not always equal to the number of nodes
   inline static int boss_comm_size;

//! Number of processes in the worker communicator
   inline static int work_comm_size;

//! Rank in the global communicator
   inline static int glob_comm_rank;

//! Rank in the node communicator
   inline static int node_comm_rank;

//! Rank in the boss communicator
   inline static int boss_comm_rank;

//! Rank in the worker communicator
   inline static int work_comm_rank;

//! Number of nodes
   inline static int n_nodes;

//! Number of sockets per node
   inline static int n_sockets_per_node;

//! Number of bosses per socket
   inline static int n_bosses_per_socket;

//! This process's node
   inline static int my_node;

//! Whether this process is a master
   inline static bool is_master;

//! Whether this process is a boss
   inline static bool is_boss;

//! Whether this process is a worker
   inline static bool is_worker;

//! Number of worker processes in this node
   inline static int workers_in_node;

//! Total number of workers
   inline static int n_workers;

//! Number of worker processes on each node (master only)
   inline static int* workers_per_node;

//! Default constructor
   MPI_Config(void) = delete;

//! Copy constructor - deleted, no copy allowed
   MPI_Config(const MPI_Config& other) = delete;

//! Constructor with arguments
   MPI_Config(int argc, char** argv);

//! Destructor
   ~MPI_Config();
};

/*!
\brief Class to handle non-blocking receives
\author Juan G Alonso Guzman

This class houses all the variables necessary to implement non-blocking communications
*/
struct MPI_Request_Info
{
//! MPI Request array
   MPI_Request* mpi_req = nullptr;

//! CPU rank array of processes that have completed requests
   int* cpu_rank = nullptr;

//! Count of CPUs that have completed requests
   int count = 0;

//! Constructor with a parameter
   MPI_Request_Info(int size);

//! Destructor
   ~MPI_Request_Info();
};

#ifdef GEO_DEBUG
//! Self-test for the class
void TestMPIConfig(void);
#endif

};

#endif

#endif
