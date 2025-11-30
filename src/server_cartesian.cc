/*!
\file server_cartesian.cc
\brief Implements a class of a data server for a uniform Cartesian grid
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_cartesian.hh"
#include "server_types.hh"
#include "reader_cartesian.hh"
#include "common/print_warn.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesian public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] file_name_pattern_in A string describing the file naming pattern
*/
template <typename HConfig>
ServerCartesian<HConfig>::ServerCartesian(const std::string& file_name_pattern_in)
                   : ServerInterface(file_name_pattern_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/01/2023
*/
template <typename HConfig>
void ServerCartesian<HConfig>::ReadData(const std::string data_file)
{
   ReadCartesianHeader(data_file.c_str(), line_width, 1);
   ReadCartesianData(data_file.c_str(), line_width, 0, 0);
   ReadCartesianGetDomain(domain_min.Data(), domain_max.Data());
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void ServerCartesian<HConfig>::ServerStart(void)
{
   ServerInterface::ServerInterfaceStart();
   block_served = new Block();

// Initialize the Cartesian library and read the data into memory
   std::string data_file = file_name_pattern + ".out";

   ReadData(data_file);
   domain_min *= unit_length_server / unit_length_fluid;
   domain_max *= unit_length_server / unit_length_fluid;

   MPI_Bcast(domain_min.Data(), 3, MPI_DOUBLE, 0, MPI::node_comm);
   MPI_Bcast(domain_max.Data(), 3, MPI_DOUBLE, 0, MPI::node_comm);
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
*/
template <typename HConfig>
void ServerCartesian<HConfig>::CleanReader(void)
{
   ReadCartesianClean();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void ServerCartesian<HConfig>::ServerFinish(void)
{
   CleanReader();
   delete block_served;

// This is a repeat of the function from ServerCartesian
   MPI_Type_free(&MPIBlockType);
   MPI_Type_free(&MPIStencilType);

   ServerInterface::ServerInterfaceFinish();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\return Number of clients that completed their tasks during this cycle

/Note "needvars" requests are only handled when SERVER_INTERP_ORDER = -1
*/
template <typename HConfig>
int ServerCartesian<HConfig>::ServerFunctions(void)
{
   if constexpr (server_interp_order == -1) {
// Handle "needvars" requests
      HandleNeedVarsRequests();
   }
// Handle "needblock" requests
   HandleNeedBlockRequests();
// Handle "stopserve" requests
   return HandleStopServeRequests();
};

/*!
\author Juan G Alonso Guzman
\date 08/03/2023
\param[in] pos    position array in reader coordinates
\param[out] vars  1D array containing all variables at pos
\param[out] found flag indicating whether variables where found or not for given position
*/
template <typename HConfig>
void ServerCartesian<HConfig>::GetBlockData(const double* pos, double* vars, int* found)
{
   ReadCartesianGetBlockData(pos, vars, found);
};

/*!
\author Juan G Alonso Guzman
\date 08/03/2023
*/
template <typename HConfig>
void ServerCartesian<HConfig>::HandleNeedVarsRequests(void)
{
   GeoVector pos_cart;
   DataFields datafields;
   int found, cpu, cpu_idx, count_needvars = 0;

// Service the "needvars" requests
   MPI_Testsome(MPI::node_comm_size, req_needvars, &count_needvars, index_needvars, MPI_STATUSES_IGNORE);

   for (cpu_idx = 0; cpu_idx < count_needvars; cpu_idx++) {
      cpu = index_needvars[cpu_idx];

// Obtain the variables requested
      pos_cart = buf_needvars[cpu].pos / unit_length_server * unit_length_fluid;
      GetBlockData(pos_cart.Data(), datafields.Array(), &found);
      if (!found) std::cerr << "Position not found\n";

// Send the variables to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(datafields.Array(), DataFields::size(), MPI_DOUBLE, cpu, MPI::tag::sendvars, MPI::node_comm);

// Post the receive for the next variables request from this worker.
      MPI_Irecv(&buf_needvars[cpu], 1, MPIInquiryType, cpu, MPI::tag::needvars, MPI::node_comm, &req_needvars[cpu]);
   };
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
\param[in] pos    position array in reader coordinates
\param[out] node  ID tag of block containing pos within its physical boundaries
*/
template <typename HConfig>
void ServerCartesian<HConfig>::GetBlock(const double* pos, int* node)
{
   ReadCartesianGetNode(pos, node);
};

/*!
\author Juan G Alonso Guzman
\date 12/01/2023
*/
template <typename HConfig>
void ServerCartesian<HConfig>::HandleNeedBlockRequests(void)
{
   GeoVector pos_cart;
   int cpu, cpu_idx, count_needblock = 0;

   // Service the "needblock" requests
   MPI_Testsome(MPI::node_comm_size, req_needblock, &count_needblock, index_needblock, MPI_STATUSES_IGNORE);

// Load the block requested. If the requestor does not know the node, figure it out.
   for (cpu_idx = 0; cpu_idx < count_needblock; cpu_idx++) {
      cpu = index_needblock[cpu_idx];

      if (buf_needblock[cpu].type) {
         pos_cart = buf_needblock[cpu].pos / unit_length_server * unit_length_fluid;
         GetBlock(pos_cart.Data(), &buf_needblock[cpu].node);
         if (buf_needblock[cpu].node == -1) throw ExServerError();
      };

      block_served->SetNode(buf_needblock[cpu].node);
      ServerInterface::LoadFromReader(block_served);
      block_served->LoadDimensions(unit_length_server);

// Send the block to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(block_served, 1, MPIBlockType, cpu, MPI::tag::sendblock, MPI::node_comm);
      if constexpr (server_interp_order > -1) {
         ServerInterface::LoadFieldsFromReader(block_served);
         MPI_Send(block_served->GetVariablesAddress(), block_served->GetVariableCount() * block_served->GetZoneCount(),
                  MPI_DOUBLE, cpu, MPI::tag::sendblock, MPI::node_comm);
      }
      if constexpr (server_interp_order > 0 && num_ghost_cells == 0) {
         ServerInterface::LoadNeighborsFromReader(block_served);
         MPI_Send(block_served->GetNeighborNodesAddress(), block_served->GetNeighborCount(),
                  MPI_INT, cpu, MPI::tag::sendblock, MPI::node_comm);
         MPI_Send(block_served->GetNeighborLevelsAddress(), block_served->GetNeighborLevelCount(),
                  MPI_INT, cpu, MPI::tag::sendblock, MPI::node_comm);
      }

// Post the receive for the next block request from this worker
      MPI_Irecv(&buf_needblock[cpu], 1, MPIInquiryType, cpu, MPI::tag::needblock, MPI::node_comm, &req_needblock[cpu]);
   };
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\return Number of clients that completed their tasks during this cycle
*/
template <typename HConfig>
int ServerCartesian<HConfig>::HandleStopServeRequests(void)
{
   int cpu, cpu_idx, count_stopserve = 0;

   // Service the "stopserve" requests. We assume that each worker sends a single request at the end of the simulation. 
   MPI_Testsome(MPI::node_comm_size, req_stopserve, &count_stopserve, index_stopserve, MPI_STATUSES_IGNORE);
   
// Cancel all "needblock" and "needstencil" receive requests from the cpus that have finished.
   for (cpu_idx = 0; cpu_idx < count_stopserve; cpu_idx++) {
      cpu = index_stopserve[cpu_idx];
      MPI_Cancel(&req_needblock[cpu]);
      MPI_Cancel(&req_needstencil[cpu]);
      MPI_Cancel(&req_needvars[cpu]);
      MPI_Request_free(&req_needblock[cpu]);
      MPI_Request_free(&req_needstencil[cpu]);
      MPI_Request_free(&req_needvars[cpu]);
   };

   return count_stopserve;
};

};
