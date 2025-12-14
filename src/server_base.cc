/*!
\file server_base.cc
\brief Implements a base class of a data server from an external source
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_base.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {


//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Vladimir Florinski
\date 10/26/2022
*/
template <typename HConfig>
void ServerBase<HConfig>::ServerStart(void)
{
   ServerInterface::ServerInterfaceStart();

// Indices of processes returned by "Testsome"
   index_needblock   = new int[MPI::node_comm_size];
   index_needstencil = new int[MPI::node_comm_size];
   index_needvars    = new int[MPI::node_comm_size];
   index_stopserve   = new int[MPI::node_comm_size];

// Request arrays
   req_needblock   = new MPI_Request[MPI::node_comm_size];
   req_needstencil = new MPI_Request[MPI::node_comm_size];
   req_needvars    = new MPI_Request[MPI::node_comm_size];
   req_stopserve   = new MPI_Request[MPI::node_comm_size];

// Message buffers
   buf_needblock   = new Inquiry[MPI::node_comm_size];
   buf_needstencil = new Inquiry[MPI::node_comm_size];
   buf_needvars    = new Inquiry[MPI::node_comm_size];

// Post initial receives for all request types
   req_needblock[0] = MPI_REQUEST_NULL;
   req_needstencil[0] = MPI_REQUEST_NULL;
   req_needvars[0] = MPI_REQUEST_NULL;
   req_stopserve[0] = MPI_REQUEST_NULL;
   for (int cpu = 1; cpu < MPI::node_comm_size; cpu++) {
      MPI_Irecv(&buf_needblock[cpu], 1, MPIInquiryType, cpu, MPI::tag::needblock, MPI::node_comm, &req_needblock[cpu]);
      MPI_Irecv(&buf_needstencil[cpu], 1, MPIInquiryType, cpu, MPI::tag::needstencil, MPI::node_comm, &req_needstencil[cpu]);
      MPI_Irecv(&buf_needvars[cpu], 1, MPIInquiryType, cpu, MPI::tag::needvars, MPI::node_comm, &req_needvars[cpu]);
      MPI_Irecv(nullptr, 0, MPI_BYTE, cpu, MPI::tag::stopserve, MPI::node_comm, &req_stopserve[cpu]);
   };
};

/*!
\author Vladimir Florinski
\date 10/26/2022
*/
template <typename HConfig>
void ServerBase<HConfig>::ServerFinish(void)
{
// Deallocate index arrays
   delete[] index_needblock;
   delete[] index_needstencil;
   delete[] index_needvars;
   delete[] index_stopserve;

// Deallocate request arrays
   delete[] req_needblock;
   delete[] req_needstencil;
   delete[] req_needvars;
   delete[] req_stopserve;
   
// Deallocate buffers
   delete[] buf_needblock;
   delete[] buf_needstencil;
   delete[] buf_needvars;

   ServerInterface::ServerInterfaceFinish();
};



/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/28/2025
*/
template <typename HConfig>
void ServerBase<HConfig>::HandleNeedVarsRequests(void)
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
\author Lucius Schoenbaum
\date 11/27/2025
\return Number of clients that completed their tasks during this cycle
*/
template <typename HConfig>
int ServerBase<HConfig>::HandleStopServeRequests(void)
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


/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/21/2025
*/
template <typename HConfig>
void ServerBase<HConfig>::HandleNeedBlockRequests(void)
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
      LoadFromReader(block_served);
      block_served->LoadDimensions(unit_length_server);

// Send the block to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(block_served, 1, MPIBlockType, cpu, MPI::tag::sendblock, MPI::node_comm);
      if constexpr (server_interpolation_order > -1) {
         LoadFieldsFromReader(block_served);
         MPI_Send(block_served->GetVariablesAddress(), block_served->GetVariableCount() * block_served->GetZoneCount(),
                  MPI_DOUBLE, cpu, MPI::tag::sendblock, MPI::node_comm);
      }
      if constexpr (server_interpolation_order > 0 && num_ghost_cells == 0) {
         LoadNeighborsFromReader(block_served);
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
*/
template <typename HConfig>
void ServerBase<HConfig>::HandleNeedStencilRequests(void)
{
   GeoVector pos_batl;
   int cpu, cpu_idx, count_needstencil = 0;

// Service the "needstencil" requests
   MPI_Testsome(MPI::node_comm_size, req_needstencil, &count_needstencil, index_needstencil, MPI_STATUSES_IGNORE);

   for (cpu_idx = 0; cpu_idx < count_needstencil; cpu_idx++) {
      cpu = index_needstencil[cpu_idx];

// Obtain the stencil requested.
      pos_batl = buf_needstencil[cpu].pos / unit_length_server * unit_length_fluid;
      GetStencil(pos_batl.Data(), &stencil);
//         spectrum_get_interpolation_stencil(pos_batl.Data(), &stencil.n_elements, stencil.blocks, &stencil.zones[0][0], stencil.weights);

// Send the stencil to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(&stencil, 1, MPIStencilType, cpu, MPI::tag::sendstencil, MPI::node_comm);

// Post the receive for the next stencil request from this worker.
      MPI_Irecv(&buf_needstencil[cpu], 1, MPIInquiryType, cpu, MPI::tag::needstencil, MPI::node_comm, &req_needstencil[cpu]);
   };
};




};
