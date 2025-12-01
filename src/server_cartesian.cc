/*!
\file server_cartesian.cc
\brief Implements a class of a data server for a uniform Cartesian grid
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_cartesian.hh"
#include "server_types.hh"
#include "common/print_warn.hh"
#include "common/physics.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesian methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Juan G Alonso Guzman
\date 08/01/2023
*/
template <typename HConfig>
void ServerCartesian<HConfig>::ReadData(const std::string data_file)
{
   reader.ReadCartesianHeader(data_file.c_str(), line_width, 1);
   reader.ReadCartesianData(data_file.c_str(), line_width, 0, 0);
   reader.ReadCartesianGetDomain(domain_min.Data(), domain_max.Data());
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void ServerCartesian<HConfig>::ServerStart(void)
{
   ServerBase::ServerStart();
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
   reader.ReadCartesianClean();
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

// todo review, is this done by ServerBase or ServerInterface?
   MPI_Type_free(&MPIBlockType);
   MPI_Type_free(&MPIStencilType);

   ServerBase::ServerFinish();
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\return Number of clients that completed their tasks during this cycle

\note "needvars" requests are only handled when SERVER_INTERP_ORDER = -1
*/
template <typename HConfig>
int ServerCartesian<HConfig>::ServerFunctions(void)
{
   if constexpr (server_interp_order == -1) {
// Handle "needvars" requests
      ServerBase::HandleNeedVarsRequests();
   }
// Handle "needblock" requests
   ServerBase::HandleNeedBlockRequests();
// Handle "stopserve" requests
   return ServerBase::HandleStopServeRequests();
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Request the bounding dimensions
*/
template <typename HConfig>
void ServerCartesian<HConfig>::LoadFromReader(BlockPtr& blockptr)
{
   reader.ReadCartesianGetBlockCorners(blockptr->node, blockptr->face_min.Data(), blockptr->face_max.Data());
}


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Request all neighbors
*/
template <typename HConfig>
void ServerCartesian<HConfig>::LoadNeighborsFromReader(BlockPtr& blockptr)
{
   reader.ReadCartesianGetNodeNeighbors(blockptr->node, blockptr->neighbor_nodes, blockptr->neighbor_levels);
}


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Load all fields into the block
*/
template <typename HConfig>
void ServerCartesian<HConfig>::LoadFieldsFromReader(BlockPtr &blockptr)
{
   reader.ReadCartesianGetBlockData(blockptr->node, blockptr->fields[0].Array());
}



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
   reader.ReadCartesianGetBlockData(pos, vars, found);
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
   reader.ReadCartesianGetNode(pos, node);
};


/*!
\author Lucius Schoenbaum
\date 11/28/2025
\param[in] pos    position array in reader coordinates
\param[out] stencil  stencil to be populated
*/
template <typename HConfig>
void ServerCartesian<HConfig>::GetStencil(const double* pos, Stencil* stencil) {
// unused because HandleNeedStencilRequests is not called.
   return;
}


};
