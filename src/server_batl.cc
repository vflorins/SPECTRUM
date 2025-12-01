/*!
\file server_batl.cc
\brief Implements a class of a data server from BATL
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_batl.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATL methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/19/2025
*/
template <typename HConfig>
void ServerBATL<HConfig>::ServerStart(void)
{
   ServerBase::ServerStart();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/19/2025
*/
template <typename HConfig>
void ServerBATL<HConfig>::ServerFinish(void)
{
   ServerBase::ServerFinish();
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
\return Number of clients that completed their tasks during this cycle

\note "needvars" requests are only handled when SERVER_INTERP_ORDER = -1. "needstencil" requests are only handled when using 1st order interpolation and either there are no ghost cells or stencils are always requested from BATL (for debugging purposes)
*/
template <typename HConfig>
int ServerBATL<HConfig>::ServerFunctions(void)
{
   if constexpr (server_interp_order == -1) {
// Handle "needvars" requests
      ServerBase::HandleNeedVarsRequests();
   }
   else if constexpr (server_interp_order == 1 && (num_ghost_cells == 0 || request_stencil_from_batl)) {
// Handle "needstencil" requests
      ServerBase::HandleNeedStencilRequests();
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
void ServerBATL<HConfig>::LoadFromReader(BlockPtr &blockptr)
{
   spectrum_get_block_corners(blockptr->node, blockptr->face_min.Data(), blockptr->face_max.Data());
}


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Request all neighbors
*/
template <typename HConfig>
void ServerBATL<HConfig>::LoadNeighborsFromReader(BlockPtr& blockptr)
{
   spectrum_get_all_neighbor_copies(blockptr->node, blockptr->neighbor_nodes, blockptr->neighbor_levels);
}


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Load all fields into the block
*/
template <typename HConfig>
void ServerBATL<HConfig>::LoadFieldsFromReader(BlockPtr& blockptr)
{
   spectrum_get_block_data(blockptr->node, blockptr->fields[0].Array());
}


/*!
\author Juan G Alonso Guzman
\date 12/01/2023
*/
template <typename HConfig>
void ServerBATL<HConfig>::ReadData(const std::string data_file)
{
// Initialize the BATL library and read the data into memory
   MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_SELF);
   spectrum_init_mpi(fcomm);

   wrapamr_read_header(data_file.c_str(), line_width, 1);
   int vars_inds_fortran[DataFields::size()];
   for (int i = 0; i < DataFields::size(); i++) vars_inds_fortran[i] = i+1;
   wrapamr_read_file_partial(data_file.c_str(), line_width, 0, 0, DataFields::size(), vars_inds_fortran);
   wrapamr_get_domain(domain_min.Data(), domain_max.Data());
};


/*!
\author Juan G Alonso Guzman
\date 07/28/2023
*/
template <typename HConfig>
void ServerBATL<HConfig>::CleanReader(void)
{
   wrapamr_clean();
};

/*!
\author Juan G Alonso Guzman
\date 08/03/2023
\param[in] pos    position array in reader coordinates
\param[out] vars  1D array containing all variables at pos
\param[out] found flag indicating whether variables where found or not for given position
*/
template <typename HConfig>
void ServerBATL<HConfig>::GetBlockData(const double* pos, double* vars, int* found)
{
   wrapamr_get_data_serial(pos, vars, found);
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
\param[in] pos    position array in reader coordinates
\param[out] node  ID tag of block containing pos within its physical boundaries
*/
template <typename HConfig>
void ServerBATL<HConfig>::GetBlock(const double* pos, int* node)
{
   spectrum_get_node(pos, node);
};


/*!
\author Lucius Schoenbaum
\date 11/28/2025
\param[in] pos    position array in reader coordinates
\param[out] stencil  stencil to be populated
*/
template <typename HConfig>
void ServerBATL<HConfig>::GetStencil(const double* pos, Stencil* stencil) {
   spectrum_get_interpolation_stencil(pos, &stencil->n_elements, stencil->blocks, &stencil->zones[0][0], stencil->weights);
}


};
