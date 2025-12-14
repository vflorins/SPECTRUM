/*!
\file background_data_base.hh
\brief Implements a base background class using grid data on distributed memory
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/


#include "background_data_base.hh"

namespace Spectrum {



/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/24/2025
*/
template <typename HConfig>
void BackgroundDataBase<HConfig>::Start(void)
{

// if a worker is a server, the interface will be started during its server initialization.
   if constexpr (!allow_server_worker) {
      ServerInterface::ServerInterfaceStart();
   }

   cache_line.Empty();
   stencil_outcomes[0] = stencil_outcomes[1] = stencil_outcomes[2] = 0;
   num_blocks_requested = 0;

   MPI_Bcast(domain_min.Data(), 3, MPI_DOUBLE, 0, MPI::node_comm);
   MPI_Bcast(domain_max.Data(), 3, MPI_DOUBLE, 0, MPI::node_comm);
// Prime "block_pri" and "block_sec" with stub blocks that always fail tests.
// These and "block_stn" must be smart pointers to avoid double free or corruption errors.
   block_pri = std::make_shared<Block>();
   block_sec = std::make_shared<Block>();
   block_stn = std::make_shared<Block>();
   block_pri->SetDimensions(domain_max, domain_min);
   block_sec->SetDimensions(domain_max, domain_min);
// todo I do not think necessary - if necessary, complications arise.
//   LoadFromReader(block_pri);
//   LoadFromReader(block_sec);
   block_pri->LoadDimensions(1.0);
   block_sec->LoadDimensions(1.0);
};




/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/22/2025
*/
template <typename HConfig>
void BackgroundDataBase<HConfig>::Finish(void)
{
   if constexpr (HConfig::buildmode == BuildMode::debug) {
      PrintStencilOutcomes();
   }
   MPI_Send(nullptr, 0, MPI_BYTE, 0, MPI::tag::stopserve, MPI::node_comm);
// if a worker is a server, the interface will be started during its server initialization.
   if constexpr (!allow_server_worker) {
      ServerInterface::ServerInterfaceFinish();
   }
};







/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 12/01/2023
\return Index of the block in "cache_line"
*/
template <typename HConfig>
int BackgroundDataBase<HConfig>::RequestBlock(void)
{
   int bidx;
   BlockPtr block_new;

// Test whether the block is cached. Either call will renew the block if it is present.
   if (_inquiry.type) {
      if (block_pri->PositionInside(_inquiry.pos)) bidx = block_pri->GetNode();
      else if (block_sec->PositionInside(_inquiry.pos)) bidx = block_sec->GetNode();
      else bidx = cache_line.PosOwner(_inquiry.pos);
   }
   else {
      if (block_pri->GetNode() == _inquiry.node) bidx = block_pri->GetNode();
      else if (block_sec->GetNode() == _inquiry.node) bidx = block_sec->GetNode();
      else bidx = cache_line.Present(_inquiry.node);
   };

// Block is not in the cache, request it from the server.
   if (bidx == -1) {
      MPI_Send(&_inquiry, 1, MPIInquiryType, 0, MPI::tag::needblock, MPI::node_comm);
      num_blocks_requested++;

// Allocate memory for block.
      block_new = std::make_shared<Block>();

// Receive the block in 4 parts (member data plus 3 dynamic arrays).
// This is called even if SERVER_INTERP_ORDER is -1 to import the block dimensions
      MPI_Recv(block_new.get(), 1, MPIBlockType, 0, MPI::tag::sendblock, MPI::node_comm, MPI_STATUS_IGNORE);

      if constexpr (server_interpolation_order > -1) {
         MPI_Recv(block_new->GetVariablesAddress(), block_new->GetVariableCount() * block_new->GetZoneCount(), MPI_DOUBLE, 0,
                  MPI::tag::sendblock, MPI::node_comm, MPI_STATUS_IGNORE);
      }
      if constexpr (server_interpolation_order > 0 && num_ghost_cells == 0) {
         MPI_Recv(block_new->GetNeighborNodesAddress(), block_new->GetNeighborCount(), MPI_INT, 0,
                  MPI::tag::sendblock, MPI::node_comm, MPI_STATUS_IGNORE);
         MPI_Recv(block_new->GetNeighborLevelsAddress(), block_new->GetNeighborLevelCount(), MPI_INT, 0,
                  MPI::tag::sendblock, MPI::node_comm, MPI_STATUS_IGNORE);
      }

// Insert the block into the cache
      block_new->ConfigureProperties();
      cache_line.AddBlock(block_new);
      bidx = block_new->GetNode();
   };

   return bidx;
};


/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\param[in] zone_lo   lower interpolation zone
\param[in] zone_hi   higher interpolation zone
\param[in] offset_lo lower interpolation offset
\param[in] offset_hi higher interpolation offset
*/
template <typename HConfig>
void BackgroundDataBase<HConfig>::InteriorInterpolationStencil(const MultiIndex zone_lo, const MultiIndex zone_hi,
                                                                    const GeoVector offset_lo, const GeoVector offset_hi,
                                                                    const GeoVector delta)
{
   stencil.zones[0] = zone_lo;
   stencil.weights[0] = offset_hi[0] * offset_hi[1] * offset_hi[2];
   stencil.derivatives[ 0] = -offset_hi[1] * offset_hi[2] / delta[0];
   stencil.derivatives[ 1] = -offset_hi[0] * offset_hi[2] / delta[1];
   stencil.derivatives[ 2] = -offset_hi[0] * offset_hi[1] / delta[2];

   stencil.zones[1] = stencil.zones[0];
   stencil.zones[1].i++;
   stencil.weights[1] = offset_lo[0] * offset_hi[1] * offset_hi[2];
   stencil.derivatives[ 3] =  offset_hi[1] * offset_hi[2] / delta[0];
   stencil.derivatives[ 4] = -offset_lo[0] * offset_hi[2] / delta[1];
   stencil.derivatives[ 5] = -offset_lo[0] * offset_hi[1] / delta[2];

   stencil.zones[2] = stencil.zones[0];
   stencil.zones[2].j++;
   stencil.weights[2] = offset_hi[0] * offset_lo[1] * offset_hi[2];
   stencil.derivatives[ 6] = -offset_lo[1] * offset_hi[2] / delta[0];
   stencil.derivatives[ 7] =  offset_hi[0] * offset_hi[2] / delta[1];
   stencil.derivatives[ 8] = -offset_hi[0] * offset_lo[1] / delta[2];

   stencil.zones[3] = stencil.zones[1];
   stencil.zones[3].j++;
   stencil.weights[3] = offset_lo[0] * offset_lo[1] * offset_hi[2];
   stencil.derivatives[ 9] =  offset_lo[1] * offset_hi[2] / delta[0];
   stencil.derivatives[10] =  offset_lo[0] * offset_hi[2] / delta[1];
   stencil.derivatives[11] = -offset_lo[0] * offset_lo[1] / delta[2];

   stencil.zones[4] = stencil.zones[0];
   stencil.zones[4].k++;
   stencil.weights[4] = offset_hi[0] * offset_hi[1] * offset_lo[2];
   stencil.derivatives[12] = -offset_hi[1] * offset_lo[2] / delta[0];
   stencil.derivatives[13] = -offset_hi[0] * offset_lo[2] / delta[1];
   stencil.derivatives[14] =  offset_hi[0] * offset_hi[1] / delta[2];

   stencil.zones[5] = stencil.zones[4];
   stencil.zones[5].i++;
   stencil.weights[5] = offset_lo[0] * offset_hi[1] * offset_lo[2];
   stencil.derivatives[15] =  offset_hi[1] * offset_lo[2] / delta[0];
   stencil.derivatives[16] = -offset_lo[0] * offset_lo[2] / delta[1];
   stencil.derivatives[17] =  offset_lo[0] * offset_hi[1] / delta[2];

   stencil.zones[6] = stencil.zones[4];
   stencil.zones[6].j++;
   stencil.weights[6] = offset_hi[0] * offset_lo[1] * offset_lo[2];
   stencil.derivatives[18] = -offset_lo[1] * offset_lo[2] / delta[0];
   stencil.derivatives[19] =  offset_hi[0] * offset_lo[2] / delta[1];
   stencil.derivatives[20] =  offset_hi[0] * offset_lo[1] / delta[2];

   stencil.zones[7] = zone_hi;
   stencil.weights[7] = offset_lo[0] * offset_lo[1] * offset_lo[2];
   stencil.derivatives[21] =  offset_lo[1] * offset_lo[2] / delta[0];
   stencil.derivatives[22] =  offset_lo[0] * offset_lo[2] / delta[1];
   stencil.derivatives[23] =  offset_lo[0] * offset_lo[1] / delta[2];

};





/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void BackgroundDataBase<HConfig>::PrintStencilOutcomes(void) const requires (HConfig::buildmode == BuildMode::debug)
{
   std::cerr << "Stencil outcomes: " << std::setw(10) << stencil_outcomes[0]
             << std::setw(10) << stencil_outcomes[1]
             << std::setw(10) << stencil_outcomes[2] << std::endl;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
*/
template <typename HConfig>
void BackgroundDataBase<HConfig>::PrintNumBlocksRequested(void) const requires (HConfig::buildmode == BuildMode::debug)
{
   std::cerr << "Number of blocks requested: " << std::setw(10) << num_blocks_requested << std::endl;
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/27/2025
\return Number of blocks in the cache
*/
template <typename HConfig>
int BackgroundDataBase<HConfig>::GetNCachedBlocks(void) const requires (HConfig::buildmode == BuildMode::debug)
{
   return cache_line.size();
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/27/2025
*/
template <typename HConfig>
void BackgroundDataBase<HConfig>::InvalidateCache(void) requires (HConfig::buildmode == BuildMode::debug)
{
   cache_line.Empty();
};





};


