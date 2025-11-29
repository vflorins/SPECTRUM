/*!
\file background_data_cartesian.cc
\brief Implements a background class using data from a uniform Cartesian grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_data_cartesian.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDataCartesian methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Request the bounding dimensions
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::LoadFromReader(BlockPtr &blockptr) const
{
   // todo review for header/include/linking issues
   if constexpr (HConfig::background == Config::Background::DataCartesian)
      ReadCartesianGetBlockCorners(blockptr->node, blockptr->face_min.Data(), blockptr->face_max.Data());
   else if constexpr (HConfig::background == Config::Background::DataBATL)
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
void BackgroundDataCartesian<HConfig>::LoadNeighborsFromReader(BlockPtr &blockptr) const
{
   if constexpr (HConfig::background == Config::Background::DataCartesian)
      ReadCartesianGetNodeNeighbors(blockptr->node, blockptr->neighbor_nodes, blockptr->neighbor_levels);
   else if constexpr (HConfig::background == Config::Background::DataBATL)
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
void BackgroundDataCartesian<HConfig>::LoadFieldsFromReader(BlockPtr &blockptr) const
{
   if constexpr (HConfig::background == Config::Background::DataCartesian)
      ReadCartesianGetBlockData(blockptr->node, blockptr->fields[0].Array());
   else if constexpr (HConfig::background == Config::Background::DataBATL)
      spectrum_get_block_data(blockptr->node, blockptr->fields[0].Array());
}





/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::Start(void) {

   ServerInterface::ServerInterfaceStart();

   BlockCache::Empty();
   stencil_outcomes[0] = stencil_outcomes[1] = stencil_outcomes[2] = 0;
   num_blocks_requested = 0;

   MPI_Bcast(domain_min.Data(), 3, MPI_DOUBLE, 0, MPI::node_comm);
   MPI_Bcast(domain_max.Data(), 3, MPI_DOUBLE, 0, MPI::node_comm);

// Prime "block_pri" and "block_sec" with stub blocks that always fail tests.
// These and "block_stn" must be smart pointers to avoid double free or corruption errors.
   MakeSharedBlock(block_pri);
   block_pri->SetDimensions(domain_max, domain_min);
   LoadFromReader(block_pri);
   block_pri->LoadDimensions(1.0);
   MakeSharedBlock(block_sec);
   block_sec->SetDimensions(domain_max, domain_min);
   LoadFromReader(block_sec);
   block_sec->LoadDimensions(1.0);
   MakeSharedBlock(block_stn);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::Finish(void)
{
   MPI_Send(nullptr, 0, MPI_BYTE, 0, MPI::tag::stopserve, MPI::node_comm);
   ServerInterface::ServerInterfaceFinish();
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\param[out] block_new pointer to block type
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::MakeSharedBlock(BlockPtr &block_new)
{
   block_new = std::make_shared<Block>();
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
void BackgroundDataCartesian<HConfig>::InteriorInterpolationStencil(const MultiIndex zone_lo, const MultiIndex zone_hi,
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

   stencil.n_elements = 8;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 12/01/2023
\param[in] pos Interpolation point position
\return Status: 0 if local interpolation was used, 1, 2, and 3 if plane interpolation was used, 4 if external interpolation was used (corner or edge)
*/
template <typename HConfig>
int BackgroundDataCartesian<HConfig>::BuildInterpolationStencil(const GeoVector& pos)
{
   int pri_idx, sec_idx, tmp_idx, xyz, iz, status, n_blocks = 0;
   MultiIndex block_size, node_idx, zone_lo, zone_hi;
   GeoVector offset_lo, offset_hi, delta;

// Get primary block data
   pri_idx = block_pri->GetNode();
   block_pri->GetZoneOffset(pos, zone_lo, offset_lo);
   delta = block_pri->GetZoneLength();

   offset_hi = 1.0 - offset_lo;
   zone_hi = zone_lo + 1;

// Default to local interpolation in the primary block
   InteriorInterpolationStencil(zone_lo, zone_hi, offset_lo, offset_hi, delta);

   if constexpr (num_ghost_cells == 0) {
      block_size = block_pri->GetBlockSize();
// Correct stencil zones and blocks for multi-block interpolation
// The weights do not change because of the uniformity of the grid
      for (iz = 0; iz < stencil.n_elements; iz++) {
         node_idx = mi_ones;

// Find which indices fall outside of the primary block's boundaries
         for (xyz = 0; xyz < 3; xyz++) {

// Use "previous" block
            if (stencil.zones[iz][xyz] == -1) {
               stencil.zones[iz][xyz] = block_size[xyz] - 1;
               node_idx[xyz]--;
               status = xyz + 1;
            }

// Use "next" block
            else if (stencil.zones[iz][xyz] == block_size[xyz]) {
               stencil.zones[iz][xyz] = 0;
               node_idx[xyz]++;
               status = xyz + 1;
            };
         };

// Assign the appropriate block for each zone
         tmp_idx = block_pri->GetNeighborNode(node_idx);
         if (tmp_idx == pri_idx) stencil.blocks[iz] = pri_idx;
         else {
            _inquiry.type = 0;
            _inquiry.node = tmp_idx;
            stencil.blocks[iz] = RequestBlock();
            sec_idx = tmp_idx;
            n_blocks++;
         };
      };
   }

   if (n_blocks == 0) status = 0;
   else if (n_blocks == 4) {
// If "block_pri" or "block_sec" is the position owner (based on the call to RequestBlock), we don't need to acccess the cache
      if (block_sec->GetNode() != sec_idx) {
         if (block_pri->GetNode() == sec_idx) block_sec = block_pri;
         else block_sec = cache_line[sec_idx];
      };
   }
   else status = 4;

   return status;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 12/01/2023
\return Index of the block in "cache_line"
*/
template <typename HConfig>
int BackgroundDataCartesian<HConfig>::RequestBlock(void)
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
      MakeSharedBlock(block_new);

// Adjust for ghost cells if necessary
// todo deprecated, happens at compile time
//#if SERVER_NUM_GHOST_CELLS > 0
//      block_new->SetGhostCells(SERVER_NUM_GHOST_CELLS);
//#endif

// Receive the block in 4 parts (member data plus 3 dynamic arrays). This is called even if SERVER_INTERP_ORDER is -1 to import the block dimensions
      MPI_Recv(block_new.get(), 1, MPIBlockType, 0, MPI::tag::sendblock, MPI::node_comm, MPI_STATUS_IGNORE);

#if SERVER_INTERP_ORDER > -1
      MPI_Recv(block_new->GetVariablesAddress(), block_new->GetVariableCount() * block_new->GetZoneCount(), MPI_DOUBLE, 0,
               MPI::tag::sendblock, MPI::node_comm, MPI_STATUS_IGNORE);
#endif
#if SERVER_INTERP_ORDER > 0 && SERVER_NUM_GHOST_CELLS == 0
      MPI_Recv(block_new->GetNeighborNodesAddress(), block_new->GetNeighborCount(), MPI_INT, 0,
               tag_sendblock, MPI_Config::node_comm, MPI_STATUS_IGNORE);
      MPI_Recv(block_new->GetNeighborLevelsAddress(), block_new->GetNeighborLevelCount(), MPI_INT, 0,
               tag_sendblock, MPI_Config::node_comm, MPI_STATUS_IGNORE);
#endif

// Insert the block into the cache
      block_new->ConfigureProperties();
      cache_line.AddBlock(block_new);
      bidx = block_new->GetNode();
   };

   return bidx;
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::PrintStencilOutcomes(void)
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
void BackgroundDataCartesian<HConfig>::PrintNumBlocksRequested(void)
{
   std::cerr << "Number of blocks requested: " << std::setw(10) << num_blocks_requested << std::endl;
};







//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesianFront templated methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/15/2025
\param[out] fields Fields
*/
template <typename HConfig>
template <typename Fields>
void BackgroundDataCartesian<HConfig>::GetVariablesFromReader(Fields& fields)
{
   int xyz;
   double rho, vars[Fields::size()] = {0.0};

// Get variables
   MPI_Send(&_inquiry, 1, MPIInquiryType, 0, MPI::tag::needvars, MPI::node_comm);
   MPI_Recv(vars, DataFields::size(), MPI_DOUBLE, 0, MPI::tag::sendvars, MPI::node_comm, MPI_STATUS_IGNORE);
   stencil_outcomes[2]++;

// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
   rho = vars[SERVER_VAR_INDEX_RHO];
#endif

// Number density, if provided
#ifdef SERVER_VAR_INDEX_DEN
   // Note: fields.Den(), formerly spdata.n_dens
   if constexpr (Fields::Den_found())
      fields.Den() = vars[SERVER_VAR_INDEX_DEN];
#endif

// Thermal pressure, if provided
// Note: fields.Prs(), formerly spdata.p_ther
#ifdef SERVER_VAR_INDEX_PTH
   if constexpr (Fields::Prs_found())
      fields.Prs() = vars[SERVER_VAR_INDEX_PTH];
#endif

   if constexpr (Fields::Vel_found()) {
// Bulk flow from mass density and momentum, if provided
// Note: fields.Vel(), formerly spdata.Uvec
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
      for (xyz = 0; xyz < 3; ++xyz)
         fields.Vel()[xyz] = vars[SERVER_VAR_INDEX_MOM + xyz] / rho;
#elif defined(SERVER_VAR_INDEX_FLO)
      for (xyz = 0; xyz < 3; ++xyz)
         fields.Vel()[xyz] = vars[SERVER_VAR_INDEX_FLO + xyz];
#else
      fields.Vel() = gv_zeros;
#endif
   }

// The magnetic field must be always provided
   for (xyz = 0; xyz < 3; ++xyz)
      fields.Mag()[xyz] = vars[SERVER_VAR_INDEX_MAG + xyz];

   if constexpr (Fields::Elc_found()) {
// Electric field, if provided
#if defined(SERVER_VAR_INDEX_ELE)
      for (xyz = 0; xyz < 3; ++xyz)
         fields.Elc()[xyz] = vars[SERVER_VAR_INDEX_ELE + xyz];
#elif defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      fields.Elc() = gv_zeros;
#endif
      // Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      fields.Elc() = -(fields.Vel() ^ fields.Mag()) / c_code;
#endif
#endif
   }

// TODO: after spdata-fields update, this section is awkward, and can be revised as needed. Old spdata section is commented out (not removed entirely).

// Region(s) indicator variable(s), if provided
//#ifdef SERVER_VAR_INDEX_REG
//   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++) spdata.region[xyz] = vars[SERVER_VAR_INDEX_REG + xyz];
//#else
//   spdata.region = gv_zeros;
//#endif

// Region(s) indicator variable(s), if provided
   if constexpr (Fields::Iv0_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv0() = vars[SERVER_VAR_INDEX_REG + 0];
#else
      fields.Iv0() = 0.0;
#endif
   }
   if constexpr (Fields::Iv1_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv1() = vars[SERVER_VAR_INDEX_REG + 1];
#else
      fields.Iv1() = 0.0;
#endif
   }
   if constexpr (Fields::Iv2_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv2() = vars[SERVER_VAR_INDEX_REG + 2];
#else
      fields.Iv2() = 0.0;
#endif
   }

};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 12/01/2023
\param[in]  pos    Position
\param[in]  block  Block containing pos
\param[out] fields Fields
*/
template <typename HConfig>
template <typename Fields>
void BackgroundDataCartesian<HConfig>::GetVariablesInterp0(const GeoVector& pos, Fields& fields)
{
   int xyz;
   double rho;

// Take the nearest cell value
   MultiIndex zone = block_pri->GetZone(pos);

// Spectrum-v1 pattern:
//   if constexpr (Fields::Den_found())
//      fields.Den() = block_pri->GetValue(zone, SERVER_VAR_INDEX_DEN);

   // todo this only needs to be rewritten in the format assuming a RequestedFields, like the ordinary Evaluate() methods

   // mass density
   if constexpr (DataFields::MassDen_found()) {
      fields.MassDen('w') = theblock[zone].Rho();
   }
   // number density
   if constexpr (DataFields::Den_found()) {
      fields.Den('w') = theblock[zone].Den();
   }
   // thermal pressure
   if constexpr (DataFields::Prs_found()) {
      fields.Prs('w') = theblock[zone].Prs();
   }
   // bulk flow velocity
   if constexpr (DataFields::Fluv_found()) {
      fields.Fluv('w') = theblock[zone].Fluv();
      // todo if Fluv is not explicitly stored/served, it can be calculated as Mom()/MassDen()
   }
   if constexpr (DataFields::Mag_found()) {
      fields.Mag() = theblock[zone].Mag();
   }


   if constexpr (Fields::Elc_found()) {
// Electric field, if provided
#if defined(SERVER_VAR_INDEX_ELE)
      for (xyz = 0; xyz < 3; xyz++)
         fields.Elc()[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_ELE + xyz);
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      fields.Elc() = gv_zeros;
#endif
// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      fields.Elc() = -(fields.Vel() ^ fields.Mag()) / c_code;
#endif
#endif
   }

// TODO: after spdata-fields update, this section is awkward, and can be revised as needed. Old spdata section is commented out (not removed entirely).

//// Region(s) indicator variable(s), if provided
//#ifdef SERVER_VAR_INDEX_REG
//   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++)
//       spdata.region[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_REG + xyz);
//#else
//   spdata.region = gv_zeros;
//#endif

// Region(s) indicator variable(s), if provided
   if constexpr (Fields::Iv0_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv0() = block_pri->GetValue(zone, SERVER_VAR_INDEX_REG + 0);
#else
      fields.Iv0() = 0.0;
#endif
   }
   if constexpr (Fields::Iv1_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv1() = block_pri->GetValue(zone, SERVER_VAR_INDEX_REG + 1);
#else
      fields.Iv1() = 0.0;
#endif
   }
   if constexpr (Fields::Iv2_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv2() = block_pri->GetValue(zone, SERVER_VAR_INDEX_REG + 2);
#else
      fields.Iv2() = 0.0;
#endif
   }


};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/25/2025
\param[in]  pos    Position
\param[out] spdata Fields, dmax, etc.
*/
template <typename HConfig>
template <typename Fields>
void BackgroundDataCartesian<HConfig>::GetVariablesInterp1(const GeoVector& pos, Fields& fields)
{
   int xyz, iz, vidx, pri_idx, sec_idx;
   // todo vars should be a Fields structured type
   double rho, var, vars[Fields::size()] = {0.0}, _Bmag2;

// Build the stencil. This is a time consuming operation.
   stencil_status = BuildInterpolationStencil(pos);
   fields.AbsMag() = 0.0;

// Internal interpolation
   if (stencil_status == 0) {
      stencil_outcomes[0]++;
      for (iz = 0; iz < stencil.n_elements; iz++) {
         for (vidx = 0; vidx < DataFields::size(); vidx++) {
            var = theblock[stencil.zones[iz]].Array()[vidx];
            // replaces::::::
//            var = block_pri->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         fields.AbsMag() += stencil.weights[iz] * sqrt(_Bmag2);
      };
   }

// Edge or corner interpolation
   else if (stencil_status == 4) {
      stencil_outcomes[2]++;
      pri_idx = block_pri->GetNode();
      sec_idx = block_sec->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else if (stencil.blocks[iz] == sec_idx) block_stn = block_sec;
         else block_stn = cache_line[stencil.blocks[iz]];
         for (vidx = 0; vidx < DataFields::size(); vidx++) {
            var = theblock[stencil.zones[iz]].Array()[vidx];
            // replaces::::
//            var = block_stn->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         fields.AbsMag() += stencil.weights[iz] * sqrt(_Bmag2);
      };
   }

// Plane interpolation
   else {
      stencil_outcomes[1]++;
      pri_idx = block_pri->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else block_stn = block_sec;
         for (vidx = 0; vidx < DataFields::size(); vidx++) {
            var = theblock[stencil.zones[iz]].Array()[vidx];
//            var = block_stn->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         fields.AbsMag() += stencil.weights[iz] * sqrt(_Bmag2);
      };
   };

// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
   rho = vars[SERVER_VAR_INDEX_RHO];
#endif

// Number density, if provided
// Note: fields.Den(), formerly spdata.n_dens
#ifdef SERVER_VAR_INDEX_DEN
   if constexpr (Fields::Den_found())
      fields.Den() = vars[SERVER_VAR_INDEX_DEN];
#endif

// Thermal pressure, if provided
// Note: fields.Prs(), formerly spdata.p_ther
#ifdef SERVER_VAR_INDEX_PTH
   if constexpr (Fields::Prs_found())
      fields.Prs() = vars[SERVER_VAR_INDEX_PTH];
#endif

// Convert the variables to SPECTRUM format

// Bulk flow from mass density and momentum, if provided
// Note: fields.Vel(), formerly spdata.Uvec
   if constexpr (Fields::Vel_found()) {
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
      for (xyz = 0; xyz < 3; xyz++)
         fields.Vel()[xyz] = vars[SERVER_VAR_INDEX_MOM + xyz] / rho;
#elif defined(SERVER_VAR_INDEX_FLO)
      for (xyz = 0; xyz < 3; xyz++)
         fields.Vel()[xyz] = vars[SERVER_VAR_INDEX_FLO + xyz];
#else
      fields.Vel() = gv_zeros;
#endif
   }

// The magnetic field must be always provided
   for (xyz = 0; xyz < 3; xyz++)
      fields.Mag()[xyz] = vars[SERVER_VAR_INDEX_MAG + xyz];

// Electric field, if provided
   if constexpr (Fields::Elc_found()) {
#if defined(SERVER_VAR_INDEX_ELE)
      for (xyz = 0; xyz < 3; xyz++)
         fields.Elc()[xyz] = vars[SERVER_VAR_INDEX_ELE + xyz];
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      fields.Elc()[xyz] = gv_zeros;
#endif
// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      fields.Elc() = -(fields.Vel() ^ fields.Mag()) / c_code;
#endif
#endif
   }

// TODO: after spdata-fields update, this section is awkward, and can be revised as needed. Old spdata section is commented out (not removed entirely).

// Region(s) indicator variable(s), if provided
//#ifdef SERVER_VAR_INDEX_REG
//   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++)
//      spdata.region[xyz] = vars[SERVER_VAR_INDEX_REG + xyz];
//#else
//   spdata.region = gv_zeros;
//#endif

// Region(s) indicator variable(s), if provided
   if constexpr (Fields::Iv0_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv0() = vars[SERVER_VAR_INDEX_REG + 0];
#else
      fields.Iv0() = 0.0;
#endif
   }
   if constexpr (Fields::Iv1_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv1() = vars[SERVER_VAR_INDEX_REG + 1];
#else
      fields.Iv1() = 0.0;
#endif
   }
   if constexpr (Fields::Iv2_found()) {
#ifdef SERVER_VAR_INDEX_REG
      fields.Iv2() = vars[SERVER_VAR_INDEX_REG + 2];
#else
      fields.Iv2() = 0.0;
#endif
   }

};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
\param[in]  t      Time
\param[in]  pos    Position
\param[out] fields Fields
\param[out] dmax dmax, field dependent
*/
template <typename HConfig>
template <typename Fields>
void BackgroundDataCartesian<HConfig>::GetVariables(double t, const GeoVector& pos, Fields& fields, double& dmax)
{
   int bidx;

// Request the block and the zone size
   _inquiry.type = 1;
   _inquiry.pos = pos;
   bidx = RequestBlock();

// If "block_pri" or "block_sec" is the position owner (based on the call to RequestBlock), we don't need to acccess the cache
   if (block_pri->GetNode() != bidx) {
      if (block_sec->GetNode() == bidx) block_pri = block_sec;
      else block_pri = cache_line[bidx];
   };
   dmax = fmin(dmax, block_pri->GetZoneLength().Smallest());

#if SERVER_INTERP_ORDER == -1
   // Get variables directly from reader program
   GetVariablesFromReader(fields);
#elif SERVER_INTERP_ORDER == 0
// Get variables using 0th order interpolation
   GetVariablesInterp0(pos, fields);
#elif SERVER_INTERP_ORDER == 1
   // Get variables using 1st order interpolation
   GetVariablesInterp1(pos, fields);
#else
#error Unsupported interpolation order!
#endif

// TODO: after spdata-fields update, this section is awkward, and can be revised as needed. Old spdata section is commented out (not removed entirely).

//// Perform unit conversion for fields and region
//#ifdef SERVER_VAR_INDEX_DEN
//   fields.region /= spdata.n_dens;
//#endif

#ifdef SERVER_VAR_INDEX_DEN
   if constexpr (Fields::Iv0_found())
      fields.Iv0() /= fields.Den();
   if constexpr (Fields::Iv1_found())
      fields.Iv1() /= fields.Den();
   if constexpr (Fields::Iv2_found())
      fields.Iv2() /= fields.Den();
#endif


   fields.Den() *= unit_number_density_server / unit_number_density_fluid;
   fields.Vel() *= unit_velocity_server / unit_velocity_fluid;
   fields.Mag() *= unit_magnetic_server / unit_magnetic_fluid;
   fields.Elc() *= unit_electric_server / unit_electric_fluid;
   fields.Prs() *= unit_pressure_server / unit_pressure_fluid;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 03/11/2024
\param[out] fields Fields
*/
template <typename HConfig>
template <typename Fields>
void BackgroundDataCartesian<HConfig>::GetGradientsInterp1(Fields& fields, DerivativeData& ddata)
{
   double var, _Bmag2, rho = 0.0;
   int vidx, xyz, uvw, iz, pri_idx, sec_idx;

   double grads[Fields::size()][3] = {0.0};
   fields.DelAbsMag() = gv_zeros;

   ddata.BACKGROUND_grad_FAIL = false;

// Internal interpolation
   if (stencil_status == 0) {
      for (iz = 0; iz < stencil.n_elements; iz++) {
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_pri->GetValue(stencil.zones[iz], vidx);
            for (uvw = 0; uvw < 3; uvw++) grads[vidx][uvw] += stencil.derivatives[3 * iz + uvw] * var;
         };
// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
         rho += stencil.weights[iz] * block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_RHO);
#endif
// gradient of B magnitude
         _Bmag2 = Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) fields.DelMag()[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
      };
   }

// Edge or corner interpolation
   else if (stencil_status == 4) {
      pri_idx = block_pri->GetNode();
      sec_idx = block_sec->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else if (stencil.blocks[iz] == sec_idx) block_stn = block_sec;
         else block_stn = cache_line[stencil.blocks[iz]];
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_stn->GetValue(stencil.zones[iz], vidx);
            for (uvw = 0; uvw < 3; uvw++) grads[vidx][uvw] += stencil.derivatives[3 * iz + uvw] * var;
         };
// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
         rho += stencil.weights[iz] * block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_RHO);
#endif
// gradient of B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) fields.DelAbsMag()[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
      };
   }

// Plane interpolation
   else {
      pri_idx = block_pri->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else block_stn = block_sec;
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_stn->GetValue(stencil.zones[iz], vidx);
            for (uvw = 0; uvw < 3; uvw++) grads[vidx][uvw] += stencil.derivatives[3 * iz + uvw] * var;
         };
// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
         rho += stencil.weights[iz] * block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_RHO);
#endif
// gradient of B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) fields.DelAbsMag()[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
      };
   };

// Convert the gradients to SPECTRUM format
   for (uvw = 0; uvw < 3; uvw++) {
      for (xyz = 0; xyz < 3; xyz++) {

         if constexpr (Fields::DelVel_found()) {
// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
            fields.DelVel()[uvw][xyz] =
                  (grads[SERVER_VAR_INDEX_MOM + xyz][uvw] - fields.Vel()[xyz] * grads[SERVER_VAR_INDEX_RHO][uvw]) / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
            fields.DelVel()[uvw][xyz] = grads[SERVER_VAR_INDEX_FLO + xyz][uvw];
#else
            fields.DelVel()[uvw][xyz] = 0.0;
#endif
         }

// The magnetic field must be always provided
         if constexpr (Fields::DelMag_found())
            fields.DelMag()[uvw][xyz] = grads[SERVER_VAR_INDEX_MAG + xyz][uvw];

// Electric field, if provided
         if constexpr (Fields::DelElc_found()) {
#ifdef SERVER_VAR_INDEX_ELE
            fields.DelElc()[uvw][xyz] = grads[SERVER_VAR_INDEX_ELE + xyz][uvw];
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
            fields.DelElc()[uvw][xyz] = 0.0;
#endif
         }

      };
   };

   if constexpr (Fields::DelElc_found()) {
// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      fields.DelElc() = -((fields.DelVel() ^ fields.Mag()) + (fields.Vel() ^ fields.DelMag())) / c_code;
#endif
#endif
   }

};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param[out] spdata Field gradients
*/
template <typename HConfig>
template <typename Fields>
void BackgroundDataCartesian<HConfig>::GetGradients(Fields& fields, DerivativeData& ddata)
{
#if SERVER_INTERP_ORDER == -1
   // Gradients must be computed numerically
   ddata.BACKGROUND_grad_FAIL = true;
   return;
#elif SERVER_INTERP_ORDER == 0
// All gradients are explicitly set to zero, and the background must not attempt to compute them using "NumericalDerivatives()"
   fields.DelVel() = gm_zeros;
   fields.DelMag() = gm_zeros;
   fields.DelAbsMag() = gv_zeros;
   fields.DelElc() = gm_zeros;
#elif SERVER_INTERP_ORDER == 1
   // Gradients can be obtained from the stencil construction
   GetGradientsInterp1(fields, ddata);
#else
#error Unsupported interpolation order!
#endif

// Perform unit conversion for gradients
   fields.DelVel() *= unit_velocity_server / unit_velocity_fluid;
   fields.DelMag() *= unit_magnetic_server / unit_magnetic_fluid;
   fields.DelElc() *= unit_electric_server / unit_electric_fluid;
};




};
