/*!
\file server_server.cc
\brief Implements a class of a data server for a uniform Cartesian grid
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_cartesian.hh"
#include "block_cartesian.hh"
#include "reader_cartesian.hh"
#include "common/print_warn.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void StencilCartesian::Print(void)
{
   for (auto iz = 0; iz < n_elements; iz++) {
      std::cerr << std::setw(10) << blocks[iz] << "  " << zones[iz] << "  " << weights[iz] << std::endl;
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesian methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
\param[out] block_ptr pointer to block type
*/
// TODO: look up all calls to this function for possible replacement with a simple "new"
void ServerCartesian::InitializeBlockPtr(BlockBase* &block_ptr)
{
   block_ptr = new BlockCartesian();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesian::CreateBlockDatatype(void)
{
   int i;

// Temporary object for displacement calculations
   BlockBase* block_tmp;
   InitializeBlockPtr(block_tmp);
   n_variables = block_tmp->GetVariableCount();

// Set up MPI data type for "BlockXXXX" (without the dynamic storage)
   MPI_Datatype block_types[] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
   int block_lengths[] = {1, 3, 3, 3, 3};
   MPI_Aint block_start, block_displ[5];

// Figure out field displacements using "block_tmp" as template
   MPI_Get_address(block_tmp, &block_start);
   MPI_Get_address(block_tmp->GetNodeAddress(), &block_displ[0]);
   MPI_Get_address(block_tmp->GetFaceMinAddress(), &block_displ[1]);
   MPI_Get_address(block_tmp->GetFaceMaxAddress(), &block_displ[2]);
   MPI_Get_address(block_tmp->GetFaceMinPhysAddress(), &block_displ[3]);
   MPI_Get_address(block_tmp->GetFaceMaxPhysAddress(), &block_displ[4]);
   for (i = 4; i >= 0; i--) block_displ[i] -= block_start;

// Commit the type
   MPI_Type_create_struct(5, block_lengths, block_displ, block_types, &MPIBlockType);
   MPI_Type_commit(&MPIBlockType);
   delete block_tmp;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
*/
void ServerCartesian::CreateStencilDatatype(void)
{
   int i;

// Set up MPI data type for "InterpolationStencil"
   MPI_Datatype stencil_types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
   int stencil_lengths[] = {1, 8, 24, 8, 24};
   MPI_Aint stencil_start, stencil_displ[5];

// Figure out field displacements using "stencil" for template
   MPI_Get_address(&stencil, &stencil_start);
   MPI_Get_address(&stencil.n_elements , &stencil_displ[0]);
   MPI_Get_address(&stencil.blocks     , &stencil_displ[1]);
   MPI_Get_address(&stencil.zones      , &stencil_displ[2]);
   MPI_Get_address(&stencil.weights    , &stencil_displ[3]);
   MPI_Get_address(&stencil.derivatives, &stencil_displ[4]);
   for (i = 4; i >= 0; i--) stencil_displ[i] -= stencil_start;

// Commit the type
   MPI_Type_create_struct(5, stencil_lengths, stencil_displ, stencil_types, &MPIStencilType);
   MPI_Type_commit(&MPIStencilType);
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesian::CreateMPIDatatypes(void)
{
// Create block datatype
   CreateBlockDatatype();
// Create stencil datatype
   CreateStencilDatatype();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesian::ServerStart(void)
{
// Call the parent version
   ServerBase::ServerStart();
   CreateMPIDatatypes();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesian::ServerFinish(void)
{
   MPI_Type_free(&MPIBlockType);
   MPI_Type_free(&MPIStencilType);

// Call the parent version
   ServerBase::ServerFinish();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesianFront methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
void ServerCartesianFront::ServerStart(void) {

// No need to call the ServerBaseFront version because it merely calls the ServerBase version
   ServerCartesian::ServerStart();

   cache_line.Empty();
   stencil_outcomes[0] = stencil_outcomes[1] = stencil_outcomes[2] = 0;
   num_blocks_requested = 0;

   MPI_Bcast(domain_min.Data(), 3, MPI_DOUBLE, 0, MPI_Config::node_comm);
   MPI_Bcast(domain_max.Data(), 3, MPI_DOUBLE, 0, MPI_Config::node_comm);

// Prime "block_pri" and "block_sec" with stub blocks that always fail tests. These and "block_stn" must be smart pointers to avoid double free or corruption errors.
   MakeSharedBlock(block_pri);
   block_pri->SetDimensions(domain_max, domain_min);
   block_pri->BlockBase::LoadDimensions(1.0);
   MakeSharedBlock(block_sec);
   block_sec->SetDimensions(domain_max, domain_min);
   block_sec->BlockBase::LoadDimensions(1.0);
   MakeSharedBlock(block_stn);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesianFront::ServerFinish(void)
{
   MPI_Send(nullptr, 0, MPI_BYTE, 0, tag_stopserve, MPI_Config::node_comm);
   ServerCartesian::ServerFinish();
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\param[out] block_new pointer to block type
*/
void ServerCartesianFront::MakeSharedBlock(BlockPtrType &block_new)
{
   block_new = std::make_shared<BlockCartesian>();
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\param[in] zone_lo   lower interpolation zone
\param[in] zone_hi   higher interpolation zone
\param[in] offset_lo lower interpolation offset
\param[in] offset_hi higher interpolation offset
*/
void ServerCartesianFront::InteriorInterpolationStencil(const MultiIndex zone_lo, const MultiIndex zone_hi,
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
int ServerCartesianFront::BuildInterpolationStencil(const GeoVector& pos)
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

#if SERVER_NUM_GHOST_CELLS == 0
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
#endif

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
int ServerCartesianFront::RequestBlock(void)
{
   int bidx;
   BlockPtrType block_new;

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
      MPI_Send(&_inquiry, 1, MPIInquiryType, 0, tag_needblock, MPI_Config::node_comm);
      num_blocks_requested++;

// Allocate memory for block.
      MakeSharedBlock(block_new);

// Adjust for ghost cells if necessary
#if SERVER_NUM_GHOST_CELLS > 0
      block_new->SetGhostCells(SERVER_NUM_GHOST_CELLS);
#endif

// Receive the block in 4 parts (member data plus 3 dynamic arrays). This is called even if SERVER_INTERP_ORDER is -1 to import the block dimensions
      MPI_Recv(block_new.get(), 1, MPIBlockType, 0, tag_sendblock, MPI_Config::node_comm, MPI_STATUS_IGNORE);

#if SERVER_INTERP_ORDER > -1
      MPI_Recv(block_new->GetVariablesAddress(), block_new->GetVariableCount() * block_new->GetZoneCount(), MPI_DOUBLE, 0,
               tag_sendblock, MPI_Config::node_comm, MPI_STATUS_IGNORE);
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
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 12/01/2023
\param[out] spdata Fields, dmax, etc.
*/
void ServerCartesianFront::GetVariablesFromReader(SpatialData& spdata)
{
   int xyz;
   double rho, vars[n_variables] = {0.0};

// Get variables
   MPI_Send(&_inquiry, 1, MPIInquiryType, 0, tag_needvars, MPI_Config::node_comm);
   MPI_Recv(vars, n_variables, MPI_DOUBLE, 0, tag_sendvars, MPI_Config::node_comm, MPI_STATUS_IGNORE);
   stencil_outcomes[2]++;

// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
   rho = vars[SERVER_VAR_INDEX_RHO];
#endif

// Number density, if provided
#ifdef SERVER_VAR_INDEX_DEN
   spdata.n_dens = vars[SERVER_VAR_INDEX_DEN];
#endif

// Thermal pressure, if provided
#ifdef SERVER_VAR_INDEX_PTH
   spdata.p_ther = vars[SERVER_VAR_INDEX_PTH];
#endif

   for (xyz = 0; xyz < 3; xyz++) {

// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
      spdata.Uvec[xyz] = vars[SERVER_VAR_INDEX_MOM + xyz] / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
      spdata.Uvec[xyz] = vars[SERVER_VAR_INDEX_FLO + xyz];
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = vars[SERVER_VAR_INDEX_MAG + xyz];

// Electric field, if provided
#if defined(SERVER_VAR_INDEX_ELE)
      spdata.Evec[xyz] = vars[SERVER_VAR_INDEX_ELE + xyz];
#elif defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)) 
      spdata.Evec[xyz] = 0.0;
#endif

   };

// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))   
   spdata.Evec = -(spdata.Uvec ^ spdata.Bvec) / c_code;
#endif
#endif

// Region(s) indicator variable(s), if provided
#ifdef SERVER_VAR_INDEX_REG
   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++) spdata.region[xyz] = vars[SERVER_VAR_INDEX_REG + xyz];
#else
   spdata.region = gv_zeros;
#endif
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 12/01/2023
\param[in]  pos    Position
\param[in]  block  Block containing pos
\param[out] spdata Fields, dmax, etc.
*/
void ServerCartesianFront::GetVariablesInterp0(const GeoVector& pos, SpatialData& spdata)
{
   int xyz;
   double rho;

// Take the nearest cell value
   MultiIndex zone = block_pri->GetZone(pos);

// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
   rho = block_pri->GetValue(zone, SERVER_VAR_INDEX_RHO);
#endif

// Number density, if provided
#ifdef SERVER_VAR_INDEX_DEN
   spdata.n_dens = block_pri->GetValue(zone, SERVER_VAR_INDEX_DEN);
#endif

// Thermal pressure, if provided
#ifdef SERVER_VAR_INDEX_PTH
   spdata.p_ther = block_pri->GetValue(zone, SERVER_VAR_INDEX_PTH);
#endif

// Convert the variables to SPECTRUM format
   for (xyz = 0; xyz < 3; xyz++) {

// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
      spdata.Uvec[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_MOM + xyz) / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
      spdata.Uvec[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_FLO + xyz);
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_MAG + xyz);

// Electric field, if provided
#if defined(SERVER_VAR_INDEX_ELE)
      spdata.Evec[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_ELE + xyz);
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      spdata.Evec[xyz] = 0.0;
#endif

   };

// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))   
   spdata.Evec = -(spdata.Uvec ^ spdata.Bvec) / c_code;
#endif
#endif

// Region(s) indicator variable(s), if provided
#ifdef SERVER_VAR_INDEX_REG
   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++) spdata.region[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_REG + xyz);
#else
   spdata.region = gv_zeros;
#endif
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 12/01/2023
\param[in]  pos    Position
\param[out] spdata Fields, dmax, etc.
*/
void ServerCartesianFront::GetVariablesInterp1(const GeoVector& pos, SpatialData& spdata)
{
   int xyz, iz, vidx, pri_idx, sec_idx;
   double rho, var, vars[n_variables] = {0.0}, _Bmag2;

// Build the stencil. This is a time consuming operation.
   stencil_status = BuildInterpolationStencil(pos);
   spdata.Bmag = 0.0;

// Internal interpolation
   if (stencil_status == 0) {
      stencil_outcomes[0]++;
      for (iz = 0; iz < stencil.n_elements; iz++) {
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_pri->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         spdata.Bmag += stencil.weights[iz] * sqrt(_Bmag2);
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
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_stn->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         spdata.Bmag += stencil.weights[iz] * sqrt(_Bmag2);
      };
   }

// Plane interpolation
   else {
      stencil_outcomes[1]++;
      pri_idx = block_pri->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else block_stn = block_sec;
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_stn->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         spdata.Bmag += stencil.weights[iz] * sqrt(_Bmag2);
      };
   };

// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
   rho = vars[SERVER_VAR_INDEX_RHO];
#endif

// Number density, if provided
#ifdef SERVER_VAR_INDEX_DEN
   spdata.n_dens = vars[SERVER_VAR_INDEX_DEN];
#endif

// Thermal pressure, if provided
#ifdef SERVER_VAR_INDEX_PTH
   spdata.p_ther = vars[SERVER_VAR_INDEX_PTH];
#endif

// Convert the variables to SPECTRUM format
   for (xyz = 0; xyz < 3; xyz++) {

// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
      spdata.Uvec[xyz] = vars[SERVER_VAR_INDEX_MOM + xyz] / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
      spdata.Uvec[xyz] = vars[SERVER_VAR_INDEX_FLO + xyz];
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = vars[SERVER_VAR_INDEX_MAG + xyz];

// Electric field, if provided
#if defined(SERVER_VAR_INDEX_ELE)
      spdata.Evec[xyz] = vars[SERVER_VAR_INDEX_ELE + xyz];
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      spdata.Evec[xyz] = 0.0;
#endif

   };

// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))   
   spdata.Evec = -(spdata.Uvec ^ spdata.Bvec) / c_code;
#endif
#endif

// Region(s) indicator variable(s), if provided
#ifdef SERVER_VAR_INDEX_REG
   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++) spdata.region[xyz] = vars[SERVER_VAR_INDEX_REG + xyz];
#else
   spdata.region = gv_zeros;
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
\param[in]  t      Time
\param[in]  pos    Position
\param[out] spdata Fields, dmax, etc.
*/
void ServerCartesianFront::GetVariables(double t, const GeoVector& pos, SpatialData& spdata)
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
   spdata.dmax = fmin(spdata.dmax, block_pri->GetZoneLength().Smallest());

#if SERVER_INTERP_ORDER == -1
// Get variables directly from reader program
   GetVariablesFromReader(spdata);
#elif SERVER_INTERP_ORDER == 0
// Get variables using 0th order interpolation
   GetVariablesInterp0(pos, spdata);
#elif SERVER_INTERP_ORDER == 1
// Get variables using 1st order interpolation
   GetVariablesInterp1(pos, spdata);
#else
#error Unsupported interpolation order!
#endif

// Perform unit conversion for fields and region
#ifdef SERVER_VAR_INDEX_DEN
   spdata.region /= spdata.n_dens;
#endif
   spdata.n_dens *= unit_number_density_server / unit_number_density_fluid;
   spdata.Uvec *= unit_velocity_server / unit_velocity_fluid;
   spdata.Bvec *= unit_magnetic_server / unit_magnetic_fluid;
   spdata.Evec *= unit_electric_server / unit_electric_fluid;
   spdata.p_ther *= unit_pressure_server / unit_pressure_fluid;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 03/11/2024
\param[out] spdata Fields, dmax, etc.
*/
void ServerCartesianFront::GetGradientsInterp1(SpatialData& spdata)
{
   double var, _Bmag2, rho = 0.0;
   int vidx, xyz, uvw, iz, pri_idx, sec_idx;

   double grads[n_variables][3] = {0.0};
   spdata.gradBmag = gv_zeros;

   LOWER_BITS(spdata._mask, BACKGROUND_grad_FAIL);

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
         for (uvw = 0; uvw < 3; uvw++) spdata.gradBmag[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
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
         for (uvw = 0; uvw < 3; uvw++) spdata.gradBmag[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
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
         for (uvw = 0; uvw < 3; uvw++) spdata.gradBmag[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
      };
   };

// Convert the gradients to SPECTRUM format
   for (uvw = 0; uvw < 3; uvw++) {
      for (xyz = 0; xyz < 3; xyz++) {

// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
         spdata.gradUvec[uvw][xyz] = (grads[SERVER_VAR_INDEX_MOM + xyz][uvw] - spdata.Uvec[xyz] * grads[SERVER_VAR_INDEX_RHO][uvw]) / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
         spdata.gradUvec[uvw][xyz] = grads[SERVER_VAR_INDEX_FLO + xyz][uvw];
#else
         spdata.gradUvec[uvw][xyz] = 0.0;
#endif

// The magnetic field must be always provided
         spdata.gradBvec[uvw][xyz] = grads[SERVER_VAR_INDEX_MAG + xyz][uvw];

// Electric field, if provided
#ifdef SERVER_VAR_INDEX_ELE
         spdata.gradEvec[uvw][xyz] = grads[SERVER_VAR_INDEX_ELE + xyz][uvw];
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
         spdata.gradEvec[uvw][xyz] = 0.0;
#endif

      };
   };

// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
   spdata.gradEvec = -((spdata.gradUvec ^ spdata.Bvec) + (spdata.Uvec ^ spdata.gradBvec)) / c_code;
#endif
#endif

};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param[out] spdata Field gradients
*/
void ServerCartesianFront::GetGradients(SpatialData& spdata)
{
#if SERVER_INTERP_ORDER == -1
// Gradients must be computed numerically
   RAISE_BITS(spdata._mask, BACKGROUND_grad_FAIL);
   return;
#elif SERVER_INTERP_ORDER == 0
// All gradients are explicitly set to zero, and the background must not attempt to compute them using "NumericalDerivatives()"
   spdata.gradUvec = gm_zeros;
   spdata.gradBvec = gm_zeros;
   spdata.gradBmag = gv_zeros;
   spdata.gradEvec = gm_zeros;
#elif SERVER_INTERP_ORDER == 1
// Gradients can be obtained from the stencil construction
   GetGradientsInterp1(spdata);
#else
#error Unsupported interpolation order!
#endif

// Perform unit conversion for gradients
   spdata.gradUvec *= unit_velocity_server / unit_velocity_fluid;
   spdata.gradBvec *= unit_magnetic_server / unit_magnetic_fluid;
   spdata.gradEvec *= unit_electric_server / unit_electric_fluid;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesianFront::PrintStencilOutcomes(void)
{
   std::cerr << "Stencil outcomes: " << std::setw(10) << stencil_outcomes[0]
                                     << std::setw(10) << stencil_outcomes[1]
                                     << std::setw(10) << stencil_outcomes[2] << std::endl;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
*/
void ServerCartesianFront::PrintNumBlocksRequested(void)
{
   std::cerr << "Number of blocks requested: " << std::setw(10) << num_blocks_requested << std::endl;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesianBack public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] file_name_pattern_in A string describing the file naming pattern
*/
ServerCartesianBack::ServerCartesianBack(const std::string& file_name_pattern_in)
                   : ServerBaseBack(file_name_pattern_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/01/2023
*/
void ServerCartesianBack::ReadData(const std::string data_file)
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
void ServerCartesianBack::ServerStart(void)
{
   ServerBaseBack::ServerStart();
   CreateMPIDatatypes();
   InitializeBlockPtr(block_served);

// Adjust for ghost cells if necessary
#if SERVER_NUM_GHOST_CELLS > 0
   block_served->SetGhostCells(SERVER_NUM_GHOST_CELLS);
#endif

// Initialize the Cartesian library and read the data into memory
   std::string data_file = file_name_pattern + ".out";

   ReadData(data_file);
   domain_min *= unit_length_server / unit_length_fluid;
   domain_max *= unit_length_server / unit_length_fluid;

   MPI_Bcast(domain_min.Data(), 3, MPI_DOUBLE, 0, MPI_Config::node_comm);
   MPI_Bcast(domain_max.Data(), 3, MPI_DOUBLE, 0, MPI_Config::node_comm);
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
*/
void ServerCartesianBack::CleanReader(void)
{
   ReadCartesianClean();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesianBack::ServerFinish(void)
{
   CleanReader();
   delete block_served;

// This is a repeat of the function from ServerCartesian
   MPI_Type_free(&MPIBlockType);
   MPI_Type_free(&MPIStencilType);

   ServerBaseBack::ServerFinish();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\return Number of clients that completed their tasks during this cycle

/Note "needvars" requests are only handled when SERVER_INTERP_ORDER = -1
*/
int ServerCartesianBack::ServerFunctions(void)
{
#if SERVER_INTERP_ORDER == -1
// Handle "needvars" requests
   HandleNeedVarsRequests();
#endif
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
void ServerCartesianBack::GetBlockData(const double* pos, double* vars, int* found)
{
   ReadCartesianGetBlockData(pos, vars, found);
};

/*!
\author Juan G Alonso Guzman
\date 08/03/2023
*/
void ServerCartesianBack::HandleNeedVarsRequests(void)
{
   GeoVector pos_cart;
   double vars[n_variables];
   int found, cpu, cpu_idx, count_needvars = 0;

// Service the "needvars" requests
   MPI_Testsome(MPI_Config::node_comm_size, req_needvars, &count_needvars, index_needvars, MPI_STATUSES_IGNORE);

   for (cpu_idx = 0; cpu_idx < count_needvars; cpu_idx++) {
      cpu = index_needvars[cpu_idx];

// Obtain the variables requested
      pos_cart = buf_needvars[cpu].pos / unit_length_server * unit_length_fluid;
      GetBlockData(pos_cart.Data(), vars, &found);
      if (!found) std::cerr << "Position not found\n";

// Send the variables to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(vars, n_variables, MPI_DOUBLE, cpu, tag_sendvars, MPI_Config::node_comm);

// Post the receive for the next variables request from this worker.
      MPI_Irecv(&buf_needvars[cpu], 1, MPIInquiryType, cpu, tag_needvars, MPI_Config::node_comm, &req_needvars[cpu]);
   };
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
\param[in] pos    position array in reader coordinates
\param[out] node  ID tag of block containing pos within its physical boundaries
*/
void ServerCartesianBack::GetBlock(const double* pos, int* node)
{
   ReadCartesianGetNode(pos, node);
};

/*!
\author Juan G Alonso Guzman
\date 12/01/2023
*/
void ServerCartesianBack::HandleNeedBlockRequests(void)
{
   GeoVector pos_cart;
   int cpu, cpu_idx, count_needblock = 0;

   // Service the "needblock" requests
   MPI_Testsome(MPI_Config::node_comm_size, req_needblock, &count_needblock, index_needblock, MPI_STATUSES_IGNORE);

// Load the block requested. If the requestor does not know the node, figure it out.
   for (cpu_idx = 0; cpu_idx < count_needblock; cpu_idx++) {
      cpu = index_needblock[cpu_idx];

      if (buf_needblock[cpu].type) {
         pos_cart = buf_needblock[cpu].pos / unit_length_server * unit_length_fluid;
         GetBlock(pos_cart.Data(), &buf_needblock[cpu].node);
         if (buf_needblock[cpu].node == -1) throw ExServerError();
      };

      block_served->SetNode(buf_needblock[cpu].node);
      block_served->LoadDimensions(unit_length_server);

// Send the block to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(block_served, 1, MPIBlockType, cpu, tag_sendblock, MPI_Config::node_comm);
#if SERVER_INTERP_ORDER > -1
      block_served->LoadVariables();
      MPI_Send(block_served->GetVariablesAddress(), block_served->GetVariableCount() * block_served->GetZoneCount(),
               MPI_DOUBLE, cpu, tag_sendblock, MPI_Config::node_comm);
#endif
#if SERVER_INTERP_ORDER > 0 && SERVER_NUM_GHOST_CELLS == 0
      block_served->LoadNeighbors();
      MPI_Send(block_served->GetNeighborNodesAddress(), block_served->GetNeighborCount(),
               MPI_INT, cpu, tag_sendblock, MPI_Config::node_comm);
      MPI_Send(block_served->GetNeighborLevelsAddress(), block_served->GetNeighborLevelCount(),
               MPI_INT, cpu, tag_sendblock, MPI_Config::node_comm);
#endif

// Post the receive for the next block request from this worker
      MPI_Irecv(&buf_needblock[cpu], 1, MPIInquiryType, cpu, tag_needblock, MPI_Config::node_comm, &req_needblock[cpu]);
   };
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\return Number of clients that completed their tasks during this cycle
*/
int ServerCartesianBack::HandleStopServeRequests(void)
{
   int cpu, cpu_idx, count_stopserve = 0;

   // Service the "stopserve" requests. We assume that each worker sends a single request at the end of the simulation. 
   MPI_Testsome(MPI_Config::node_comm_size, req_stopserve, &count_stopserve, index_stopserve, MPI_STATUSES_IGNORE);
   
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
