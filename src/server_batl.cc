/*!
\file server_batl.cc
\brief Implements a class of a data server from BATL
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_batl.hh"
#include "block_batl.hh"
#include "common/print_warn.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATL methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
\param[out] block_ptr pointer to block type
*/
void ServerBATL::InitializeBlockPtr(BlockBase* &block_ptr)
{
   block_ptr = new BlockBATL();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATLFront methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\param[out] block_new pointer to block type
*/
void ServerBATLFront::MakeSharedBlock(BlockPtrType &block_new)
{
   block_new = std::make_shared<BlockBATL>();
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 09/18/2020
*/
int ServerBATLFront::RequestStencil(const GeoVector& pos)
{
   int iz;

   _inquiry.type = 1;
   _inquiry.pos = pos;
   MPI_Send(&_inquiry, 1, MPIInquiryType, 0, tag_needstencil, MPI_Config::node_comm);
   MPI_Recv(&stencil, 1, MPIStencilType, 0, tag_sendstencil, MPI_Config::node_comm, MPI_STATUS_IGNORE);

// The "blocks" field of the returned stencil currently contains nodes. We need to convert them to blocks, requesting the blocks from the server, if necessary.
   _inquiry.type = 0;
   for (iz = 0; iz < stencil.n_elements; iz++) {
      _inquiry.node = stencil.blocks[iz];
      stencil.blocks[iz] = RequestBlock();
   };

   return 4;
};

/*!
\author Vladimir Florinski
\date 07/24/2020
\param[in] pos   Interpolation point position
\param[in] plane Direction of projection (0-2)
\param[in] half  Which half of the 3D stencil to work on (0 or 1)
\return Status: 0 if interpolation succeeded, -1 if not (neighbour blocks are at different levels)
*/
int ServerBATLFront::BuildInterpolationPlane(const GeoVector& pos, int plane, int half)
{
   int idx0, idx1, idx2, bidx_self, bidx_hori, bidx_vert, bidx_diag, offset;
   MultiIndex block_size, zones[4], node_idx, level_idx;
   GeoVector offset_lo, offset_hi, delta;

// Offset based on whether this is the first or the second half of the 3D stencil.
   offset = (half == 0 ? 0 : 4);

// Either 0 (primary) or 4 (secondary)
   bidx_self = stencil.blocks[offset];
   BlockPtrType block_self = cache_line[bidx_self];

   block_self->GetZoneOffset(pos, zones[0], offset_lo);
   offset_hi = 1.0 - offset_lo;
   delta = block_self->GetZoneLength();
   block_size = block_self->GetBlockSize();

// Index order. Center arrangement in a plane stencil (idx0 is out of the plane of the screen) is shown below:
//
//        idx2
//          ^
//          |
//
//          2---------3
//          |         |
//          |         |
//          |         |
//          0---------1   --> idx1

   idx0 = plane;
   idx1 = (plane + 1) % 3;
   idx2 = (plane + 2) % 3;

// Figure out the improper zones (no account for boundaries)
   zones[1] = zones[0];
   zones[1][idx1] += 1;

   zones[2] = zones[0];
   zones[2][idx2] += 1;

   zones[3] = zones[1];
   zones[3][idx2] += 1;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Compile a list of up to four neighbor blocks
//----------------------------------------------------------------------------------------------------------------------------------------------------

// The "idx0" component is always 1 because we work in a single plane.
   level_idx[idx0] = 1;
   _inquiry.type = 0;

// Check the left neighbor level. Return if it is different from the starting block's level. Otherwise, load the block if necessary.
   if (zones[0][idx1] < 0) {
      level_idx[idx1] = 0;
      level_idx[idx2] = 1;
      if (block_self->GetNeighborLevel(level_idx) != 0) return -1;
      node_idx = LevelToNode(level_idx);
      _inquiry.node = block_self->GetNeighborNode(node_idx);
      bidx_hori = RequestBlock();
   }

// Check the right neighbor level. Return if it is different from the starting block's level. Otherwise, load the block if necessary.
   else if (zones[1][idx1] >= block_size[idx1]) {
      level_idx[idx1] = 2;
      level_idx[idx2] = 1;
      if (block_self->GetNeighborLevel(level_idx) != 0) return -1;
      node_idx = LevelToNode(level_idx);
      _inquiry.node = block_self->GetNeighborNode(node_idx);
      bidx_hori = RequestBlock();
   };

// Check the bottom neighbor level. Return if it is different from the starting block's level. Otherwise, load the block if necessary.
   if (zones[0][idx2] < 0) {
      level_idx[idx1] = 1;
      level_idx[idx2] = 0;
      if (block_self->GetNeighborLevel(level_idx) != 0) return -1;
      node_idx = LevelToNode(level_idx);
      _inquiry.node = block_self->GetNeighborNode(node_idx);
      bidx_vert = RequestBlock();
   }

// Check the top neighbor level. Return if it is different from the starting block's level. Otherwise, load the block if necessary.
   else if (zones[2][idx2] >= block_size[2]) {
      level_idx[idx1] = 1;
      level_idx[idx2] = 2;
      if (block_self->GetNeighborLevel(level_idx) != 0) return -1;
      node_idx = LevelToNode(level_idx);
      _inquiry.node = block_self->GetNeighborNode(node_idx);
      bidx_vert = RequestBlock();
   };

// Check the lower left neighbor level. Return if it is different from the starting block's level. Otherwise, load the block if necessary.
   if ((zones[0][idx1] < 0) && (zones[0][idx2] < 0)) {
      level_idx[idx1] = 0;
      level_idx[idx2] = 0;
      if (block_self->GetNeighborLevel(level_idx) != 0) return -1;
      node_idx = LevelToNode(level_idx);
      _inquiry.node = block_self->GetNeighborNode(node_idx);
      bidx_diag = RequestBlock();
   }

// Check the lower right neighbor level. Return if it is different from the starting block's level. Otherwise, load the block if necessary.
   else if ((zones[1][idx1] >= block_size[1]) && (zones[1][idx2] < 0)) {
      level_idx[idx1] = 2;
      level_idx[idx2] = 0;
      if (block_self->GetNeighborLevel(level_idx) != 0) return -1;
      node_idx = LevelToNode(level_idx);
      _inquiry.node = block_self->GetNeighborNode(node_idx);
      bidx_diag = RequestBlock();
   }

// Check the upper left neighbor level. Return if it is different from the starting block's level. Otherwise, load the block if necessary.
   else if ((zones[2][idx1] < 0) && (zones[2][idx2] >= block_size[2])) {
      level_idx[idx1] = 0;
      level_idx[idx2] = 2;
      if (block_self->GetNeighborLevel(level_idx) != 0) return -1;
      node_idx = LevelToNode(level_idx);
      _inquiry.node = block_self->GetNeighborNode(node_idx);
      bidx_diag = RequestBlock();
   }

// Check the upper right neighbor level. Return if it is different from the starting block's level. Otherwise, load the block if necessary.
   else if ((zones[3][idx1] >= block_size[1]) && (zones[3][idx2] >= block_size[2])) {
      level_idx[idx1] = 2;
      level_idx[idx2] = 2;
      if (block_self->GetNeighborLevel(level_idx) != 0) return -1;
      node_idx = LevelToNode(level_idx);
      _inquiry.node = block_self->GetNeighborNode(node_idx);
      bidx_diag = RequestBlock();
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Sort out the order of the block in the plane stencil
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Left and right margins, including corners.
   if (zones[0][idx1] < 0) {
      if (zones[0][idx2] < 0) {
         stencil.blocks[0 + offset] = bidx_diag;
         stencil.blocks[1 + offset] = bidx_vert;
         stencil.blocks[2 + offset] = bidx_hori;
         stencil.blocks[3 + offset] = bidx_self;
      }
      else if (zones[2][idx2] >= block_size[idx2]) {
         stencil.blocks[0 + offset] = bidx_hori;
         stencil.blocks[1 + offset] = bidx_self;
         stencil.blocks[2 + offset] = bidx_diag;
         stencil.blocks[3 + offset] = bidx_vert;
      }
      else {
         stencil.blocks[0 + offset] = bidx_hori;
         stencil.blocks[1 + offset] = bidx_self;
         stencil.blocks[2 + offset] = bidx_hori;
         stencil.blocks[3 + offset] = bidx_self;
      };
   }

   else if (zones[1][idx1] >= block_size[idx1]) {
      if (zones[0][idx2] < 0) {
         stencil.blocks[0 + offset] = bidx_vert;
         stencil.blocks[1 + offset] = bidx_diag;
         stencil.blocks[2 + offset] = bidx_self;
         stencil.blocks[3 + offset] = bidx_hori;
      }
      else if (zones[2][idx2] >= block_size[idx2]) {
         stencil.blocks[0 + offset] = bidx_self;
         stencil.blocks[1 + offset] = bidx_hori;
         stencil.blocks[2 + offset] = bidx_vert;
         stencil.blocks[3 + offset] = bidx_diag;
      }
      else {
         stencil.blocks[0 + offset] = bidx_self;
         stencil.blocks[1 + offset] = bidx_hori;
         stencil.blocks[2 + offset] = bidx_self;
         stencil.blocks[3 + offset] = bidx_hori;
      };
   }

// Top and botttom margins. The corner cases have already been covered.
   else if (zones[0][idx2] < 0) {
      stencil.blocks[0 + offset] = bidx_vert;
      stencil.blocks[1 + offset] = bidx_vert;
      stencil.blocks[2 + offset] = bidx_self;
      stencil.blocks[3 + offset] = bidx_self;
   }
   else if (zones[2][idx2] >= block_size[idx2]) {
      stencil.blocks[0 + offset] = bidx_self;
      stencil.blocks[1 + offset] = bidx_self;
      stencil.blocks[2 + offset] = bidx_vert;
      stencil.blocks[3 + offset] = bidx_vert;
   }

// Central part of the block
   else {
      stencil.blocks[0 + offset] = bidx_self;
      stencil.blocks[1 + offset] = bidx_self;
      stencil.blocks[2 + offset] = bidx_self;
      stencil.blocks[3 + offset] = bidx_self;
   };

// Compute the proper zones by wrapping coordinates. Note that the out of plane coordinates are invalid and must be set in the calling function.
   stencil.zones[0 + offset] = (zones[0] + block_size) % block_size;
   stencil.zones[1 + offset] = (zones[1] + block_size) % block_size;
   stencil.zones[2 + offset] = (zones[2] + block_size) % block_size;
   stencil.zones[3 + offset] = (zones[3] + block_size) % block_size;

// Two-dimensional weights are products of opposite offsets
   stencil.weights[0 + offset] = offset_hi[idx1] * offset_hi[idx2];
   stencil.weights[1 + offset] = offset_lo[idx1] * offset_hi[idx2];
   stencil.weights[2 + offset] = offset_hi[idx1] * offset_lo[idx2];
   stencil.weights[3 + offset] = offset_lo[idx1] * offset_lo[idx2];

// Derivatives in the direction normal to the plane are equal to the weights at this point
   stencil.derivatives[3 * (0 + offset) + idx0] = stencil.weights[0 + offset];
   stencil.derivatives[3 * (1 + offset) + idx0] = stencil.weights[1 + offset];
   stencil.derivatives[3 * (2 + offset) + idx0] = stencil.weights[2 + offset];
   stencil.derivatives[3 * (3 + offset) + idx0] = stencil.weights[3 + offset];

// Derivatives in the first direction in the plane
   stencil.derivatives[3 * (0 + offset) + idx1] = -offset_hi[idx2] / delta[idx1];
   stencil.derivatives[3 * (1 + offset) + idx1] =  offset_hi[idx2] / delta[idx1];
   stencil.derivatives[3 * (2 + offset) + idx1] = -offset_lo[idx2] / delta[idx1];
   stencil.derivatives[3 * (3 + offset) + idx1] =  offset_lo[idx2] / delta[idx1];

// Derivatives in the second direction in the plane
   stencil.derivatives[3 * (0 + offset) + idx2] = -offset_hi[idx1] / delta[idx2];
   stencil.derivatives[3 * (1 + offset) + idx2] = -offset_lo[idx1] / delta[idx2];
   stencil.derivatives[3 * (2 + offset) + idx2] =  offset_hi[idx1] / delta[idx2];
   stencil.derivatives[3 * (3 + offset) + idx2] =  offset_lo[idx1] / delta[idx2];

   return 0;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 12/01/2023
\param[in] pos Interpolation point position
\return Status: 0 if local interpolation was used, 1, 2, and 3 if plane interpolation was used, 4 if external interpolation was used
*/
int ServerBATLFront::BuildInterpolationStencil(const GeoVector& pos)
{
#if REQUEST_STENCIL_FROM_BATL == 0

   int pri_idx, sec_idx, xyz, iz;
   MultiIndex zone_lo, zone_hi;
   GeoVector offset_lo, offset_hi, delta;

// Get primary block data
   pri_idx = block_pri->GetNode();
   block_pri->GetZoneOffset(pos, zone_lo, offset_lo);
   delta = block_pri->GetZoneLength();

   offset_hi = 1.0 - offset_lo;
   zone_hi = zone_lo + 1;

#if SERVER_NUM_GHOST_CELLS == 0

   int nbr_level, plane, idx0, idx1, idx2, outcome, ipl, oop_zone_pri, oop_zone_sec;
   MultiIndex block_size, quadrant, planes, node_idx[3], level_idx;
   double offset_pri, offset_sec, del;
   quadrant = block_pri->GetQuadrant(pos);
   block_size = block_pri->GetBlockSize();

// Build three neighbor node indices for up to 3 possible secondary blocks. These are face neighbors, never edge or vertex neighbors. Also, find if the stencil extends beyond the primary block.
   for (xyz = 0; xyz < 3; xyz++) {
      node_idx[xyz][xyz] = (quadrant[xyz] == 1 ? 0 : 3);
      node_idx[xyz][(xyz + 1) % 3] = quadrant[(xyz + 1) % 3];
      node_idx[xyz][(xyz + 2) % 3] = quadrant[(xyz + 2) % 3];
      level_idx[xyz] = (zone_lo[xyz] < 0 ? 0 : (zone_hi[xyz] >= block_size[xyz] ? 2 : 1));
   };

// Interior positions - local interpolation in the primary block.
   if ((level_idx.i == 1) && (level_idx.j == 1) && (level_idx.k == 1)) {
      for (iz = 0; iz < 8; iz++) stencil.blocks[iz] = pri_idx;
      InteriorInterpolationStencil(zone_lo, zone_hi, offset_lo, offset_hi, delta);
      return 0;
   };

// Figure out up to three possible planes (xyz, xy, xz, x, yz, y, z)
   for (xyz = 0; xyz < 3; xyz++) planes[xyz] = -1;

   if (level_idx.i != 1) {
      planes[0] = 0;
      if (level_idx.j != 1) {
         planes[1] = 1;
         if (level_idx.k != 1) planes[2] = 2;
      }
      else if (level_idx.k != 1) planes[1] = 2;
   }
   else if (level_idx.j != 1) {
      planes[0] = 1;
      if (level_idx.k != 1) planes[1] = 2;
   }
   else if (level_idx.k != 1) planes[0] = 2;
   else {
      std::cerr << "Interpolation plane error\n";
      return -1;
   };

// Try up to three planes
// TODO Do we need to request an entire block for each "ipl"? Two out of 3 might be wrong, and RequestBlock() is expensive.
   for (ipl = 0; ipl < 3; ipl++) {
      if (planes[ipl] == -1) break;
      plane = planes[ipl];

// First see if we can interpolate in the plane of the primary, go directly to the next plane on failure.
      stencil.blocks[0] = pri_idx;
      outcome = BuildInterpolationPlane(pos, plane, 0);
      if (outcome == -1) continue;

// Find the secondary block
      _inquiry.type = 0;
      _inquiry.node = block_pri->GetNeighborNode(node_idx[plane]);
      sec_idx = RequestBlock();

// See if we can interpolate in the plane of the secondary. If successful, terminate the loop. The value of "plane" can be used later.
      stencil.blocks[4] = sec_idx;
      outcome = BuildInterpolationPlane(pos, plane, 1);
      if (outcome == -1) continue;
      else break;
   };

// The plane interpolator has failed (level change in all three planes). We must request a stencil using the external interpolator.
   if (outcome == -1) return RequestStencil(pos);

   idx0 = plane;
   idx1 = (plane + 1) % 3;
   idx2 = (plane + 2) % 3;

// Use "node_idx" to figure out the position and level of the neighbor
   level_idx = NodeToLevel(node_idx[idx0]);
   nbr_level = block_pri->GetNeighborLevel(level_idx);

// Secondary block is below/above
   if (node_idx[idx0][idx0] == 0) {
      offset_pri = offset_hi[idx0];
      del = delta[idx0];
   }
   else {
      offset_pri = offset_lo[idx0];
      del = -delta[idx0];
   };

// Secondary block is coarser/finer/same
   if (nbr_level == -1) {
      offset_pri *= 4.0 / 3.0;
      del *= 1.5;
   }
   else if (nbr_level == 1) {
      offset_pri *= 2.0 / 3.0;
      del *= 0.75;
   };
   offset_sec = 1.0 - offset_pri;

// Correct the out of plane zone components and scale the weight for each plane.
   oop_zone_pri = (node_idx[idx0][idx0] == 0 ? 0 : block_size[idx0] - 1);
   oop_zone_sec = block_size[idx0] - oop_zone_pri - 1;
   for (iz = 0; iz < 4; iz++) {
      stencil.zones[iz][idx0] = oop_zone_pri;
      stencil.weights[iz] *= offset_sec;
      stencil.derivatives[3 * iz + idx0] /= del;
      stencil.derivatives[3 * iz + idx1] *= offset_sec;
      stencil.derivatives[3 * iz + idx2] *= offset_sec;

      stencil.zones[iz + 4][idx0] = oop_zone_sec;
      stencil.weights[iz + 4] *= offset_pri;
      stencil.derivatives[3 * iz + idx0] /= del;
      stencil.derivatives[3 * iz + idx1] *= offset_pri;
      stencil.derivatives[3 * iz + idx2] *= offset_pri;
   };

   stencil.n_elements = 8;
   return (nbr_level == 0 ? 1 : (nbr_level == -1 ? 2 : 3));

#elif SERVER_NUM_GHOST_CELLS > 0

// Interpolation is always internal
   InteriorInterpolationStencil(zone_lo, zone_hi, offset_lo, offset_hi, delta);
   return 0;

#endif

#else

// Load stencil directly from BATL
   return RequestStencil(pos);

#endif
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATLBack public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/25/2022
\param[in] file_name_pattern_in A string describing the file naming pattern
*/
ServerBATLBack::ServerBATLBack(const std::string& file_name_pattern_in)
              : ServerCartesianBack(file_name_pattern_in),
                ServerBaseBack(file_name_pattern_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/01/2023
*/
void ServerBATLBack::ReadData(const std::string data_file)
{
// Initialize the BATL library and read the data into memory
   MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_SELF);
   spectrum_init_mpi(fcomm);

   wrapamr_read_header(data_file.c_str(), line_width, 1);
   int vars_inds_fortran[n_variables];
   for (int i = 0; i < n_variables; i++) vars_inds_fortran[i] = i+1;
   wrapamr_read_file_partial(data_file.c_str(), line_width, 0, 0, n_variables, vars_inds_fortran);
   wrapamr_get_domain(domain_min.Data(), domain_max.Data());
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
*/
void ServerBATLBack::CleanReader(void)
{
   wrapamr_clean();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
\return Number of clients that completed their tasks during this cycle

\note "needvars" requests are only handled when SERVER_INTERP_ORDER = -1. "needstencil" requests are only handled when using 1st order interpolation and either there are no ghost cells or stencils are always requested from BATL (for debugging purposes)
*/
int ServerBATLBack::ServerFunctions(void)
{
#if SERVER_INTERP_ORDER == -1
// Handle "needvars" requests
   HandleNeedVarsRequests();
#elif SERVER_INTERP_ORDER == 1 && (SERVER_NUM_GHOST_CELLS == 0 || REQUEST_STENCIL_FROM_BATL == 1)
// Handle "needstencil" requests
   HandleNeedStencilRequests();
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
void ServerBATLBack::GetBlockData(const double* pos, double* vars, int* found)
{
   wrapamr_get_data_serial(pos, vars, found);
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
\param[in] pos    position array in reader coordinates
\param[out] node  ID tag of block containing pos within its physical boundaries
*/
void ServerBATLBack::GetBlock(const double* pos, int* node)
{
   spectrum_get_node(pos, node);
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
*/
void ServerBATLBack::HandleNeedStencilRequests(void)
{
   GeoVector pos_batl;
   int cpu, cpu_idx, count_needstencil = 0;

// Service the "needstencil" requests
   MPI_Testsome(MPI_Config::node_comm_size, req_needstencil, &count_needstencil, index_needstencil, MPI_STATUSES_IGNORE);

   for (cpu_idx = 0; cpu_idx < count_needstencil; cpu_idx++) {
      cpu = index_needstencil[cpu_idx];

// Obtain the stencil requested.
      pos_batl = buf_needstencil[cpu].pos / unit_length_server * unit_length_fluid;
      spectrum_get_interpolation_stencil(pos_batl.Data(), &stencil.n_elements, stencil.blocks, &stencil.zones[0][0], stencil.weights);

// Send the stencil to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(&stencil, 1, MPIStencilType, cpu, tag_sendstencil, MPI_Config::node_comm);

// Post the receive for the next stencil request from this worker.
      MPI_Irecv(&buf_needstencil[cpu], 1, MPIInquiryType, cpu, tag_needstencil, MPI_Config::node_comm, &req_needstencil[cpu]);
   };
};

};
