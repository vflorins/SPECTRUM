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
\date 01/04/2024
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::Start(void) {

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
   block_pri->SetDimensions(domain_max, domain_min);
   LoadFromReader(block_pri);
   block_pri->LoadDimensions(1.0);
   block_sec = std::make_shared<Block>();
   block_sec->SetDimensions(domain_max, domain_min);
   LoadFromReader(block_sec);
   block_sec->LoadDimensions(1.0);
   block_stn = std::make_shared<Block>();
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/22/2025
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::Finish(void)
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
      block_new = std::make_shared<Block>();

// Receive the block in 4 parts (member data plus 3 dynamic arrays).
// This is called even if SERVER_INTERP_ORDER is -1 to import the block dimensions
      MPI_Recv(block_new.get(), 1, MPIBlockType, 0, MPI::tag::sendblock, MPI::node_comm, MPI_STATUS_IGNORE);

      if constexpr (server_interp_order > -1) {
         MPI_Recv(block_new->GetVariablesAddress(), block_new->GetVariableCount() * block_new->GetZoneCount(), MPI_DOUBLE, 0,
                  MPI::tag::sendblock, MPI::node_comm, MPI_STATUS_IGNORE);
      }
      if constexpr (server_interp_order > 0 && num_ghost_cells == 0) {
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
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/28/2025
\param[out] fields Fields
*/
template <typename HConfig>
template <typename Fields, typename RequestedFields>
status_t BackgroundDataCartesian<HConfig>::Evaluate_FromReader(Fields& fields)
{
   DataFields datafields;
// Get variables
   MPI_Send(&_inquiry, 1, MPIInquiryType, 0, MPI::tag::needvars, MPI::node_comm);
   MPI_Recv(datafields.Array(), DataFields::size(), MPI_DOUBLE, 0, MPI::tag::sendvars, MPI::node_comm, MPI_STATUS_IGNORE);
   stencil_outcomes[2]++;

// Mass density
   if constexpr (RequestedFields::MassDen_found()) {
      fields.MassDen('w') = datafields.MassDen();
   }
// Number density
   if constexpr (RequestedFields::MassDen_found()) {
      fields.Den('w') = datafields.Den();
   }
// Thermal pressure
   if constexpr (RequestedFields::Prs_found()) {
      fields.Prs('w') = datafields.Prs();
   }
// Bulk flow from mass density and momentum
   if constexpr (RequestedFields::Fluv_found()) {
      if constexpr (DataFields::Mom_found() && DataFields::MassDen_found())
         fields.Fluv('w') = datafields.Mom()/datafields.MassDen();
      else if constexpr (DataFields::Fluv_found())
         fields.Fluv('w') = datafields.Fluv();
      else
         fields.Fluv('w') = gv_zeros;
   }
// Magnetic field
   if constexpr (RequestedFields::Mag_found()) {
         fields.Mag('w') = datafields.Mag();
   }
// Electric field
   if constexpr (RequestedFields::Elc_found()) {
      if constexpr (DataFields::Fluv_found() && !(DataFields::Flum_found() && DataFields::MassDen_found()))
         fields.Elc('w') = gv_zeros;
      else if constexpr (DataFields::Ele_found())
         fields.Ele('w') = datafields.Ele();
      else if constexpr (DataFields::Fluv_found() && DataFields::Mag_found())
         fields.Ele('w') = -(fields.Fluv() ^ fields.Mag()) / c_code;
   }
   if constexpr (RequestedFields::Iv0_found()) {
      fields.Iv0('w') = datafields.Iv0();
   }
   if constexpr (RequestedFields::Iv1_found()) {
      fields.Iv1('w') = datafields.Iv1();
   }
   if constexpr (RequestedFields::Iv2_found()) {
      fields.Iv2('w') = datafields.Iv2();
   }
   return 0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/28/2025
\param[in]  pos    Position
\param[out] fields Fields
*/
template <typename HConfig>
template <typename Fields, typename RequestedFields>
status_t BackgroundDataCartesian<HConfig>::Evaluate_Interp0(const GeoVector& pos, Fields& fields)
{
// Neatens expressions involving the block we will be using
#define B (*block_pri)
// Take the nearest cell value
   MultiIndex zone = B.GetZone(pos);

   // mass density
   if constexpr (RequestedFields::MassDen_found()) {
      fields.MassDen('w') = B[zone].Rho();
   }
   // number density
   if constexpr (RequestedFields::Den_found()) {
      fields.Den('w') = B[zone].Den();
   }
   // thermal pressure
   if constexpr (RequestedFields::Prs_found()) {
      fields.Prs('w') = B[zone].Prs();
   }
   // bulk flow velocity
   if constexpr (RequestedFields::Fluv_found()) {
      if constexpr (DataFields::Fluv_found())
         fields.Fluv('w') = B[zone].Fluv();
      else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found())
         fields.Fluv('w') = B[zone].Mom()/B[zone].MassDen();
      else
         fields.Fluv('w') = gv_zeros;
   }
   if constexpr (RequestedFields::Mag_found()) {
      fields.Mag('w') = B[zone].Mag();
   }
// Electric field
   if constexpr (RequestedFields::Ele_found()) {
      if constexpr (DataFields::Ele_found())
         fields.Elc('w') = B[zone].Elc();
      else if constexpr (DataFields::Fluv_found() && DataFields::Mag_found())
         fields.Ele('w') = -(fields.Fluv() ^ fields.Mag()) / c_code;
      else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found())
         fields.Ele('w') = -((fields.Mom()/fields.MassDen()) ^ fields.Mag()) / c_code;
      else
         fields.Ele('w') = gv_zeros;
   }
// Region indicator variables
   if constexpr (Fields::Iv0_found()) {
      if (DataFields::Iv0_found())
         fields.Iv0('w') = B[zone].Iv0();
      else
         fields.Iv0('w') = 0.0;
   }
   if constexpr (Fields::Iv1_found()) {
      if (DataFields::Iv1_found())
         fields.Iv1('w') = B[zone].Iv1();
      else
         fields.Iv1('w') = 0.0;
   }
   if constexpr (Fields::Iv2_found()) {
      if (DataFields::Iv2_found())
         fields.Iv2('w') = B[zone].Iv2();
      else
         fields.Iv2('w') = 0.0;
   }
#undef B
};


/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/28/2025
\param[in]  pos    Position
\param[out] spdata Fields
*/
// TODO: the same code block is used three times
template <typename HConfig>
template <typename Fields, typename RequestedFields>
status_t BackgroundDataCartesian<HConfig>::Evaluate_Interp1(const GeoVector& pos, Fields& fields)
{
// Neatens expressions involving the block we will be using
   int iz;
   MultiIndex zone;
   int pri_idx, sec_idx;

// Build the stencil. This is a time consuming operation.
   stencil_status = BuildInterpolationStencil(pos);
// Reconstruction as weighted stencil values.
// If AbsMag is a stored value,
// AbsMag is reconstructed from the AbsMag stencil
// and not from the reconstructed Mag.
// AbsMag does not need to be a stored value of DataFields.
   if constexpr (RequestedFields::MassDen_found()) {
      fields.MassDen('w') = 0.0;
   }
   if constexpr (RequestedFields::Den_found()) {
      fields.Den('w') = 0.0;
   }
   if constexpr (RequestedFields::Prs_found()) {
      fields.Prs('w') = 0.0;
   }
   if constexpr (RequestedFields::Fluv_found()) {
      fields.Fluv('w') = gv_zeros;
   }
   if constexpr (RequestedFields::Mag_found()) {
      fields.Mag('w') = gv_zeros;
   }
   if constexpr (RequestedFields::AbsMag_found()) {
      fields.AbsMag('w') = 0.0;
   }
   if constexpr (RequestedFields::Ele_found()) {
      fields.Ele('w') = gv_zeros;
   }
   if constexpr (RequestedFields::Iv0_found()) {
      fields.Iv0('w') = 0.0;
   }
   if constexpr (RequestedFields::Iv1_found()) {
      fields.Iv1('w') = 0.0;
   }
   if constexpr (RequestedFields::Iv2_found()) {
      fields.Iv2('w') = 0.0;
   }
// Internal interpolation
   if (stencil_status == 0) {
      stencil_outcomes[0]++;
      for (iz = 0; iz < stencil.n_elements; iz++) {
#define B (*block_pri)
         zone = stencil.zones[iz];
         if constexpr (RequestedFields::MassDen_found()) {
            fields.MassDen('w') += stencil.weights[iz] * B[zone].MassDen();
         }
         if constexpr (RequestedFields::Den_found()) {
            fields.Den('w') += stencil.weights[iz] * B[zone].Den();
         }
         if constexpr (RequestedFields::Prs_found()) {
            fields.Prs('w') += stencil.weights[iz] * B[zone].Prs();
         }
         if constexpr (RequestedFields::Fluv_found()) {
            if constexpr (DataFields::Fluv_found())
               fields.Fluv('w') += stencil.weights[iz] * B[zone].Fluv();
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found())
               fields.Fluv('w') = stencil.weights[iz] * B[zone].Mom()/B[zone].MassDen();
            else
               fields.Fluv('w') = gv_zeros;
         }
         if constexpr (RequestedFields::Mag_found()) {
            fields.Mag('w') += stencil.weights[iz] * B[zone].Mag();
         }
         if constexpr (RequestedFields::AbsMag_found()) {
            fields.AbsMag('w') += stencil.weights[iz] * B[zone].AbsMag();
         }
         if constexpr (RequestedFields::Ele_found()) {
            if constexpr (DataFields::Ele_found())
               fields.Ele('w') += stencil.weights[iz] * B[zone].Ele();
            else if constexpr (DataFields::Fluv_found() && DataFields::Mag_found())
               fields.Ele('w') = stencil.weights[iz] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found())
               fields.Ele('w') = stencil.weights[iz] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            else
               fields.Ele('w') = gv_zeros;
         }
         if constexpr (RequestedFields::Iv0_found()) {
            fields.Iv0('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv0();
         }
         if constexpr (RequestedFields::Iv1_found()) {
            fields.Iv1('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv1();
         }
         if constexpr (RequestedFields::Iv2_found()) {
            fields.Iv2('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv2();
         }
#undef B
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
#define B (*block_stn)
         zone = stencil.zones[iz];
         if constexpr (RequestedFields::MassDen_found()) {
            fields.MassDen('w') += stencil.weights[iz] * B[zone].MassDen();
         }
         if constexpr (RequestedFields::Den_found()) {
            fields.Den('w') += stencil.weights[iz] * B[zone].Den();
         }
         if constexpr (RequestedFields::Prs_found()) {
            fields.Prs('w') += stencil.weights[iz] * B[zone].Prs();
         }
         if constexpr (RequestedFields::Fluv_found()) {
            if constexpr (DataFields::Fluv_found())
               fields.Fluv('w') += stencil.weights[iz] * B[zone].Fluv();
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found())
               fields.Fluv('w') = stencil.weights[iz] * B[zone].Mom()/B[zone].MassDen();
            else
               fields.Fluv('w') = gv_zeros;
         }
         if constexpr (RequestedFields::Mag_found()) {
            fields.Mag('w') += stencil.weights[iz] * B[zone].Mag();
         }
         if constexpr (RequestedFields::AbsMag_found()) {
            fields.AbsMag('w') += stencil.weights[iz] * B[zone].AbsMag();
         }
         if constexpr (RequestedFields::Ele_found()) {
            if constexpr (DataFields::Ele_found())
               fields.Ele('w') += stencil.weights[iz] * B[zone].Ele();
            else if constexpr (DataFields::Fluv_found() && DataFields::Mag_found())
               fields.Ele('w') = stencil.weights[iz] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found())
               fields.Ele('w') = stencil.weights[iz] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            else
               fields.Ele('w') = gv_zeros;
         }
         if constexpr (RequestedFields::Iv0_found()) {
            fields.Iv0('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv0();
         }
         if constexpr (RequestedFields::Iv1_found()) {
            fields.Iv1('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv1();
         }
         if constexpr (RequestedFields::Iv2_found()) {
            fields.Iv2('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv2();
         }
#undef B
      };
   }

// Plane interpolation
   else {
      stencil_outcomes[1]++;
      pri_idx = block_pri->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else block_stn = block_sec;
#define B (*block_stn)
         zone = stencil.zones[iz];
         if constexpr (RequestedFields::MassDen_found()) {
            fields.MassDen('w') += stencil.weights[iz] * B[zone].MassDen();
         }
         if constexpr (RequestedFields::Den_found()) {
            fields.Den('w') += stencil.weights[iz] * B[zone].Den();
         }
         if constexpr (RequestedFields::Prs_found()) {
            fields.Prs('w') += stencil.weights[iz] * B[zone].Prs();
         }
         if constexpr (RequestedFields::Fluv_found()) {
            if constexpr (DataFields::Fluv_found())
               fields.Fluv('w') += stencil.weights[iz] * B[zone].Fluv();
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found())
               fields.Fluv('w') = stencil.weights[iz] * B[zone].Mom()/B[zone].MassDen();
            else
               fields.Fluv('w') = gv_zeros;
         }
         if constexpr (RequestedFields::Mag_found()) {
            fields.Mag('w') += stencil.weights[iz] * B[zone].Mag();
         }
         if constexpr (RequestedFields::AbsMag_found()) {
            fields.AbsMag('w') += stencil.weights[iz] * B[zone].AbsMag();
         }
         if constexpr (RequestedFields::Ele_found()) {
            if constexpr (DataFields::Ele_found())
               fields.Ele('w') += stencil.weights[iz] * B[zone].Ele();
            else if constexpr (DataFields::Fluv_found() && DataFields::Mag_found())
               fields.Ele('w') = stencil.weights[iz] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found())
               fields.Ele('w') = stencil.weights[iz] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            else
               fields.Ele('w') = gv_zeros;
         }
         if constexpr (RequestedFields::Iv0_found()) {
            fields.Iv0('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv0();
         }
         if constexpr (RequestedFields::Iv1_found()) {
            fields.Iv1('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv1();
         }
         if constexpr (RequestedFields::Iv2_found()) {
            fields.Iv2('w') += stencil.weights[iz] * B[stencil.zones[iz]].Iv2();
         }
#undef B
      };
   };
};




/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/28/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundDataCartesian<HConfig>::Gradients_Interp1(Coordinates& coords, Fields& fields)
{
   double var, Bmag2, rho = 0.0;
   int vidx, xyz, uvw, iz, pri_idx, sec_idx;

   double grads[Fields::size()][3] = {0.0};
   fields.DelAbsMag() = gv_zeros;

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
         Bmag2 = Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) fields.DelMag()[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(Bmag2);
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
         Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) fields.DelAbsMag()[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(Bmag2);
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
         Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                  + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) fields.DelAbsMag()[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(Bmag2);
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
// Perform unit conversion
   fields.DelVel() *= unit_velocity_server / unit_velocity_fluid;
   fields.DelMag() *= unit_magnetic_server / unit_magnetic_fluid;
   fields.DelElc() *= unit_electric_server / unit_electric_fluid;

   return 0;
};











/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/28/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundDataCartesian<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   // todo static assert every field in RequestedFields is also in DataFields

   int bidx;
   status_t status;

// Request the block and the zone size
// todo awk - _inquiry type not well-typed [sic]
   _inquiry.type = 1;
   _inquiry.pos = coords.Pos();
   bidx = RequestBlock();

// If "block_pri" or "block_sec" is the position owner (based on the call to RequestBlock), we don't need to acccess the cache
   if (block_pri->GetNode() != bidx) {
      if (block_sec->GetNode() == bidx) block_pri = block_sec;
      else block_pri = cache_line[bidx];
   };

   if constexpr (server_interp_order == -1) {
// Get fields directly from reader program
      status = Evaluate_FromReader<Fields, RequestedFields>(fields);
   }
   else if constexpr (server_interp_order == 0) {
// Get fields using 0th order interpolation
      status = Evaluate_Interp0<Fields, RequestedFields>(coords.Pos(), fields);
   }
   else if constexpr (server_interp_order == 1) {
// Get fields using 1st order interpolation
      status = Evaluate_Interp1<Fields, RequestedFields>(coords.Pos(), fields);
   }

   if constexpr (RequestedFields::Prs_found()) {
      fields.Prs() *= unit_pressure_server / unit_pressure_fluid;
   }
   if constexpr (RequestedFields::Den_found()) {
      fields.Den() *= unit_number_density_server / unit_number_density_fluid;
      if constexpr (RequestedFields::Iv0_found())
         fields.Iv0() /= fields.Den();
      if constexpr (RequestedFields::Iv1_found())
         fields.Iv1() /= fields.Den();
      if constexpr (RequestedFields::Iv2_found())
         fields.Iv2() /= fields.Den();
      // ...adjust region modifications as needed
   }
   if constexpr (RequestedFields::Fluv_found()) {
      fields.Fluv() *= unit_velocity_server / unit_velocity_fluid;
   }
   if constexpr (RequestedFields::Mag_found()) {
      fields.Mag() *= unit_magnetic_server / unit_magnetic_fluid;
   }
   if constexpr (RequestedFields::Elc_found()) {
      fields.Elc() *= unit_electric_server / unit_electric_fluid;
   }
   return status;
};







/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/27/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundDataCartesian<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   if constexpr (server_interp_order == -1) {
// Gradients must be computed numerically
      return FALLBACK_NUMERIC_DERIVATIVES;
   }
   else if constexpr (server_interp_order == 0) {
// All gradients are explicitly set to zero, and the background must not attempt to compute them using "NumericalDerivatives()"
      if constexpr (RequestedFields::DelFluv_found()) {
         fields.DelFluv() = gm_zeros;
      }
      if constexpr (RequestedFields::DelMag_found()) {
         fields.DelMag() = gm_zeros;
      }
      if constexpr (RequestedFields::DelAbsMag_found()) {
         fields.DelAbsMag() = gv_zeros;
      }
      if constexpr (RequestedFields::DelElc_found()) {
         fields.DelElc() = gm_zeros;
      }
   }
   else if constexpr (server_interp_order == 1) {
// Gradients can be obtained from the stencil construction
      Gradients_Interp1<Coordinates, Fields, RequestedFields>(coords, fields);
   }
   return 0;
};



/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/27/2025
*/
template <typename HConfig>
template <typename Coordinates>
status_t BackgroundDataCartesian<HConfig>::EvaluateDmax(Coordinates& coords, double& dmax)
{
   dmax = fmin(dmax0, block_pri->GetZoneLength().Smallest());
   return 0;
};



/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::PrintStencilOutcomes(void) const requires (HConfig::buildmode == BuildMode::debug)
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
void BackgroundDataCartesian<HConfig>::PrintNumBlocksRequested(void) const requires (HConfig::buildmode == BuildMode::debug)
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
int BackgroundDataCartesian<HConfig>::GetNCachedBlocks(void) const requires (HConfig::buildmode == BuildMode::debug)
{
   return cache_line.size();
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/27/2025
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::InvalidateCache(void) requires (HConfig::buildmode == BuildMode::debug)
{
   cache_line.Empty();
};


};
