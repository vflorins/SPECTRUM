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
\date 11/28/2025
*/
template <typename HConfig>
void BackgroundDataCartesian<HConfig>::Start(void) {
   BackgroundDataBase::Start();
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
   BackgroundDataBase::Finish();
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
   constexpr auto block_size = Block::block_size;
   MultiIndex node_idx, zone_lo, zone_hi;
   GeoVector offset_lo, offset_hi, delta;

// Get primary block data
   pri_idx = block_pri->GetNode();
   block_pri->GetZoneOffset(pos, zone_lo, offset_lo);
   delta = block_pri->GetZoneLength();

   offset_hi = 1.0 - offset_lo;
   zone_hi = zone_lo + 1;

// Default to local interpolation in the primary block
   BackgroundDataBase::InteriorInterpolationStencil(zone_lo, zone_hi, offset_lo, offset_hi, delta);

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
            stencil.blocks[iz] = BackgroundDataBase::RequestBlock();
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
// TODO: return value missing, review
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
// AbsMag does not need to be a stored value of DataFields;
// if it isn't, it will be computed by the Fields type.
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
               fields.Fluv('w') += stencil.weights[iz] * B[zone].Mom()/B[zone].MassDen();
            else
               ;
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
               fields.Ele('w') += stencil.weights[iz] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found())
               fields.Ele('w') += stencil.weights[iz] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            else
               ;
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
               fields.Fluv('w') += stencil.weights[iz] * B[zone].Mom()/B[zone].MassDen();
            else
               ;
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
               fields.Ele('w') += stencil.weights[iz] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found())
               fields.Ele('w') += stencil.weights[iz] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            else
               ;
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
               fields.Fluv('w') += stencil.weights[iz] * B[zone].Mom()/B[zone].MassDen();
            else
               ;
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
               fields.Ele('w') += stencil.weights[iz] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found())
               fields.Ele('w') += stencil.weights[iz] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            else
               ;
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
// TODO: the same code block is used three times
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundDataCartesian<HConfig>::Gradients_Interp1(Coordinates& coords, Fields& fields)
{
   int iz;
   MultiIndex zone;
   int uvw;
   int pri_idx, sec_idx;

// todo in the Spectrum-v1 code, rho (mass density) is computed here, 'on-the-fly'.
//  But there should be a code invariant, namely, quantities in Fields that are not
//  derivatives should be up-to-date at the time of this call. So rho should be
//  already available as fields.MassDen(). Check or ask around.

// todo I tried, but cannot convince myself that the gradients are computed
//  correctly, but I attempted to follow the code as I received it in spite of this.

// todo I cannot guarantee that the formatting is handled correctly,
//  e.g. matrices may need to be transposed.
//  If so, you can use a GeoMatrix call at the end.

   if constexpr (RequestedFields::DelMassDen_found()) {
      fields.DelMassDen('w') = gv_zeros;
   }
   if constexpr (RequestedFields::DelDen_found()) {
      fields.DelDen('w') = gv_zeros;
   }
   if constexpr (RequestedFields::DelPrs_found()) {
      fields.DelPrs('w') = gv_zeros;
   }
   if constexpr (RequestedFields::DelFluv_found()) {
      fields.DelFluv('w') = gm_zeros;
   }
   if constexpr (RequestedFields::DelMag_found()) {
      fields.DelMag('w') = gm_zeros;
   }
   if constexpr (RequestedFields::DelAbsMag_found()) {
      fields.DelAbsMag('w') = gv_zeros;
   }
   if constexpr (RequestedFields::DelEle_found()) {
      fields.DelEle('w') = gm_zeros;
   }

// Internal interpolation
   if (stencil_status == 0) {
      for (iz = 0; iz < stencil.n_elements; iz++) {
#define B (*block_pri)
         zone = stencil.zone[iz];
         // todo one big uvw loop , maybe ... it would be less intuitive
         if constexpr (RequestedFields::DelMassDen_found()) {
            if constexpr (DataFields::MassDen_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelMassDen('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].MassDen()[uvw];
            }
         }
         if constexpr (RequestedFields::DelDen_found()) {
            if constexpr (DataFields::Den_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelDen('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Den()[uvw];
            }
         }
         if constexpr (RequestedFields::DelPrs_found()) {
            if constexpr (DataFields::Prs_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelPrs('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Prs()[uvw];
            }
         }
         if constexpr (RequestedFields::DelAbsMag_found()) {
            if constexpr (DataFields::AbsMag_found() || DataFields::Mag_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelAbsMag('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].AbsMag()[uvw];
            }
         }
         if constexpr (RequestedFields::DelFluv_found()) {
            if constexpr (DataFields::Fluv_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelFluv('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Fluv()[uvw];
            }
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found()){
               fields.DelFluv('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Mom()[uvw]/B[zone].MassDen();
            }
            else
               ;
            fields.DelFluv('w') *= unit_velocity_server / unit_velocity_fluid;
         }
         if constexpr (RequestedFields::DelMag_found()) {
            if constexpr (DataFields::Mag_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelMag('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Mag()[uvw];
            }
            fields.DelMag('w') *= unit_magnetic_server / unit_magnetic_fluid;
         }
         if constexpr (RequestedFields::DelEle_found()) {
            if constexpr (DataFields::Ele_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelEle('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Ele()[uvw];
            }
            else if constexpr (DataFields::Fluv_found() && DataFields::Mag_found()) {
               fields.DelEle('w')[uvw] += stencil.derivatives[3 * iz + uvw] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            }
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found()) {
               fields.DelEle('w')[uvw] += stencil.derivatives[3 * iz + uvw] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            }
            else
               ;
            fields.DelEle('w') *= unit_electric_server / unit_electric_fluid;
         }
#undef B
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
#define B (*block_stn)
         zone = stencil.zone[iz];
         // todo one big uvw loop , maybe ... it would be less intuitive
         if constexpr (RequestedFields::DelMassDen_found()) {
            if constexpr (DataFields::MassDen_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelMassDen('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].MassDen()[uvw];
            }
         }
         if constexpr (RequestedFields::DelDen_found()) {
            if constexpr (DataFields::Den_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelDen('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Den()[uvw];
            }
         }
         if constexpr (RequestedFields::DelPrs_found()) {
            if constexpr (DataFields::Prs_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelPrs('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Prs()[uvw];
            }
         }
         if constexpr (RequestedFields::DelAbsMag_found()) {
            if constexpr (DataFields::AbsMag_found() || DataFields::Mag_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelAbsMag('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].AbsMag()[uvw];
            }
         }
         if constexpr (RequestedFields::DelFluv_found()) {
            if constexpr (DataFields::Fluv_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelFluv('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Fluv()[uvw];
            }
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found()){
               fields.DelFluv('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Mom()[uvw]/B[zone].MassDen();
            }
            else
               ;
            fields.DelFluv('w') *= unit_velocity_server / unit_velocity_fluid;
         }
         if constexpr (RequestedFields::DelMag_found()) {
            if constexpr (DataFields::Mag_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelMag('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Mag()[uvw];
            }
            fields.DelMag('w') *= unit_magnetic_server / unit_magnetic_fluid;
         }
         if constexpr (RequestedFields::DelEle_found()) {
            if constexpr (DataFields::Ele_found()) {
               for (uvw = 0; uvw < 3; ++uvw)
                  fields.DelEle('w')[uvw] += stencil.derivatives[3 * iz + uvw] * B[zone].Ele()[uvw];
            }
            else if constexpr (DataFields::Fluv_found() && DataFields::Mag_found()) {
               fields.DelEle('w')[uvw] += stencil.derivatives[3 * iz + uvw] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            }
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found()) {
               fields.DelEle('w')[uvw] += stencil.derivatives[3 * iz + uvw] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            }
            else
               ;
            fields.DelEle('w') *= unit_electric_server / unit_electric_fluid;
         }
#undef B
      };

   }

// Plane interpolation
   else {
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
               fields.Fluv('w') += stencil.weights[iz] * B[zone].Mom()/B[zone].MassDen();
            else
               ;
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
               fields.Ele('w') += stencil.weights[iz] * -(B[zone].Fluv() ^ B[zone].Mag()) / c_code;
            else if constexpr (DataFields::MassDen_found() && DataFields::Mom_found() && DataFields::Mag_found())
               fields.Ele('w') += stencil.weights[iz] * -((B[zone].Mom()/B[zone].MassDen()) ^ B[zone].Mag()) / c_code;
            else
               ;
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
   bidx = BackgroundDataBase::RequestBlock();

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
// Note: if server_interp_order == -1, this method is never called.
   if constexpr (server_interp_order == 0) {
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
status_t BackgroundDataCartesian<HConfig>::EvaluateDmax(Coordinates& coords, double* dmax)
{
   *dmax = fmin(dmax0, block_pri->GetZoneLength().Smallest());
   return 0;
};


};
