/*!
\file server_types.cc
\brief Defines data structures for server-worker communication in a 3D grid computational environment
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_types.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RequestInfo methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Juan G Alonso Guzman
\date 04/27/2021
\param[in] size Size of communicator in which non-blocking receives will be called
*/
RequestInfo::RequestInfo(int size)
{
   mpi_req = new MPI_Request[size];
   cpu_rank = new int[size];
};

/*!
\author Juan G Alonso Guzman
\date 04/27/2021
*/
RequestInfo::~RequestInfo()
{
   delete[] mpi_req;
   delete[] cpu_rank;
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Inquiry methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Stencil methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <int n_elements>
void Stencil<n_elements>::Print(void) const
{
   for (auto iz = 0; iz < n_elements; iz++) {
      std::cerr << std::setw(10) << blocks[iz] << "  " << zones[iz] << "  " << weights[iz] << std::endl;
   };
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Block methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/21/2025
*/
template <typename HConfig>
Block<HConfig>::Block(void)
{
   fields = new DataFields[block_size.Prod()];
   neighbor_nodes = new int[max_neighbors];
   neighbor_levels = new int[max_neighbor_levels];
};

/*!
\author Vladimir Florinski
\date 06/23/2020
\param[in] other Block to copy
*/
template <typename HConfig>
Block<HConfig>::Block(const Block<HConfig>& other)
{
   operator =(other);
};

/*!
\author Vladimir Florinski
\date 01/12/2023
\param[in] other Block to move
*/
template <typename HConfig>
Block<HConfig>::Block(Block<HConfig>&& other) noexcept
{
   operator =(other);
};


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/21/2025
*/
template <typename HConfig>
Block<HConfig>::~Block(void)
{
   delete[] fields;
   delete[] neighbor_levels;
   delete[] neighbor_nodes;
};

/*!
\author Vladimir Florinski
\date 10/26/2022
\param[in] other Block to copy
*/
template <typename HConfig>
Block<HConfig>& Block<HConfig>::operator =(const Block& other)
{
   if (&other == this) return *this;

   node = other.node;
   face_min = other.face_min;
   face_max = other.face_max;
   ConfigureProperties();

   memcpy(neighbor_nodes, other.neighbor_nodes, max_neighbors * SZINT);
   memcpy(neighbor_levels, other.neighbor_levels, max_neighbor_levels * SZINT);
   memcpy(fields, other.fields, block_size.Prod() * DataFields::size() * SZDBL);

   return *this;
};

/*!
\author Vladimir Florinski
\date 01/12/2023
\param[in] other Block to move
*/
template <typename HConfig>
Block<HConfig>& Block<HConfig>::operator =(Block<HConfig>&& other)
{
   if (&other == this) return *this;

   node = other.node;
   face_min = other.face_min;
   face_max = other.face_max;
   ConfigureProperties();

   delete[] neighbor_levels;
   delete[] neighbor_nodes;
   delete[] fields;

   neighbor_levels = other.neighbor_levels;
   neighbor_nodes = other.neighbor_nodes;
   fields = other.fields;

   other.neighbor_levels = nullptr;
   other.neighbor_nodes = nullptr;
   other.fields = nullptr;

   return *this;
};




/*!
\author Vladimir Florinski
\date 08/13/2020
\param[in] pos Coordinates of a point
\return Index of the zone that owns "pos"
*/
template <typename HConfig>
MultiIndex Block<HConfig>::GetZone(const GeoVector& pos) const
{
// This is the simplest implementation for cuboid blocks.
   MultiIndex zone = (pos - face_min) / zone_length;
   return zone;
};



/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] pos Coordinates of a point
\return Which quadrant of the block the point is in (1 or 2)
*/
template <typename HConfig>
MultiIndex Block<HConfig>::GetQuadrant(const GeoVector& pos) const
{
   MultiIndex quadrant;
   for(auto xyz = 0; xyz < 3; xyz++) {
      quadrant[xyz] = (pos[xyz] < center[xyz] ? 1 : 2);
   };
   return quadrant;
};





/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/28/2025
\param[in] level_idx Position in the "neighbor_levels" array
\return Relative refinement level of the neighbor block
*/
template <typename HConfig>
int Block<HConfig>::GetNeighborLevel(const MultiIndex& level_idx) const
{
   return neighbor_levels[level_idx.i + max_neighbors_per_dim_mns * (level_idx.j + max_neighbors_per_dim_mns * level_idx.k)];
};





/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] node_idx Position in the "neighbor_nodes" array
\return Node of the neighbor block
*/
template <typename HConfig>
int Block<HConfig>::GetNeighborNode(const MultiIndex& node_idx) const
{
   return neighbor_nodes[node_idx.i + max_neighbors_per_dim * (node_idx.j + max_neighbors_per_dim * node_idx.k)];
};





/*!
\author Vladimir Florinski
\date 06/15/2020
\param[in]  pos    Coordinates of a point
\param[out] zone   LN zone
\param[out] offset Displacement from the LN zone center in units of "zone_length"
*/
template <typename HConfig>
void Block<HConfig>::GetZoneOffset(const GeoVector& pos, MultiIndex& zone, GeoVector& offset) const
{
   offset = (pos - cent_min) / zone_length;
   zone = offset;

   for (auto xyz = 0; xyz < 3; xyz++) {
      if (offset[xyz] < 0.0) zone[xyz]--;
   };
   offset = offset - zone;
};


/*!
\author Vladimir Florinski
\date 04/28/2020
\param[in] pos Position
\return True if the vector is in the interior of the block
*/
template <typename HConfig>
bool Block<HConfig>::PositionInside(const GeoVector& pos) const
{
   return ((pos[0] >= face_min_phys[0]) && (pos[0] <= face_max_phys[0])
           && (pos[1] >= face_min_phys[1]) && (pos[1] <= face_max_phys[1])
           && (pos[2] >= face_min_phys[2]) && (pos[2] <= face_max_phys[2]));
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/21/2025
Call after setting face_min, face_max from the data reader.
*/
template <typename HConfig>
void Block<HConfig>::LoadDimensions(double unit_length_block)
{
// Scale units
   face_min *= unit_length_block / unit_length_fluid;
   face_max *= unit_length_block / unit_length_fluid;

// Copy to physical faces
   face_min_phys = face_min;
   face_max_phys = face_max;

// Adjust for ghost cells
   if (num_ghost_cells > 0) {
      GeoVector incr = (face_max - face_min);
      incr[0] *= ghost_to_phys_ratio[0];
      incr[1] *= ghost_to_phys_ratio[1];
      incr[2] *= ghost_to_phys_ratio[2];
      face_min -= incr;
      face_max += incr;
   };
   // todo what is the difference between this and configureproperties?
   //  why not call that here?
   // todo I am just going to add it myself...
   ConfigureProperties();
};



/*!
\author Vladimir Florinski
\date 06/23/2020
*/
template <typename HConfig>
void Block<HConfig>::ConfigureProperties(void)
{
// Zone sizes
   zone_length = (face_max - face_min) / block_size;

// FIXME
   if ((zone_length[0] <= 0.0) ||(zone_length[1] <= 0.0) || (zone_length[2] <= 0.0)) {
      std::cerr << "ConfigureProperties Error: Zone length for node " << node << " is negative.\n";
      std::cerr << face_max << "  " << face_min << "\n";
   };

   cent_min = face_min + 0.5 * zone_length;
   cent_max = face_max - 0.5 * zone_length;
   center = 0.5 * (face_min + face_max);
};



/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
*/
template <typename HConfig>
void Block<HConfig>::PrintFields(void) const requires (HConfig::build_type == Config::BuildMode::debug)
{
   int i, j, k, n;
   std::cerr << std::setprecision(8);
   for(k = 0; k < block_size.k; k++) {
      for(j = 0; j < block_size.j; j++) {
         for(i = 0; i < block_size.i; i++) {
            for(n = 0; n < DataFields::size(); n++)
               std::cerr << std::setw(16) << fields[fields_index({i,j,k})].Array()[n];
            std::cerr << std::endl;
         };
      };
   };
};



/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void Block<HConfig>::PrintNeighbors(void) const requires (HConfig::build_type == Config::BuildMode::debug)
{
   const int width = 10;
   const int shift = 3;

   std::cerr << std::endl;
   std::cerr << "Printing neighbors for block " << node << std::endl;
   std::cerr << std::endl;

// Print in an isometric perspective
   for(auto k = max_neighbors_per_dim_mns; k >= 0; k--) {
      for(auto j = max_neighbors_per_dim_mns; j >= 0; j--) {
         for(auto s = 0; s < j * shift; s++) std::cerr << " ";
         for(auto i = 0; i < max_neighbors_per_dim; i++) {
            std::cerr << std::setw(width) << GetNeighborNode(i,j,k);
         };
         std::cerr << std::endl;
      };
      std::cerr << std::endl;
   };
   std::cerr << std::endl;
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void Block<HConfig>::LoadCoords(void) requires (HConfig::build_type == Config::BuildMode::debug)
{
   for(auto k = 0; k < block_size.k; k++) {
      for(auto j = 0; j < block_size.j; j++) {
         for(auto i = 0; i < block_size.i; i++) {
            for(auto xyz = 0; xyz < 3; xyz++)
               fields[fields_index({i, j, k})].Array()[xyz] = cent_min[xyz] + i * zone_length[xyz];
         };
      };
   };
};


};
