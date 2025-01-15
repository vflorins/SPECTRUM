/*!
\file block_base.cc
\brief Implements a simple class to operate on a 3D grid block
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "block_base.hh"
#include "common/physics.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/23/2020
\param[in] other Block to copy
*/
BlockBase::BlockBase(const BlockBase& other)
{
   operator =(other);
};

/*!
\author Vladimir Florinski
\date 01/12/2023
\param[in] other Block to move
*/
BlockBase::BlockBase(BlockBase&& other)
{
   operator =(other);
};

/*!
\author Vladimir Florinski
\date 10/27/2022
*/
BlockBase::~BlockBase(void)
{
// Although the base class does not allocate storage for the arrays, we can still delete them so that the derived classes don't have to
   delete[] neighbor_levels;
   delete[] neighbor_nodes;
   delete[] variables;
};

/*!
\author Vladimir Florinski
\date 10/26/2022
\param[in] other Block to copy
*/
BlockBase& BlockBase::operator =(const BlockBase& other)
{
   if (&other == this) return *this;

   node = other.node;
   face_min = other.face_min;
   face_max = other.face_max;
   ConfigureProperties();

   memcpy(neighbor_nodes, other.neighbor_nodes, max_neighbors * SZINT);
   memcpy(neighbor_levels, other.neighbor_levels, max_neighbor_levels * SZINT);
   memcpy(variables, other.variables, block_size.Prod() * n_variables * SZDBL);

   return *this;
};

/*!
\author Vladimir Florinski
\date 01/12/2023
\param[in] other Block to move
*/
BlockBase& BlockBase::operator =(BlockBase&& other)
{
   if (&other == this) return *this;

   node = other.node;
   face_min = other.face_min;
   face_max = other.face_max;
   ConfigureProperties();

   delete[] neighbor_levels;
   delete[] neighbor_nodes;
   delete[] variables;

   neighbor_levels = other.neighbor_levels;
   neighbor_nodes = other.neighbor_nodes;
   variables = other.variables;

   other.neighbor_levels = nullptr;
   other.neighbor_nodes = nullptr;
   other.variables = nullptr;

   return *this;
};

/*!
\author Juan G Alonso Guzman
\date 07/30/2023
\param[in] num_ghost_cells_in Number of ghost cells per side
*/
void BlockBase::SetGhostCells(int num_ghost_cells_in)
{
   num_ghost_cells = num_ghost_cells_in;

   ghost_to_phys_ratio[0] = (double)num_ghost_cells / (double)block_size.i;
   ghost_to_phys_ratio[1] = (double)num_ghost_cells / (double)block_size.j;
   ghost_to_phys_ratio[2] = (double)num_ghost_cells / (double)block_size.k;

// Update block_size to account for ghost cells if necessary
   if (num_ghost_cells > 0) {
      block_size = block_size + 2 * num_ghost_cells;

// Re-allocate memory for variables array
      delete[] variables;
      variables = new double[n_variables * block_size.Prod()];
   };
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023

\note Can be called from server processes only
*/
void BlockBase::LoadDimensions(double unit_length_block)
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
};

/*!
\author Vladimir Florinski
\date 06/23/2020
*/
void BlockBase::ConfigureProperties(void)
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
\author Vladimir Florinski
\date 08/13/2020
\param[in] pos Coordinates of a point
\return Index of the zone that owns "pos"
*/
MultiIndex BlockBase::GetZone(const GeoVector& pos) const
{
// This is the simplest implementation for cuboid blocks.
   MultiIndex zone = (pos - face_min) / zone_length;
   return zone;
};

/*!
\author Vladimir Florinski
\date 06/15/2020
\param[in]  pos    Coordinates of a point
\param[out] zone   LN zone
\param[out] offset Displacement from the LN zone center in units of "zone_length"
*/
void BlockBase::GetZoneOffset(const GeoVector& pos, MultiIndex& zone, GeoVector& offset) const
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
bool BlockBase::PositionInside(const GeoVector& pos) const
{
   return ((pos[0] >= face_min_phys[0]) && (pos[0] <= face_max_phys[0])
        && (pos[1] >= face_min_phys[1]) && (pos[1] <= face_max_phys[1])
        && (pos[2] >= face_min_phys[2]) && (pos[2] <= face_max_phys[2]));
};

/*!
\author Vladimir Florinski
\date 06/03/2020
\param[in] vidx Variable
\param[in] zone Zone's multi-index
\return Location in the array
*/
int BlockBase::DataLoc(int vidx, const MultiIndex& zone) const
{
   return vidx + n_variables * (zone.i + block_size.i * (zone.j + block_size.j * zone.k));
};

/*!
\author Vladimir Florinski
\date 10/26/2022
\param[in] vidx Variable
\param[in] i First index
\param[in] j Second index
\param[in] k Third index
\return Location in the array
*/
int BlockBase::DataLoc(int vidx, int i, int j, int k) const
{
   return vidx + n_variables * (i + block_size.i * (j + block_size.j * k));
};

/*!
\author Vladimir Florinski
\date 06/03/2020
\param[in] zone Zone's multi-index
\return True if the zone is in the block
*/
bool BlockBase::ZoneInsideBlock(const MultiIndex& zone) const
{
   if (   (zone.i < 0) || (zone.i >= block_size.i)
      || (zone.j < 0) || (zone.j >= block_size.j)
      || (zone.k < 0) || (zone.k >= block_size.k)) return false;
   else return true;
};

/*!
\author Vladimir Florinski
\date 06/19/2020
\param[in] zone Zone
\param[in] vidx Variable to retrieve
\return The value of the variable
*/
double BlockBase::GetValue(const MultiIndex& zone, int vidx) const
{
   return variables[DataLoc(vidx, zone)];
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\param[in] pos  Position
\param[in] vidx Variable to retrieve
\return The value of the variable
*/
double BlockBase::GetValue(const GeoVector& pos, int vidx) const
{
   MultiIndex zone = GetZone(pos);
   return variables[DataLoc(vidx, zone)];
};

};
