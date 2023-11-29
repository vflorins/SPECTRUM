/*!
\file block_cartesian.cc
\brief Implements a class to operate on Cartesian blocks
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "block_cartesian.hh"
#include "reader_cartesian.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockCartesian methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/30/2023
*/
void BlockCartesian::AllocateMemory(void)
{
// Undo the base class allocation, if any
   delete[] variables;
   delete[] neighbor_nodes;
   delete[] neighbor_levels;

// Allocate new storage
   variables = new double[n_variables * block_size.Prod()];
   neighbor_nodes = new int[max_neighbors];
   neighbor_levels = new int[max_neighbor_levels];
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
BlockCartesian::BlockCartesian(void)
{
// Set these based on the constants in the header file
   n_variables = n_variables_cartesian;
   block_size = block_size_cartesian;
   max_neighbors_per_dim = max_neighbors_per_dim_cartesian;
   max_neighbors_per_dim_mns = max_neighbors_per_dim_cartesian - 1;
   max_neighbors = max_neighbors_cartesian;
   max_neighbor_levels = max_neighbor_levels_cartesian;
   AllocateMemory();
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] pos Coordinates of a point
\return Which quadrant of the block the point is in (1 or 2)
*/
MultiIndex BlockCartesian::GetQuadrant(const GeoVector& pos) const
{
   MultiIndex quadrant;
   for(auto xyz = 0; xyz < 3; xyz++) {
      quadrant[xyz] = (pos[xyz] < center[xyz] ? 1 : 2);
   };
   return quadrant;
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] level_idx Position in the "neighbor_levels" array
\return Relative refinement level of the neighbor block
*/
int BlockCartesian::GetNeighborLevel(const MultiIndex& level_idx) const
{
   return 0;
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] i x-position in the "neighbor_levels" array
\param[in] j y-position in the "neighbor_levels" array
\param[in] k z-position in the "neighbor_levels" array
\return Relative refinement level of the neighbor block
*/
int BlockCartesian::GetNeighborLevel(int i, int j, int k) const
{
   return 0;
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] node_idx Position in the "neighbor_nodes" array
\return Node of the neighbor block
*/
int BlockCartesian::GetNeighborNode(const MultiIndex& node_idx) const
{
   return neighbor_nodes[node_idx.i + max_neighbors_per_dim * (node_idx.j + max_neighbors_per_dim * node_idx.k)];
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] i x-position in the "neighbor_nodes" array
\param[in] j y-position in the "neighbor_nodes" array
\param[in] k z-position in the "neighbor_nodes" array
\return Node of the neighbor block
*/
int BlockCartesian::GetNeighborNode(int i, int j, int k) const
{
   return neighbor_nodes[i + max_neighbors_per_dim * (j + max_neighbors_per_dim * k)];
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023

\note Can be called from server processes only
*/
void BlockCartesian::LoadDimensions(double unit_length_block)
{
   ReadCartesianGetBlockCorners(node, face_min.Data(), face_max.Data());
   BlockBase::LoadDimensions(unit_length_block);
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023

\note Can be called from server processes only
*/
void BlockCartesian::LoadNeighbors(void)
{
   ReadCartesianGetNodeNeighbors(node, neighbor_nodes, neighbor_levels);
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023

\note Can be called from server processes only
*/
void BlockCartesian::LoadVariables(void)
{
   ReadCartesianGetBlockData(node, variables);
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
*/
void BlockCartesian::PrintVariables(void) const
{
   int i, j, k, n;
   std::cerr << std::setprecision(8);
   for(k = 0; k < block_size.k; k++) {
      for(j = 0; j < block_size.j; j++) {
         for(i = 0; i < block_size.i; i++) {
            for(n = 0; n < n_variables; n++) std::cerr << std::setw(16) << variables[DataLoc(n,i,j,k)];
            std::cerr << std::endl;
         };
      };
   };
};

#ifdef GEO_DEBUG

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void BlockCartesian::PrintNeighbors(void)
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
void BlockCartesian::LoadCoords(void)
{
   for(auto k = 0; k < block_size.k; k++) {
      for(auto j = 0; j < block_size.j; j++) {
         for(auto i = 0; i < block_size.i; i++) {
            for(auto xyz = 0; xyz < 3; xyz++) variables[DataLoc(xyz, i, j, k)] = cent_min[xyz] + i * zone_length[xyz];
         };
      };
   };
};

#endif

};
