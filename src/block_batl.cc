/*!
\file block_batl.cc
\brief Implements a class to operate on BATL blocks
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "block_batl.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockBATL methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
*/
BlockBATL::BlockBATL(void)
{
// Set these based on the constants in the header file
   n_variables = n_variables_batl;
   block_size = block_size_batl;
   max_neighbors_per_dim = max_neighbors_per_dim_batl;
   max_neighbors_per_dim_mns = max_neighbors_per_dim_batl - 1;
   max_neighbors = max_neighbors_batl;
   max_neighbor_levels = max_neighbor_levels_batl;
   AllocateMemory();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
\param[in] level_idx Position in the "neighbor_levels" array
\return Relative refinement level of the neighbor block
*/
int BlockBATL::GetNeighborLevel(const MultiIndex& level_idx) const
{
   return neighbor_levels[level_idx.i + max_neighbors_per_dim_mns * (level_idx.j + max_neighbors_per_dim_mns * level_idx.k)];
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
\param[in] i x-position in the "neighbor_levels" array
\param[in] j y-position in the "neighbor_levels" array
\param[in] k z-position in the "neighbor_levels" array
\return Relative refinement level of the neighbor block
*/
int BlockBATL::GetNeighborLevel(int i, int j, int k) const
{
   return neighbor_levels[i + max_neighbors_per_dim_mns * (j + max_neighbors_per_dim_mns * k)];
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024

\note Can be called from server processes only
*/
void BlockBATL::LoadDimensions(double unit_length_block)
{
   spectrum_get_block_corners(node, face_min.Data(), face_max.Data());
   BlockBase::LoadDimensions(unit_length_block);
};

/*!
\author Vladimir Florinski
\date 06/10/2020

\note Can be called from server processes only
*/
void BlockBATL::LoadNeighbors(void)
{
   spectrum_get_all_neighbor_copies(node, neighbor_nodes, neighbor_levels);
};

/*!
\author Vladimir Florinski
\date 06/08/2020

\note Can be called from server processes only
*/
void BlockBATL::LoadVariables(void)
{
   spectrum_get_block_data(node, variables);
};

};
